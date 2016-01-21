#!/usr/bin/env python

# for tgi cluster:
#/gapp/x64linux/opt/pythonbrew/venvs/Python-2.7.6/gemini/bin/python
# for uva cluster:

import pysam
import sys
import argparse
from argparse import RawTextHelpFormatter
import string
from string import *
from collections import defaultdict

__author__ = "Ryan Smith (ryanpsmith@wustl.edu) with code by Colby Chiang (cc2qe@virginia.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2014-12-15 11:43 $"


# ============================================
# functions
# ============================================

def filter_merged(bamfile, is_sam, out_bam, mei_names):
    """Main BAM parsing loop"""

    # set input file
    if bamfile == None:
        if is_sam:
            in_bam = pysam.Samfile('-', 'r')
        else:
            in_bam = pysam.Samfile('-', 'rb')
    else:
        if is_sam:
            in_bam = pysam.Samfile(bamfile, 'r')
        else:
            in_bam = pysam.Samfile(bamfile, 'rb')

    out_bam = pysam.Samfile(out_bam, 'wh', template=in_bam)

    # header = "@HD\tVN:1.3\tSO:queryname\n"
    # header+="\n".join(in_bam.text.split("\n")[1:])

    #get meinames from file
    meis = set(line.strip() for line in mei_names)


    #assuming name-sorted bam
    filtered = []
    group = Namegroup(in_bam, meis)



    #sys.stderr.write("Filtering...\n")
    for al in in_bam:
        if al.rname < 0:
            continue
        #if qname doesnt match current group:
        if not group.add(al):
            #pass the group to filtering function
            hit = group.process()
            #if the group hit an MEI, append to filtered list.
            if hit:
                filtered.extend(group.anchors)
            #create a new namegroup
            group = Namegroup(in_bam, meis, al)

    for al in filtered:
        out_bam.write(al)

    #sys.stderr.write("Done!\n")

def mismatch(al):
    """Returns number of mismatches in given alignment"""

    cigar = al.cigar
    #these flags are all mismatches (I,D,N,S,H,X)
    mismatches = set([1,2,3,4,5,8])
    NM = sum([sgmt[1] for sgmt in cigar if sgmt[0] in mismatches])
    length = len(al.seq)
    percent = NM/(float(length))
    return percent
        
def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
bamgroupreads.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: Group BAM file by read IDs without sorting")
    parser.add_argument('-i', '--input', metavar='BAM', required=False, help='Input BAM file')
    parser.add_argument('-S', required=False, action='store_true', help='Input is SAM format')
    parser.add_argument('-o', required=True, help='Output MEI BAM')
    parser.add_argument('-m', '--mei_names', required=True, type=argparse.FileType('r'), help='List of MEI reference names')
    # parse the arguments
    args = parser.parse_args()
    
    # bail if no BAM file
    if args.input is None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)
    
    # send back the user input
    return args


# ============================================
# classes
# ============================================
class Namegroup(object):
    """Class for containing reads of the same name"""

    def __init__(self, in_bam, meis, al=False):

        self.bam = in_bam
        self.als = []
        self.anchors = []
        self.mei = False
        self.mei_rnames = meis

        if al:
            self.name = al.qname.split("_")[0]
            self.als.append(al)

        else:
            self.name = False

    def add(self,al):
        name, num = al.qname.split("_")
        tags = al.opt("TY").split(",")
        #if self.name hasn't been set, set with current al.
        if self.name == False:
            self.name = name
        #if the name doesn't match the group, return false
        if name != self.name:
            return False

        self.als.append(al)
        return True

    def process(self):
        """Returns true if group hit an MEI"""

        for al in self.als:
            #get name and read #
            name, num = al.qname.split("_")
            tags = al.opt("TY").split(",")
            #some untig als have non-spec rname indices
            if al.rname < 0:
                continue
            if not al.is_secondary: #ignore secondary hits for now
                #check for both UUs aligned to MEI
                if 'UU' in tags:
                    #if its an MEI alignment, and another doesn't exist,
                    #keep it.
                    if self.bam.getrname(al.rname) in self.mei_rnames:
                        if not self.mei:
                            self.mei = al
                        #if we've already seen an mei UU alignment,
                        #both sides of the UU pair have mapped to the MEI.
                        else:
                            return False
                    else:
                        self.anchors.append(al)

                elif any(x in tags for x in ['SL','SR','RU']):
                    self.mei = al
                #check for anchor aligned tags
                elif any(x in tags for x in ['ASL','ASR','UR']):
                    self.anchors.append(al)

        # check for anchor ref chrom discrepancies 
        # if there are multiple anchor reads in the group.
        #     discard if there are mismatches.
        if len(self.anchors) > 1:
            rname = self.anchors[0].rname
            for anch in self.anchors[1:]:
                if anch.rname != rname:
                    return False

        #some (known) cases left to account for:
        # UU/UU-Split.
        #       First UU isn't informative, only the split realignment matters.
        #       Should probably remove the UU-only anchor (leaving only the split anchor)
        # UU-anchor/UU-mei:
        #       This one is probably fine for now, may want to relabel though.
        # RU/UR-split:

        #if we have anchor(s) and an mei, time to add the mei info to the anchors.
        if len(self.anchors) > 0 and self.mei:
            #if we have mei and anchor, return true.
            for anc in self.anchors:

                tags = anc.opt("TY").split(",")
                mei_name = self.bam.getrname(self.mei.rname)
                mei_pos = self.mei.pos
                mei_cigar = self.mei.cigarstring
                mei_qual = self.mei.mapq
                mei_tags = ";".join(self.mei.opt("TY").split(","))
                if self.mei.is_reverse:
                    mei_ori = "-"
                else:
                    mei_ori = "+"
                anc.setTag("ME", ",".join([mei_name,self.mei.qname,mei_tags,str(mei_pos),mei_ori,mei_cigar,str(mei_qual)]))
                #anc.rnext = self.mei.rname
                #anc.mpos = self.mei.pos
                anc.tlen = 0
            return True

        return False

# ============================================
# driver
# ============================================

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def main():
    args = get_args()
    filter_merged(args.input, args.S, args.o, args.mei_names)

if __name__ == "__main__":
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise

