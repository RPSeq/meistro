#!/usr/bin/env python

import pysam
import sys
import argparse
from argparse import RawTextHelpFormatter
import string
from string import *
from collections import defaultdict

__author__ = "Ryan Smith (ryanpsmith@wustl.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2016-04-13 10:52 $"

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

    #get meinames from file (#should probably get this from difference between the headers of the anchors and merged bams instead)
    mei_names = set(line.strip() for line in mei_names)

    #assuming name-sorted bam
    filtered = []
    group = Namegroup(in_bam, mei_names)

    #sys.stderr.write("Filtering...\n")
    for al in in_bam:
        if al.rname < 0:
            continue
        #if qname doesnt match current group:
        if not group.add(al):
            #pass the group to filtering function
            hits = group.process()
            #if the group hit an MEI, append to filtered list.
            if hits:
                for hit in hits:
                    out_bam.write(hit)
            #create a new namegroup
            group = Namegroup(in_bam, mei_names, al)

    #sys.stderr.write("Done!\n")
  
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

    def __init__(self, in_bam, mei_names, al=False):

        self.bam = in_bam
        self.als = []
        self.anchors = []
        self.meis = []
        self.mei_rnames = mei_names
        self.uu_mei = [False, False]
        self.filtered_anchors = []

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

    def update_RA(self, RAtag, mei, anchor):

        ori = "+"
        if mei.is_reverse:
            ori = "-"

        a_ori = "+"
        if anchor.is_reverse:
            a_ori = "-"

        try:
            # DON"T WORRY: anchors can have 2 items in their TY tags,
            # but any single MEI realignment can only have one TY tag. 
            # (if an anchor has two TY tags, it could have two mei reals, each with one TY tag.)
            RA_type = mei.opt("TY")
        except:
            sys.stderr.write("Reads must have TY tags; make sure the input .bam was processed by extract_candidates.py")

        #set pair direction tag based on anchor ori. now these will match the split tags.
        if RA_type == "UU" or RA_type == "RU":
            RA_type = "PR"
            if anchor.is_reverse:
                RA_type = "PL"

        if not RAtag:
            RAtag = RA_type+","+",".join([self.bam.getrname(mei.rname), 
                                                str(mei.pos), 
                                                mei.cigarstring, 
                                                a_ori+ori])
        else:
            #when we use the Mobster Mosaik moblist ref, L1HS has a polyA tail starting at 6017 bp.
            # this seems to "soak up" a lot of actual polyA sequences, so lets just ignore it for now
            # and see if the SSW check_polyA function works properly.
            if "polyA" in RAtag:
                if self.bam.getrname(mei.rname) == "L1HS" and mei.pos >= 6000:
                    return RAtag
                    
            RAtag += ";" + RA_type+","+",".join([self.bam.getrname(mei.rname), str(mei.pos), mei.cigarstring, a_ori+ori])

        return RAtag

    def process(self):
        """Returns true if group hit an MEI"""

        for al in self.als:
            #get name and read pair number
            name, num = al.qname.split("_")
            num = int(num)
            #get tags
            tags = al.opt("TY").split(",")
            #some untig als have non-spec rname indices
            if al.rname < 0:
                continue
            if not al.is_secondary: #ignore secondary hits for now
                if self.bam.getrname(al.rname) in self.mei_rnames:
                    self.meis.append(al)
                    if tags[0] == "UU":
                        self.uu_mei[num-1] = True
                else:
                    self.anchors.append(al)

        #if both UUs aligned to mei, discard
        if all(self.uu_mei):
            return False

        for anchor in self.anchors:
            mnum = False
            uu_hit = False
            try:
                #account for the polyA RA tag that already exists if this readgroup aliged to polyA.
                RAtag = anchor.opt("RA")
            except:
                RAtag = False

            anum = int(anchor.qname.split("_")[1])
            #if anchor is a UU only
            tags = anchor.opt("TY").split(",")

            #lots of code duplication in this next section.
            # also, i should delete an mei from self.meis when it's added to the RAtag,
            # this could speed things up a bit, cutting down on the number of iterations through self.meis
            if "UU" in tags:
                #look for UU realignment in meis
                for mei in self.meis:
                    if mei.opt("TY") == "UU":
                        mnum = int(mei.qname.split("_")[1])
                        #if its the mate of the UU, report it and add the new tag.
                        if mnum != anum:
                            uu_hit = True
                            RAtag = self.update_RA(RAtag, mei, anchor)
                            break

                if not uu_hit:
                    #delete the UU tag if there are other tags, no longer relevant
                    if len(tags) > 1:
                        #this only works because read pairs are ALWAYS checked before clippers
                        #in extract_candidates. thus, if a read is a split and a UU, the UU tag will come first.
                        #SHIT NOW, polyA hits of UUs MAY be deleted! what do I do about this? (think i fixed it.)
                        #del tags[0]
                        for i in range(len(tags)):
                            if tags[i] == "UU":
                                del tags[i]
                                break                        

            if "UR" in tags:
                for mei in self.meis:
                    if mei.opt("TY") == "RU":
                        mnum = int(mei.qname.split("_")[1])
                        if mnum != anum:
                            RAtag = self.update_RA(RAtag, mei, anchor)
                            break

            if 'ASL' in tags:
                hit = False
                for mei in self.meis:
                    if mei.opt("TY") == 'SL':
                        mnum = int(mei.qname.split("_")[1])
                        if mnum == anum:
                            RAtag = self.update_RA(RAtag, mei, anchor)
                            break

            elif 'ASR' in tags:
                for mei in self.meis:
                    if mei.opt("TY") == 'SR':
                        mnum = int(mei.qname.split("_")[1])
                        if mnum == anum:
                            RAtag = self.update_RA(RAtag, mei, anchor)
                            break

            if RAtag:
                #don't add to output if ONLY polyA signal (these will be in a separate file)
                # if len(RAtag.split(";")) == 1 and RAtag.split(",")[1] == "polyA":
                #     break
                anchor.setTag("RA", RAtag)
                self.filtered_anchors.append(anchor)

        if len(self.filtered_anchors) > 0:
            return self.filtered_anchors

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

