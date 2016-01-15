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

def filter_merged(bamfile, is_sam, out_bam):
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

    #assuming name-sorted bam
    filtered = []
    group = Namegroup(in_bam)

    sys.stderr.write("Filtering...\n")
    for al in in_bam:
        #if qname doesnt match current group:
        if not group.add(al):
            #pass the group to filtering function
            hit = group.process()
            #if the group hit an MEI, append to filtered list.
            if hit:
                filtered.extend(group.als)
            #create a new namegroup
            group = Namegroup(in_bam, al)

    for al in filtered:
        out_bam.write(al)

    sys.stderr.write("Done!\n")

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

    def __init__(self, in_bam, al=False):

        self.bam = in_bam
        self.als = []

        if al:
            self.name = al.qname.split("_")[0]
            self.als.append(al)

        else:
            self.name = False

    def add(self,al):
        if self.name == False:
            self.name = al.qname.split("_")[0]
        elif al.qname.split("_")[0] != self.name:
            return False

        self.als.append(al)
        return True

    def process(self):
        """Returns true if group hit an MEI"""

        #Tags used for characterizing an MEI namegroup.
        #If an alignment in the group matches a tag, the tag is set to that alignment.
        #MUU and AUU are MEI Unique-Unique and Anchor Unique-Unique
        MUU1,MUU2,AUU1,AUU2=False,False,False,False
        #UR and RU are Repeat-Unique and Unique-Repeat
        UR1,UR2,RU1,RU2=False,False,False,False
        #ASL and ASR are (Anchor, Split Left) and (Anchor, Split Right)
        ASL1,ASL2,ASR1,ASR2=False,False,False,False
        #SL and SR are Split Left and Split Right
        SL1,SL2,SR1,SR2=False,False,False,False


        for al in self.als:
            #get name and read #
            name, num = al.qname.split("_")
            tags = al.opt("TY").split(",")
            if al.rname < 0:
                continue
            if not al.is_secondary: #ignore secondary hits for now

                # Here we update all the tag info
                if "UU" in tags:
                    if "moblist" in self.bam.getrname(al.rname):
                        if num == "1":
                            MUU1 = al
                        elif num == "2":
                            MUU2 = al
                    else:
                        if num == "1":
                            AUU1 = al
                        elif num == "2":
                            AUU2 = al

                if "UR" in tags:
                    if num == "1":
                        UR1 = al
                    elif num == "2":
                        UR2 = al

                if "RU" in tags:
                    if num == "1":
                        RU1 = al
                    elif num == "2":
                        RU2 = al

                if "ASL" in tags:
                    if num == "1":
                        ASL1 = al
                    elif num == "2":
                        ASL2 = al

                if "ASR" in tags:
                    if num == "1":
                        ASR1 = al
                    elif num == "2":
                        ASR2 = al

                if "SL" in tags:
                    if num == "1":
                        SL1 = al
                    elif num == "2":
                        SL2 == al

                if "SR" in tags:
                    if num == "1":
                        SR1 = al
                    elif num == "2":
                        SR2 = al


        # ==================================== #
        # main logic for filtering mei groups! #
        # ==================================== #

        #If both sides of a UU pair map to MEI, remove.
        if MUU1 and MUU2:
            #UNLESS 
            # if (SR1 and SL2) or (SR2 and SL1):
            #     return True
            return False

        #If UR pair, keep
        if (UR1 and RU2) or (UR2 and RU1):
            return True

        #split left, keep
        if (SL1 and ASL1) or (SL2 and ASL2):
            return True

        #split right, keep
        if (SR1 and ASR1) or (SR2 and ASR2):
            return True




        

# ============================================
# driver
# ============================================

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def main():
    args = get_args()
    filter_merged(args.input, args.S, args.o)

if __name__ == "__main__":
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise

