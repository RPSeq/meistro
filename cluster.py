import pysam
import sys
import argparse
from argparse import RawTextHelpFormatter
from collections import defaultdict


__author__ = "Ryan Smith (ryanpsmith@wustl.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2016-2-8 13:45 $"

class cluster(object):
    #subclasses for the 4 different sub-clusters in a MEI read cluster.
    #doesn't have to be this way
    def __init__(self):
        self.als = []
        self.chrom = False
        self.start = False
        self.end = False

    def add_al(self, al):
        if True:
            self.al.add(read)

class split_clust(object):
    def __init__(self):
        self.als = []
        self.ori = False
        self.chrom = False
        self.start = False
        self.end = False

    def add_al(self,al):

class pair_clust(object):
    def __init__(self):
        self.als = []
        self.ori = False
        self.start = False
        self.end = False

    def add_al(self, al):

def scan(bamfile, pA_file, is_sam, out_file="-"):
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

    out_file = open(out_file, "w")

    pA_file = pysam.Samfile(pA_file, 'rb')
    #fetch like this:#
    ##################
    # for pa in pA_file.fetch("3", 0,  115670485):
    #     count+=1
    # print(str(count))
    ###################

    #assuming name-sorted bam

    #need to have a "current cluster" for each MEI family. 
    #1- how should I define families? Easy way: first two chars of MEI refname.
    #with the Mobster ref, this allows AL, L1, SV, and HE
    #we also need to have Sl, SR, PL, PR, and account for the polyA signal.
    for al in in_bam:
        #tag = [(SL/SR/UR/UU),mei_name,start,cigar,ori]
        chrom = in_bam.getrname(al.rname)
        try:
            tags = [x.split(",") for x in al.opt("RA").split(";")]
        except:
            sys.stderr.write("Error: MEI anchor reads must have RA tag added by filter_merged.py\n")
            exit(1)
        pos = al.


def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
cluster.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: Process cluster.bam in MEIstro pipeline")
    parser.add_argument('-i', '--input', metavar='BAM', required=False, help='Input filtered BAM file')
    parser.add_argument('-pa', '--polyA', metavar='BAM', required=True, help='Input polyA BAM file')
    parser.add_argument('-S', required=False, action='store_true', help='Input is SAM format')
    parser.add_argument('-o', required=False, help='Output file')

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
# driver
# ============================================

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def main():
    args = get_args()
    cluster(args.input, args.polyA, args.S, args.o)

if __name__ == "__main__":
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise