import pysam
import sys
import argparse
from argparse import RawTextHelpFormatter
from collections import defaultdict


__author__ = "Ryan Smith (ryanpsmith@wustl.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2016-2-8 13:45 $"


#we will have one active cluster for each family (as we see the reads)
#scanning across. als are collected into their appropriate cluster.
#if evidence of an MEI, query for polyA signal in the appropriate region/ori
#when cluster is finished collecting, need to merge subclusters (can be PR,PL,SR,SL (and polyA?))
class cluster(object):
    def __init__(self, mei):
        self.PR = []
        self.PL = []
        self.SR = []
        self.SL = []
        self.pA = []
        self.PR_end = 0
        self.PL_end = 0
        self.SR_end = 0
        self.SL_end = 0
        self.pA_end = 0
        self.ori = False
        self.mei = mei

class split_cluster(object):
    def __init__(self):
        self.anchors = []
        self.side = False
        self.end = 0
        self.ori = ori

class pair_cluster(object):
    def __init__(self):
        self.anchors = []
        self.side = False
        self.end = 0

class mei_tag(object):
    '''Encapsulates an mei realignment tag'''
    def __init__(self, tag):
        tag = tag.split(",")
        self.type = tag[0]
        self.mei = tag[1]
        self.start = tag[2]
        self.cigar = tag[3]
        self.ori = tag[4]

class anchor(object):
    '''Encapsulates an mei realignment anchor'''
    def __init__(self, al, al_bam):
        self.al = al
        self.bam = al_bam
        self.name = al.qname
        self.chrom = self.bam.getrname(al.rname)
        self.start = int(al.pos)
        self.end = int(al.aend)
        self.cigar = al.cigarstring
        self.ori = "+"
        if al.is_reverse:
            self.ori = "-"
        #load mei_tags for an anchor
        try:
            self.tags = [mei_tag(tag) for tag in al.opt("RA").split(";")]
        except:
            sys.stderr.write("cluster.py Error: Reads must have RA tags added from filter_merged.py\n")
            exit(1)

def scan(bamfile, pA_file, is_sam, out_file="/dev/stdout"):
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

    #out_file = open(out_file, "w")

    pA_file = pysam.Samfile(pA_file, 'rb')
    #fetch like this:#
    ##################
    # for pa in pA_file.fetch("3", 0,  115670485):
    #     count+=1
    # print(str(count))
    ###################

    ######class testing#####
    # in_bam.next()
    # for i in range(1000):
    #     #tag = [(SL/SR/UR/UU),mei_name,start,cigar,ori]
    #     al = in_bam.next()
    #     anc = anchor(al, in_bam)
    #     for tag in anc.tags:
    #         sys.stdout.write("\t".join([anc.name, anc.chrom, str(anc.start),str(anc.end), anc.ori, tag.mei, tag.type, tag.ori])+"\n")
    ######class testing######

    for al in in_bam:
        anc = anchor(al, in_bam)
        for tag in anc.tags:
            if tag.type = "SL":
                #something

            elif tag.type = "SR":
                #something

            elif tag.type = "UR":
                #something

            elif tag.type = "UU":
                #something

        


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
    scan(args.input, args.polyA, args.S, args.o)

if __name__ == "__main__":
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise