import pysam
import sys
import argparse
from argparse import RawTextHelpFormatter


__author__ = "Ryan Smith (ryanpsmith@wustl.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2016-2-8 13:45 $"

def cluster(bamfile, is_sam, out_file="-"):
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

    out_file = pysam.Samfile("-", 'wh', template=in_bam)

    #assuming name-sorted bam
    filtered = []
    clust = anchor_cluster()

    #sys.stderr.write("Filtering...\n")
    for al in in_bam:
        anch = anchor(al, in_bam)
        if not clust.add(anch):
            if clust.get_breakends():
                filtered.append(clust)
            clust = anchor_cluster(anch)

    for clust in filtered:
        for item in clust.anchors:
            out_file.write(item.al)



class anchor_cluster(object):
    '''Class representing a cluster of MEI anchors'''

    def __init__(self, anchor=False):
        self.ID = False         #cluster ID
        self.anchors = []    #collector for anchors in cluster
        
        self.SL = []
        self.SR = []

        self.polyA_L = []
        self.polyA_R = []

        self.RU_L = []
        self.RU_R = []

        self.BND_L = 0
        self.BND_R = float("inf")

        if anchor:
            self.add(anchor)

    def add(self, anchor):
        if not self.ID:
            self.ID = anchor.ID
        elif self.ID:
            if anchor.ID != self.ID:
                return False

        self.anchors.append(anchor)

        if anchor.SR:
            self.SR.append(anchor)
        if anchor.SL:
            self.SL.append(anchor)

        if anchor.polyA_R:
            self.polyA_R.append(anchor)
        if anchor.polyA_L:
            self.polyA_L.append(anchor)

        if anchor.UR:
            tag, mei, pos, cigar, ori = anchor.UR
            if ori == "+":
                self.RU_L.append(anchor)
            elif ori == "-":
                self.RU_R.append(anchor)

        return True

    def get_breakends(self, min_num=0):
        #if we don't actually have any MEI reads (lots of polyA-only pileups to be discarded)
        if sum([len(self.SL), len(self.SR), len(self.RU_L), len(self.RU_R)]) <= 2:
            return False

        #if len(self.polyA_L) > 0 and len(self.polyA_R) > 0:
         #   return False

        if self.RU_R:
            for anchor in self.RU_R:
                if anchor.start < self.BND_R:
                    self.BND_R = anchor.start
        if self.RU_L:
            for anchor in self.RU_L:
                if anchor.end > self.BND_L:
                    self.BND_L = anchor.end
        if self.SR:
            for anchor in self.SR:
                if anchor.end > self.BND_L:
                    self.BND_L = anchor.end
        if self.SL:
            for anchor in self.SL:
                if anchor.start < self.BND_R:
                    self.BND_R = anchor.start
        if self.polyA_L:
            for anchor in self.polyA_L:
                if anchor.start < self.BND_R:
                    self.BND_R = anchor.start
        if self.polyA_R:
            for anchor in self.polyA_R:
                if anchor.end > self.BND_L:
                    self.BND_L = anchor.end
        return True

class anchor(object):

    def __init__(self, al, al_bam):
        self.al = al
        self.bam = al_bam
        self.chrom = self.bam.getrname(al.rname)
        self.start = int(al.pos)
        self.end = int(al.aend)

        self.SR = False
        self.SL = False
        self.UR = False

        self.polyA_L = False
        self.polyA_R = False

        self.ori = "+"
        if al.is_reverse:
            self.ori = "-"

        try:
            self.tags = [x.split(",") for x in al.opt("RA").split(";")]
        except:
            sys.stderr.write("Error: MEI anchor reads must have RA tag added by filter_merged.py\n")
            exit(1)

        try:
            self.ID = int(al.opt("CL"))
        except:
            sys.stderr.write("Error: MEI anchor reads must have CL tag added via bedtools cluster\n(currently in meistro.sh pipeline)\n")
            exit(1)

        self.parse_tags()

    def parse_tags(self):
        #tag = [(SL/SR/UR/UU),mei_name,start,cigar,ori]
        for tag in self.tags:
            if tag[0][0] == 'S':    #if its a split-read tag
                if tag[0][1] == 'L':
                    if tag[1] == 'polyA': #check if polyA or mei
                        self.polyA_L = tag
                    else:
                        self.SL = tag
                elif tag[0][1] == 'R':
                    if tag[1] == 'polyA':
                        self.polyA_R = tag
                    else:
                        self.SR = tag

            elif tag[0] == "UR" or tag[0] == "UU":
                self.UR = tag


def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
cluster.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: Process cluster.bam in MEIstro pipeline")
    parser.add_argument('-i', '--input', metavar='BAM', required=False, help='Input BAM file')
    parser.add_argument('-S', required=False, action='store_true', help='Input is SAM format')
    parser.add_argument('-o', required=False, help='Output ? file')

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
    cluster(args.input, args.S, args.o)

if __name__ == "__main__":
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise