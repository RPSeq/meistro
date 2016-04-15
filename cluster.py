import pysam
import sys
import argparse
from argparse import RawTextHelpFormatter
from collections import defaultdict

__author__ = "Ryan Smith (ryanpsmith@wustl.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2016-2-8 13:45 $"


# we will have one active cluster for each family (as we see the reads)
# scanning across. als are collected into their appropriate sub_cluster.
# when cluster is finished collecting, need to merge subclusters (can be PR,PL,SR,SL (and polyA?))

# clust type can be pair or split. but, i have a generic L or R orientation for splits
# and do not have the same for pairs; however there is no reason i can't.
# so, here i define a new convention: + paired anchors are R, as the MEI segment is on the RIGHT from the anchor segment.
# and thus, - paired anchors are L, as the MEI segment is on the LEFT relative to the anchor.
# now, as long as the orientations for a cluster are concordant, (+-,++) are PR and (-+, --) are PL.
# it might be a good idea to modify filter_merged.py to change RU and UU to PR and PL.
# BUT: SR and SL do not have a native anchor orientation; that is, a PL can still be in the positive orientation,
# but the 5' side of the segment is clipped and remapped to an MEI. how do i account for this in the generic sub_cluster type?


from collections import defaultdict
hash = defaultdict(dict)

sides = ["PR","PL","SR","Sl"]
meis = ["AL", "L1", "SV", "HE"]
for mei in meis:
    for side in sides:
        hash[mei][side] = []
#{
#'SV': {'PR': [], 'SR': [], 'PL': [], 'Sl': []}, 
#'HE': {'PR': [], 'SR': [], 'PL': [], 'Sl': []}, 
#'AL': {'PR': [], 'SR': [], 'PL': [], 'Sl': []}, 
#'L1': {'PR': [], 'SR': [], 'PL': [], 'Sl': []}
#}

class ReadCluster(object):
    def __init__(self, cluster):

class Cluster(object):
    def __init__(self, cluster):
        self.mei_hash = {"PR":[], "SR":[], "PL":[], "SL":[]}
        self.polyA_hash = {"SR":[], "SL":[]}

        #collect the reads in their appropriate group
        for anch in cluster:
            for tag in anch.tags:
                #separate MEIs als from polyA als.
                if tag.mei != "polyA":
                    self.mei_hash[tag.RA_type].append(anchor)
                else:
                    self.polyA_hash[tag.RA_type].append(anchor)

    def process(self):

        for RA_type, anchors in self.mei_hash:
            for anch in anchors:



class meiTag(object):
    """Encapsulates an mei realignment tag"""
    def __init__(self, tag):
        tag = tag.split(",")
        self.RA_type = tag[0]
        self.mei = tag[1]
        self.start = tag[2]
        self.cigar = tag[3]
        self.ori = tag[4]

class anchor(object):
    """Encapsulates an mei realignment anchor"""
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
            self.tag_str = al.opt("RA")
            self.tags = [meiTag(tag) for tag in self.tag_str.split(";")]
        except:
            sys.stderr.write("cluster.py Error: Reads must have RA tags added from filter_merged.py\n")
            exit(1)


def scan(bamfile, is_sam):
    """Main BAM parsing loop"""
    #max_dist should be somthing like insert size + 2 stdev
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

    for group in cluster_generator(in_bam, PRECLUSTER_DISTANCE):
        cluster = ME_cluster(group)

def cluster_generator(bamfile, max_dist):
    """Generator function that clusters bam entries and yields a list for each cluster."""

    prev = anchor(bamfile.next(), bamfile)  #grab first alignment
    cluster = [prev]    #initialize first cluster

    for al in bamfile:

        curr = anchor(al, bamfile)
        #if the cluster range is exceeded, yield current clust and initialize next cluster
        if (curr.start - prev.end > max_dist) or curr.chrom != prev.chrom:

            yield cluster
            cluster = [curr]

        else:
            #otherwise append to cluster
            cluster.append(curr)

        prev = curr

    yield cluster



def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
cluster.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: Process cluster.bam in MEIstro pipeline")
    parser.add_argument('-i', '--input', metavar='BAM', required=False, help='Input filtered BAM file')
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

    #set global clustering parameters
    PRECLUSTER_DISTANCE = 1300
    # PAIR_CLUST_DIST = 
    # PAIR_MERGE_DIST = 
    # PAIR_MERGE_OLAP =

    # SPLIT_CLUST_DIST = 
    # SPLIT_MERGE_DIST = 
    # SPLIT_MERGE_OLAP =

    scan(args.input, args.S)

if __name__ == "__main__":
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise