import pysam
import sys
from argparse import RawTextHelpFormatter, ArgumentParser
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

class Cluster(object):
    def __init__(self, cluster=False):
        self.anchor_types = ["PR","PL","SR","SL"]
        self.meis = ["Al", "L1", "SV", "HE"]
        self.hash = defaultdict(dict)
        self.polyA_hash = {"SR": [], "SL": []}

        #create hash structure
        for mei in self.meis:
            for side in self.anchor_types:
                self.hash[mei][side] = {"+-": [], "-+": [],
                                        "++": [], "--": []}

        for side in self.polyA_hash:
            self.polyA_hash[side] = {"+-": [], "-+": [],
                                    "++": [], "--": []}
        #{
        #'SV': {'PR': [], 'SR': [], 'PL': [], 'Sl': []}, 
        #'HE': {'PR': [], 'SR': [], 'PL': [], 'Sl': []}, 
        #'AL': {'PR': [], 'SR': [], 'PL': [], 'Sl': []}, 
        #'L1': {'PR': [], 'SR': [], 'PL': [], 'Sl': []}
        #}

        #load the anchors into the hash
        if cluster: 
            self.load(cluster)
            self.generate_sub_clusters()

    def load(self, cluster):
        for anch in cluster:
            for tag in anch.tags:
                tagstr = tag.RA_type
                #separate MEIs als from polyA als.
                if tag.mei != "polyA":
                    self.hash[tag.mei[:2]][tagstr][tag.ori].append(anch)
                else:
                    self.polyA_hash[tagstr][tag.ori].append(anch)

    def generate_sub_clusters(self):
        for mei_fam, sides in self.hash.iteritems():
            for side, oris in sides.iteritems():
                for ori, anchs in oris.iteritems():
                    self.hash[mei_fam][side][ori] = self.sub_cluster(anchs, side)

        for side, oris in self.polyA_hash.iteritems():
            for ori, anchs in oris.iteritems():
                self.polyA_hash[side][ori] = self.sub_cluster(anchs, side)

    def filter(self):
        """Define filter for clusters"""
        mei_hit = 0
        polyA_hit = 0
        for mei_fam, sides in self.hash.iteritems():
            for side, oris in sides.iteritems():
                for ori, clusts in oris.iteritems():
                    if clusts:
                        for clust in clusts:
                            mei_hit +=1

        for side, clusters in self.polyA_hash.iteritems():
            for clust in clusters:
                polyA_hit +=1

        if mei_hit > 5:
            return True

        else:
            return False
=======
                    self.hash[tag.mei[:2]][tagstr].append(anch)
                else:
                    self.polyA_hash[tagstr].append(anch)

    def generate_sub_clusters(self):
        for mei_fam, sides in self.hash.iteritems():
            for side, anchs in sides.iteritems():
                self.hash[mei_fam][side] = self.sub_cluster(anchs, side)

        for side, anchs in self.polyA_hash.iteritems():
            self.polyA_hash[side] = self.sub_cluster(anchs, side)
>>>>>>> bbe89d81f9a1b31c4a9831d0b13f1f959f563fde

    def sub_cluster(self, anch_list, side):

        if not anch_list: return False

        #if MEI on right side: 
        if side[1] == "R":
            #re-sort anchors by end position
            anch_list.sort(key=lambda x: x.end)

        if side[0] == "S":
            clusters = self.clust(anch_list, SPLIT_CLUST_DIST, side[1])

        elif side[0] == "P":
            clusters = self.clust(anch_list, PAIR_CLUST_DIST, side[1])

        return clusters

<<<<<<< HEAD
=======

>>>>>>> bbe89d81f9a1b31c4a9831d0b13f1f959f563fde
    def clust(self, clust, distance, side):
        clusters = []
        prev = None

        for anch in clust:
            if not prev:
                prev = anch
                curr_clust = [prev]
                continue

            curr = anch

            if not self.check_dist(curr, prev, distance, side):
                clusters.append(curr_clust)
                curr_clust = [curr]

            else:
                curr_clust.append(curr)

            prev = curr

        clusters.append(curr_clust)

        return clusters

    def check_dist(self, curr, prev, distance, side):
        if side == "L":
            if curr.start - prev.start > distance:
                return False

        if side == "R":
            if curr.end - prev.end > distance:
                return False

        return True


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
#####################################################################
#####################################################################

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

    count = 0
    for group in cluster_generator(in_bam, PRECLUSTER_DISTANCE):
        count += 1
        clust = Cluster(group)
<<<<<<< HEAD
        if clust.filter():
            myhash = clust.hash
            for mei_fam, sides in myhash.iteritems():
                for side, oris in sides.iteritems():
                    for ori, clusters in oris.iteritems():
                        if clusters:
                            for clust in clusters:
                                if clust:
                                    for anch in clust:
                                        print(anch.chrom, anch.start, anch.end, ori, mei_fam, side)
=======
        print(clust.hash)
>>>>>>> bbe89d81f9a1b31c4a9831d0b13f1f959f563fde

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
    parser = ArgumentParser(
        formatter_class=RawTextHelpFormatter, add_help=False)

    parser.add_argument('-i', metavar='BAM', 
                        help='Input BAM (stdin)')

    parser.add_argument('-S', action='store_true', help='Input is SAM format')  

    parser.add_argument('-h','--help', action='help') 


    # parse the arguments
    args = parser.parse_args()
    
    # bail if no BAM file
    if args.i is None:
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
    global PRECLUSTER_DISTANCE
    global PAIR_CLUST_DIST
    global SPLIT_CLUST_DIST

    PRECLUSTER_DISTANCE = 1300
    PAIR_CLUST_DIST = 300
    SPLIT_CLUST_DIST = 7

    # PAIR_MERGE_DIST = 
    # PAIR_MERGE_OLAP =

    # SPLIT_MERGE_DIST = 
    # SPLIT_MERGE_OLAP =

    scan(args.i, args.S)

if __name__ == "__main__":
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise