import pysam
import sys
from argparse import RawTextHelpFormatter, ArgumentParser
from collections import defaultdict


__author__ = "Ryan Smith (ryanpsmith@wustl.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2016-2-8 13:45 $"


# here i define a new convention: + paired anchors are R, as the MEI segment is on the RIGHT from the anchor segment.
# and thus, - paired anchors are L, as the MEI segment is on the LEFT relative to the anchor.
# now, as long as the orientations for a cluster are concordant, (+-,++) are PR and (-+, --) are PL.
class Cluster(object):
    def __init__(self, cluster):
        #mei_family, 
        #   anchor_type,
        #       orientation : anchor_list
        self.collector = defaultdict(lambda: 
                            defaultdict(lambda: 
                                defaultdict(list)))

        self.ev_count = 0
        self.polyA_count = 0
        self.mei_count = 0

        #map concordant split oris to same value to simplify collection
        self.split_conc_oris = {"+-":"+--+", "-+":"+--+", 
                                "++":"++--","--":"++--"}

        self.pair_conc = {"+-":"-+", "-+":"+-", 
                        "++":"--", "--":"++"}

        #collect the anchors into the hash
        self.collect(cluster)

        #group them into subclusters
        self.load_sub_clusters()


    def collect(self, anchors):
        for anch in anchors:

            self.ev_count += 1
            if anch.me == "polyA":
                self.polyA_count += 1

            else:
                self.mei_count += 1

            ori = anch.ori

            if anch.read_type == "S":
                ori = self.split_conc_oris[ori]

            #for now, the mei_refname's first two chars are the mei_family hash keys
            self.collector[anch.me_fam][anch.RA_tag][ori].append(anch)


    def load_sub_clusters(self):
        #traverse the collector, generating subclusters as we go
        #need to load this info into better data structure now
        for me_fam, sides in self.collector.iteritems():

            for RA_tag, oris in sides.iteritems():

                for ori, anchs in oris.iteritems():
                    self.collector[me_fam][RA_tag][ori] = self.gen_subclusters(anchs, RA_tag)


    def gen_subclusters(self, clust, RA_tag):
        """Generate subclusters from a list of anchors"""

        #sort anchors by end pos if MEI on right.
        #(incoming reads are sorted by start pos)
        if RA_tag[1] == "R":
            #sort anchors by end position
            clust.sort(key=lambda x: (x.end, x.ori))

        clusters = []
        prev = None

        for anch in clust:

            if not prev:
                prev = anch
                curr_clust = SubCluster(prev, self)
                continue

            curr = anch

            if not curr_clust.add(prev):
                clusters.append(curr_clust)
                curr_clust = SubCluster(curr, self)

            prev = curr

        clusters.append(curr_clust)

        return clusters


    def filter(self):
        """Define filter for clusters"""
        if self.mei_count > 10:

            return True

        return False


class DoubleCluster(object):

    def __init__(self, SubClusterA, SubClusterB, dist):

        side = SubClusterA.side
        self.read_type == SubClusterA.read_type

        if SubClusterA.side == SubClusterB.side:
            sys.stderr.write("MergedCluster Error: \
                SubClusters to merge have same side (L/R)")
            exit(1)

        if SubClusterA.read_type != SubClusterB.read_type:
            sys.stderr.write("MergedCluster Error: \
                SubClusters to merge have different read type (S/R)")
            exit(1)

        if self.read_type == "P":
            if SubClusterA.cluster.pair_conc[SubClusterA.ori] == SubClusterB.ori:
                return False

        if self.read_type == "S":
            if SubClusterA.ori != SubClusterB.ori:
                return False

        if side == "L":
            self.L_Cluster = SubClusterA
            self.R_Cluster = SubClusterB

        elif side == "R":
            self.R_Cluster = SubClusterA
            self.L_Cluster = SubClusterB

        self.dist = dist
        self.clust_type == self.read_type

    def add(self, clust):
    # def get_tsd(self)

class SubCluster(object):

    def __init__(self, anchor, cluster):
        self.cluster = cluster
        self.anchors = [anchor]
        self.side = anchor.side
        self.me_fam = anchor.me_fam
        self.RA_tag = anchor.RA_tag
        self.read_type = anchor.read_type
        self.side = anchor.side

        ori = anchor.ori
        if self.read_type == "S":
            ori = self.cluster.split_conc_oris[ori]

        if self.side == "L":
            pos = anchor.start
        elif self.side == "R":
            pos = anchor.end

        self.pos = pos
        self.start = pos
        self.end = pos

        if self.read_type == "S":
            self.CLUST_DIST = SPLIT_MERGE_DIST
            self.MERGE_DIST = SPLIT_MERGE_OLAP
            self.MERGE_OLAP = SPLIT_CLUST_DIST

        elif self.read_type == "P":
            self.CLUST_DIST = PAIR_MERGE_DIST
            self.MERGE_DIST = PAIR_MERGE_OLAP
            self.MERGE_OLAP = PAIR_CLUST_DIST


    def add(self, anchor):
        if self.RA_tag != anchor.RA_tag or self.me_fam != anchor.me_fam:
            return False

        if anchor.clust_pos() - self.end > self.CLUST_DIST:
            return False

        else:
            self.anchors.append(anchor)
            self.end = anchor.clust_pos()
            return True

    def clust_pos(self):
        """Returns inner clustering position"""

        if self.side == "R":   return self.end

        elif self.side == "L": return self.start


    def distance(self, subcluster):
        """Returns the distance between two SubClusters (- dist means overlap)"""

        if self.side == "R":
            return subcluster.clust_pos() - self.clust_pos()

        elif self.side == "L":
            return self.clust_pos() - subcluster.clust_pos()

    def ori_concordant(self, subcluster):

        if self.read_type == "S" and self.ori == subcluster.ori:
            return True

        elif self.read_type == "P" and self.cluster.pair_conc[subcluster.ori] == self.ori:
            return True


    def merge(self, subcluster):
        """Returns False, False if discordant, 
        False, dist if too far/excess olap,
        and True, DoubleCluster object if can be merged."""

        #not concordant for merging
        if self.side == cluster.side or \
        self.read_type != cluster.read_type or \
        self.me_fam != cluster.me_fam or \
        not self.ori_concordant(cluster):
            return False, False

        #get distance between inside ends
        dist = self.distance(cluster)

        #too much overlap or too far:
        if dist < (0-self.MERGE_OLAP) or dist > self.MERGE_DIST:
            return False, dist

        else:
            return True, DoubleCluster(self, cluster, dist)


class Anchor(object):
    """Encapsulates an mei realignment anchor"""

    def __init__(self, al, al_bam):
        self.al = al
        self.bam = al_bam
        self.name = al.qname
        self.chrom = self.bam.getrname(al.rname)
        self.start = int(al.pos)
        self.end = int(al.aend)
        self.cigar = al.cigarstring
        self.mapq = al.mapq

        #read mei tag for an anchor
        try:
            self.tag_str = al.opt("RA")
        except:
            sys.stderr.write("cluster.py Error: Reads must have RA tag added from filter_merged.py")
            exit(1)

        if len(self.tag_str.split(";")) > 1:
            sys.stderr.write("cluster.py Error: Multiple RA tags found for anchor")
            exit(1)

        tag = self.tag_str.split(",")

        self.RA_tag = tag[0]
        self.read_type = tag[0][0]
        self.side = tag[0][1]
        self.me = tag[1]
        self.me_fam = tag[1][:2]
        self.me_start = int(tag[2])
        self.me_cigar = tag[3]
        self.ori = tag[4]

    def clust_pos(self):
        if self.side == "L":
            return self.start

        elif self.side == "R":
            return self.end


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
        # if clust.filter(): 
            # sys.stdout.write(clust.__str__(count))
            # sys.stdout.write("\n")
        #if clust.filter():  clust.merge_pair_clusters()

def cluster_generator(bamfile, max_dist):
    """Generator function that clusters \
    bam entries and yields a list for each cluster."""

    prev = Anchor(bamfile.next(), bamfile)  #grab first alignment
    cluster = [prev]    #initialize first cluster

    for al in bamfile:

        curr = Anchor(al, bamfile)
        #if the cluster range is exceeded:
        if (curr.start - prev.end > max_dist) or curr.chrom != prev.chrom:

            yield cluster       #yield current clust
            cluster = [curr]    #initialize next clust

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

    parser.add_argument('-S', action='store_true', 
                        help='Input is SAM format')  

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

    #   PARAMETERS: THIS WAY IS REALLY BAD. NEED TO REWORK PARAMETER SETTING!!
    global PRECLUSTER_DISTANCE
    global PAIR_MERGE_DIST
    global PAIR_MERGE_OLAP
    global SPLIT_MERGE_DIST
    global SPLIT_MERGE_OLAP
    global PAIR_CLUST_DIST
    global SPLIT_CLUST_DIST

    #   distance for bedtools cluster-like generic clustering
    PRECLUSTER_DISTANCE = 1800

    #   how close do reads need to be to be in same sub-cluster?
    #   measured from rightmost-base if upstream of MEI (PR/SR)
    #   leftmost-base if downstream of MEI (PL/SL)
    PAIR_CLUST_DIST = 400
    SPLIT_CLUST_DIST = 8

    #   how close do concordant sub-clusters need to be to be merged?
    #   what is the max overlap for merging pair clusters? 
    #   (this effectively sets the max allowed TSDup / TSDel sizes)
    #   NOTE: TSDel is far far less common than no TS alteration,
    #   and TSDup is by FAR the most common)
    PAIR_MERGE_DIST = 800
    PAIR_MERGE_OLAP = 60

    SPLIT_MERGE_DIST = 20
    SPLIT_MERGE_OLAP = 60

    scan(args.i, args.S)

if __name__ == "__main__":
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise