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

        #collect the anchors into the hash
        self.collect(cluster)

        #group them into subclusters
        self.load_sub_clusters()

    def collect(self, cluster):
        for anch in cluster:

            for tag in anch.tags:
                self.ev_count += 1

                if tag.mei == "polyA":
                    self.polyA_count += 1

                else:
                    self.mei_count += 1

                ori = tag.ori

                if tag.RA_type[0] == "S":
                    ori = self.split_conc_oris[ori]

                #for now, the mei_refname's first two chars are the mei_family hash keys
                self.collector[tag.mei[:2]][tag.RA_type][ori].append(anch)

    def load_sub_clusters(self):
        #traverse the collector, generating subclusters as we go
        #need to load this info into better data structure now
        for mei_fam, sides in self.collector.iteritems():

            for side, oris in sides.iteritems():

                for ori, anchs in oris.iteritems():
                    self.collector[mei_fam][side][ori] = self.gen_subclusters(anchs, side, ori)



    def gen_subclusters(self, clust, side, ori):
        """Generate subclusters from a list of anchors"""

        #sort anchors by end pos if MEI on right.
        if side[1] == "R":
            #sort anchors by end position
            clust.sort(key=lambda x: x.end)

        clusters = []
        prev = None

        for anch in clust:

            if not prev:
                prev = anch
                curr_clust = SubCluster(prev, side, ori)
                continue

            curr = anch

            if not curr_clust.add(prev, side):
                clusters.append(curr_clust)
                curr_clust = SubCluster(curr, side, ori)

            prev = curr

        clusters.append(curr_clust)

        return clusters


    def filter(self):
        """Define filter for clusters"""
        if self.mei_count > 10:

            return True

        return False

class MergedCluster(object):
    def __init__(self, ClusterA, ClusterB, dist):

        if ClusterA.dir == ClusterB.dir or ClusterA.dir:
            sys.stderr.write("Error: Clusters to merge have same side (L/R")
            exit(1)

        side = ClusterA.side
        if side == "L":
            self.L_Cluster = Cluster.A
            self.R_Cluster = Cluster.B

        elif side == "R":
            self.R_Cluster = Cluster.A
            self.L_Cluster = Cluster.B

        self.r_type == ClusterA.r_type
        self.dist = dist

    # def get_tsd(self)

class SubCluster(object):
    def __init__(self, Cluster, anchor, side, ori):
        self.cluster = Cluster
        self.anchors = [anchor]
        self.side = side
        self.read_type = side[0]
        self.dir = side[1]
        self.ori = ori

        pos = anchor.start
        if self.dir == "R":
            pos = anchor.end

        self.pos = pos
        self.start = pos
        self.end = pos

        if self.r_type == "S":
            self.clust_dist = SPLIT_MERGE_DIST
            self.merge_dist = SPLIT_MERGE_OLAP
            self.merge_olap = SPLIT_CLUST_DIST

        elif self.r_type == "P":
            self.clust_dist = PAIR_MERGE_DIST
            self.merge_dist = PAIR_MERGE_OLAP
            self.merge_olap = PAIR_CLUST_DIS

    def add(self, anchor, side):
        if self.side != side:
            return False

        a_pos = self.get_a_pos(anchor)

        if a_pos - self.end > self.clust_dist:
            return False

        else:
            self.anchors.append(anchor)
            self.end = a_pos
            return True

    def get_a_pos(self, anchor):

        a_pos = anchor.end
        if self.side == "L":
            a_pos = anchor.start

        return a_pos

    # def get_ori(self):
    #     anch = self.anchors[0]
    #     for tag in anch.tags:
    #         if tag.RA_type == self.side:
    #             return tag.ori

    def clust_pos(self):
        """Returns cluster coord from appropriate R or L side"""

        if self.dir == "R":   return self.end

        elif self.dir == "L": return self.start

    def clust_dist(self, cluster):
        """Returns the distance between two SubClusters"""

        if self.dir == "R":
            return cluster.clust_pos() - self.clust_pos

        elif self.dir == "L":
            return self.clust_pos - cluster.clust_pos()

    def merge(self, cluster):
        """Returns False, False if discordant, 
        False, dist if too far/excess olap,
        and True, MergedCluster object if can be merged."""

        #not concordant for merging
        if self.dir == cluster.dir or self.read_type != cluster.read_type:
            return False, False

        #get distance between inside ends
        dist = self.clust_dist(cluster)

        #get split/pair cutoff values
        m_olap, m_dist = PAIR_MERGE_OLAP, PAIR_MERGE_DIST
        if self.read_type == "P":
            m_olap, m_dist = SPLIT_MERGE_OLAP, SPLIT_MERGE_DIST

        #too much overlap or too far:
        if dist < (0-m_olap) or dist > m_dist:
            return False, dist

        else:
            return True, MergedCluster(self, cluster, dist)





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
        self.mapq = al.mapq

        if al.is_reverse:
            self.ori = "-"
        #load mei_tags for an anchor
        try:
            self.tag_str = al.opt("RA")
            self.tags = [meiTag(tag) for tag in self.tag_str.split(";")]
        except:
            sys.stderr.write("cluster.py Error: Reads must have RA \
                                tags added from filter_merged.py\n")
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
        # if clust.filter(): 
            # sys.stdout.write(clust.__str__(count))
            # sys.stdout.write("\n")
        #if clust.filter():  clust.merge_pair_clusters()

def cluster_generator(bamfile, max_dist):
    """Generator function that clusters \
    bam entries and yields a list for each cluster."""

    prev = anchor(bamfile.next(), bamfile)  #grab first alignment
    cluster = [prev]    #initialize first cluster

    for al in bamfile:

        curr = anchor(al, bamfile)
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
    SPLIT_CLUST_DIST = 7

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