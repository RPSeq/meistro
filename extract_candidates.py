Skip to content
This repository
Search
Pull requests
Issues
Gist
 @RPSeq
 Unwatch 1
  Star 0
 Fork 0 RPSeq/meistro
 Code  Issues 0  Pull requests 0  Wiki  Pulse  Graphs  Settings
Branch: swalign Find file Copy pathmeistro/extract_candidates.py
29f5edd  a day ago
@RPSeq RPSeq Try using swalign package for polyA detection
1 contributor
RawBlameHistory     478 lines (356 sloc)  13.4 KB
#!/usr/bin/env python

import sys
from argparse import RawTextHelpFormatter, ArgumentParser
from itertools import izip
from string import maketrans

#installed modules
import pysam
from intervaltree import IntervalTree
import swalign

__author__ = "Ryan Smith (ryanpsmith@wustl.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2016-04-11 11:43 $"

# ====================
# SAM Class
# ====================
class sam_al(object):
    '''Class representing a SAM file alignment entry'''

    def __init__(self, sam, in_sam=False):
        #manual overloading based on arg types
        if type(sam) == pysam.AlignedSegment and in_sam:
            self.read_pysam(sam, in_sam)

        elif type(sam) == str or type(sam) == list:
            self.read(sam)

        else:
            exit("Error creating sam_al.\nUsage:sam_al(samlist), \
                sam_al(samstr), sam_al(pysam.al, pysam.infile)")

    def read(self, sam):
        if type(sam)==str:
            sam = sam.strip("\n").split("\t")
        self.qname = sam[0]
        self.flag = int(sam[1])
        self.rname = sam[2]
        self.pos = int(sam[3])
        self.mapq = int(sam[4])
        self.cigar = sam[5]
        self.rnext = sam[6]
        self.pnext = sam[7]
        self.tlen = int(sam[8])
        self.seq = sam[9]
        self.qual = sam[10]
        self.tags = {}

        for i in range(10,len(sam)):
            tag, ttype, val = sam[i].split(":")
            tags[tag]=(val)

        return

    def read_pysam(self, al, in_sam):
        self.qname = al.qname
        self.flag = al.flag
        self.rname = in_sam.getrname(al.tid)
        self.pos = al.pos
        self.mapq = al.mapq
        self.cigarstring = al.cigarstring
        self.cigar = al.cigar
        self.rnext = in_sam.getrname(al.mrnm)
        self.pnext = al.pnext
        self.tlen = al.tlen
        self.seq = al.seq
        self.qual = al.qual
        self.tags = dict(al.tags)
        self.is_secondary = al.is_secondary
        self.is_duplicate = al.is_duplicate
        self.is_proper_pair = al.is_proper_pair
        self.is_read1 = al.is_read1
        self.is_read2 = al.is_read2
        self.qstart = al.qstart
        self.qend = al.qend
        self.is_reverse = al.is_reverse
        self.is_paired = al.is_paired
        return

    def __str__(self, pair_tag=True):
        """Returns the sam record as a single string"""

        name = self.qname

        if pair_tag:

            if self.is_read1:
                name = self.qname+"_1"

            elif self.is_read2:
                name = self.qname+"_2"

        outlist = [name, str(self.flag), self.rname, 
                    str(self.pos), str(self.mapq), self.cigarstring, 
                    self.rnext, str(self.pnext), str(self.tlen), 
                    self.seq, self.qual]

        for tag, val in self.tags.viewitems():

            if type(val)==int:
                ttype = 'i'

            else:
                ttype = 'Z'

            outlist.append("{0}:{1}:{2}".format(tag, ttype, val))

        return "\t".join(outlist)+"\n"

#main loop function
def extract_candidates(bamfile, 
                        is_sam, 
                        anchors_out, 
                        fastq_out,
                        clip_len, 
                        single_only, 
                        max_opp_clip):
    # set input file
    if bamfile == None:

        if is_sam:
            in_bam = pysam.Samfile("-", "r")

        else:
            in_bam = pysam.Samfile('-', 'rb')
    else:

        if is_sam:
            in_bam = pysam.Samfile(bamfile, 'r')

        else:
            in_bam = pysam.Samfile(bamfile, "rb")

    header = "@HD\tVN:1.3\tSO:unsorted\n"   #write header to SAM output
    header+="\n".join(in_bam.text.split("\n")[1:]) 
    anchors_out = open(anchors_out, 'w')
    anchors_out.write(header)

    if fastq_out == "-":                #allow - output argument
        fastq_out = "/dev/stdout"

    fastq_out = open(fastq_out, 'w')    #open fastq output

    batchsize = 5
    anchor_batch = []
    fq_batch = []

    #iterate over the als
    for al in in_bam:
        anchor = False
        is_clip = False
        fastq_p, fastq_s = False, False

        #check if the batches need to be printed
        if len(anchor_batch) >= batchsize:
            anchors_out.write("".join(anchor_batch))
            del anchor_batch[:]

        if len(fq_batch) >= batchsize:
            fastq_out.write("".join(fq_batch))
            del fq_batch[:]

        #skip secondary or duplicate alignments
        if al.is_secondary or al.is_duplicate:
            continue

        #check if part of discordant pair.
        #should add zscore test for reads mapped to close/far together 
        #(and not use proper pair flag)
        conditions = [
            al.mapq == 0 and al.opt('MQ') > 0,
            al.is_unmapped != al.mate_is_unmapped,
            al.is_reverse == al.mate_is_reverse,
            not al.is_proper_pair,
            al.rname != al.mrnm
            ]

        if any(conditions):
            #use check pairs to determine which side to align (or both)
            remap, anchor = check_pairs(al, in_bam)

            if remap:
                fq_batch.append(remap)
            if anchor: 
                al = anchor

        al, is_clip = check_clip(al, in_bam, clip_len, max_opp_clip, anchor)

        if is_clip:
            fastq_s = fastq_str(al, is_clip)    # get fastq string
            fq_batch.append(fastq_s)            # append to fq output batch

            #pass to polyA function
            ssw_al = check_polyA(fastq_s.split("\n")[1])


            if ssw_al:  #if we got a polyA hit,
                al_tags = al.opt("TY").split(",")  #get the al's tags
                cigar, ori = ssw_al         

                #print cigar, ori       # get the polyA result

                if 'ASL' in al_tags:               #add the SR or SL -polyA tag.
                    pAtag = "SL"

                elif 'ASR' in al_tags:
                    pAtag = "SR"

                #generate the new tag and update al.
                newtag = pAtag+","+"polyA,0,"+cigar+","+ori
                al.setTag("RA",newtag)

        if anchor or is_clip:
            anchor_batch.append(str(sam_al(al, in_bam)))

    #after finishing, write the remaining items in the batches.
    anchors_out.write("".join(anchor_batch))
    fastq_out.write("".join(fq_batch))

    anchors_out.close()
    fastq_out.close()
            
# ============================================
# functions
# ============================================

def reverse_complement(sequence):
    """returns the reverse complement of the input DNA sequence"""

    complement = maketrans("ACTGactg", "TGACtgac")  #DNA translation table
    return sequence[::-1].translate(complement)     #reverse and translate

def hamming(str1, str2):
    """Returns the hamming distance between two strings of equal length"""
    assert len(str1) == len(str2)
    return sum(c1 != c2 for c1, c2 in izip(str1,str2))

def fastq_str(al, is_clip=False):
    """Returns a fastq string of the given BAM record"""
    seq = al.seq
    quals = al.qual
    name = al.qname

    tags=al.opt("TY") #get TY tag

    #if read is a clipper,
    if is_clip:
        #get tags (there will be at least one, "," delimited)
        tags = tags.split(",")

        for i in range(len(tags)):
            #if a tag is ASR or ASL, strip the A (fastq is NOT the anchor)
            #and pull the clipped portion of the read
            if tags[i] == "ASL":
                tags=tags[i][1:]
                seq = seq[:al.qstart]
                quals = quals[:al.qstart]

            elif tags[i] == "ASR":
                tags=tags[i][1:]
                seq = seq[al.qend:]
                quals = quals[al.qend:]

    #reverse the sequence if al is reversed
    if al.is_reverse:
        seq = reverse_complement(seq)
        quals = quals[::-1] 

    #add read # tags
    if al.is_read1:
        name += "_1"

    elif al.is_read2:
        name += "_2"

    #   return the fastq string, appending tags to read name (rname:tags)
    #   (BWA can add these FASTQ comments (rname<\t>comment) 
    #   to the output bam, but MOSAIK cannot.)
    return "@"+name+":"+tags+"\n"+seq+"\n+\n"+quals+"\n"

def check_pairs(al1, in_bam):
    """Determines if a read from a pair is a UU, UR, or RU."""

    anc, fq = False, False
    
    mate_mapq = al1.opt('MQ')   #get mate mapq

    if al1.mapq > 0 and mate_mapq > 0:  #if both are unique:
        al1.setTag("TY","UU")
        fq = fastq_str(al1)
        anc = al1               #realign both

    elif al1.mapq == 0 and mate_mapq > 0:   #if this al is not unique,
        al1.setTag("TY","RU")
        fq = fastq_str(al1)     #realign this al only
        anc = False
    
    elif al1.mapq > 0 and mate_mapq == 0:   #if other al is not unique,
        al1.setTag("TY","UR")
        fq = False
        anc = al1               #realign other al only

    return fq, anc

def check_clip(al, 
            in_bam, 
            clip_len, 
            max_opp_clip, 
            disc_pair):

    """Checks a given alignment for clipped bases on one end""" 
    cigar = al.cigar

    if (al.mapq == 0 or len(cigar)) < 2:
        return al, False

    try:
        tag = al.opt("TY")
    except:
        tag = False

    #if side=="L":
    #Clipped on L side
    if cigar[0][0] == 4 and cigar[0][1] >= clip_len:
        #if opposite is not clipped more than max opposite clip len, write for realignment
        if cigar[-1][0] != 4 or (cigar[-1][0] == 4 and cigar[-1][1] <= max_opp_clip):

            if tag:
                al.setTag("TY",tag+",ASL")

            else:
                al.setTag("TY","ASL")

            return al, True

    #elif side=="R":
    #Clipped on R side
    elif cigar[-1][0] == 4 and cigar[-1][1] >= clip_len:
    #if opposite is not clipped more than max opposite clip len, write for realignment
        if cigar[0][0] != 4 or (cigar[0][0] == 4 and cigar[0][1] <= max_opp_clip):

            if tag:
                al.setTag("TY",tag+",ASR")

            else:
                al.setTag("TY","ASR")

            return al, True

    return al, False

def check_polyA(q_seq):

    match = 2
    mismatch = -1
    q_len = len(q_seq)
    max_score = float(match*q_len)
    polyA_ref = "A"*q_len

    matrix = swalign.NucleotideScoringMatrix(match, mismatch)
    #default gap penalty=-1, gap_extend=-1, extend_decay=0.0
    sw = swalign.LocalAlignment(matrix)

    #get best forward/reverse al
    best = None
    for strand in '+-':
        if strand == '-':
            aln = sw.align(polyA_ref, swalign.revcomp(q_seq), rc=True)
        else:
            aln = sw.align(polyA_ref, q_seq)

        if not best or aln.score > best.score:
            best = aln

    strand = "+"
    if best.rc:
        strand = "-"

    if best.score/max_score >= 0.8:
        return best.cigar_str, strand

    return False

def zscore(val, mean, stdev):
    return abs((float(val)-float(mean))/float(stdev))

def filter_excludes(variants, exclude_file):
    """Returns variants that do not intersect features in a given BED file"""
    #excludes is a dict of IntervalTrees (one for each chromosome)
    excludes = defaultdict(IntervalTree)
    filtered = [] #list to be populated with filtered items

    for line in exclude_file:
        if line[0] != "#":  #skip headers
            line = line.strip().split("\t")
            chrom = line[0]
            start = int(line[1])
            stop = int(line[2])

            #add chrom:interval to excludes
            #intervaltree.addi(start, stop, data)
            excludes[chrom].addi(start, stop, 1)

    #could probably speed this up using a mapping function instead
    for variant in variants:
        #get interval tree for var chrom, and query with the position.
        if len(excludes[variant.chrom][variant.pos]) == 0:
            filtered.append(variant)

    return filtered

def get_args():
    parser = ArgumentParser(
        formatter_class=RawTextHelpFormatter, add_help=False)

    parser.add_argument('-a', metavar='SAM_OUT', required=True, 
                        help='Output anchors SAM (required)')

    parser.add_argument('-f', metavar='FASTQ_OUT', required=True,
                        help='Output FASTQ for realignment (required)')

    parser.add_argument('-i', metavar='BAM', 
                        help='Input BAM (stdin)')

    parser.add_argument('-c', metavar='MIN_CLIP', type=int, default=25,
                        help='Minimum clip length (25)')

    parser.add_argument('-oc', metavar='MAX_OPP_CLIP', type=int, default=7, 
                        help='Maximum opposite clip length (7)')

    parser.add_argument('-s', action='store_true', help='Input single-ended')

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

    extract_candidates(args.i, args.S, args.a, 
                        args.f, args.c, args.s, args.oc)

if __name__ == "__main__":
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE error for streaming
            raise




Status API Training Shop Blog About
© 2016 GitHub, Inc. Terms Privacy Security Contact Help