#!/usr/bin/env python

# for tgi cluster:
#/gapp/x64linux/opt/pythonbrew/venvs/Python-2.7.6/gemini/bin/python
# for uva cluster:

import pysam
import sys
import argparse
from itertools import izip
from intervaltree import IntervalTree
from argparse import RawTextHelpFormatter
from string import maketrans
from ssw_wrap import Aligner

__author__ = "Ryan Smith (ryanpsmith@wustl.edu) with code by Colby Chiang (cc2qe@virginia.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2014-12-15 11:43 $"

# ====================
# SAM Class
# ====================
class sam_al(object):
    '''Class representing a SAM file alignment entry'''

    def __init__(self, sam, in_sam=False):
        #manual overloading based on arg types
        if type(sam) == pysam.AlignedRead and in_sam:
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

    def sam_str(self, tag=False):
        """Returns the sam record as a single string"""
        name = self.qname
        if tag:
            if self.is_read1:
                name = self.qname+"_1"
            elif self.is_read2:
                name = self.qname+"_2"
        outlist = [name, str(self.flag), self.rname, 
                    str(self.pos), str(self.mapq), self.cigarstring, 
                    self.rnext, str(self.pnext), str(self.tlen), self.seq, self.qual]

        for tag, val in self.tags.viewitems():
            if type(val)==int:
                ttype = 'i'
            else:
                ttype = 'Z'
            outlist.append("{0}:{1}:{2}".format(tag, ttype, val))
        return "\t".join(outlist)+"\n"

    def write(self, output):
        line = self.sam_str()
        output.write(line)
        return

#main loop function
def extract_candidates(bamfile, 
                        is_sam, 
                        anchors_out, 
                        fastq_out,
                        clip_len, 
                        single_only, 
                        max_opp_clip=7):
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

    #generate header for anchors output file
    header = "@HD\tVN:1.3\tSO:unsorted\n"
    header+="\n".join(in_bam.text.split("\n")[1:])
    
    #open and print header for anchor and polyA bams
    anchors_out = open(anchors_out, 'w')
    anchors_out.write(header)


    #allow - input argument
    if fastq_out == "-":
        fastq_out = "/dev/stdout"
    fastq_out = open(fastq_out, 'w')



    #these collect output strings,
    #to be written to output file in batches to reduce IO
    anchor_batch = []
    fq_batch = []


    #create striped smith-waterman aligner object
    #calibrated so gaps are not allowed, only mismatches.
    polyA_ssw = Aligner("A"*clip_len,
                match=4,
                mismatch=8,
                gap_open=900,
                gap_extend=600,
                report_secondary=False,
                report_cigar=True)

    #this sets the batch size
    batchsize = 1000000

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
        if (al.rname != al.mrnm) or (al.is_reverse == al.mate_is_reverse) or \
            (al.is_unmapped != al.mate_is_unmapped) or \
            (not al.is_proper_pair) or (al.mapq == 0 and al.opt('MQ') > 0):

            #use check pairs to determine which side to align (or both)
            fastq_p, anchor = check_pairs(al, in_bam)

            if fastq_p:
                fq_batch.append(fastq_p)
            if anchor: 
                al = anchor

        al, is_clip = check_clip(al, in_bam, clip_len, max_opp_clip, anchor)

        if is_clip:
            fastq_s = fastq_str(al, is_clip)    # get fastq string
            fq_batch.append(fastq_s)            # append to fq output batch
            seq = fastq_s.split("\n")[1]        # get just the seq line
            ssw_al = check_polyA(seq, polyA_ssw) #pass to polyA function

            #if we got a polyA hit,
            if ssw_al:
                al_tags = al.opt("TY").split(",")   #get the al's tags
                cigar, ori = ssw_al                 # get the polyA result
                #need to add a new tag
                if 'ASL' in al_tags:                #add the SR or SL -polyA tag.
                    pAtag = "SL"
                elif 'ASR' in al_tags:
                    pAtag = "SR"

                #generate the new tag and update al.
                newtag = pAtag+","+"polyA,0,"+cigar+","+ori
                al.setTag("PA",newtag)

        if anchor or is_clip:
            anchor_batch.append(sam_al(al, in_bam).sam_str(1))


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

    complement = maketrans("ACTGactg", "TGACtgac") #define translation table for DNA
    return sequence[::-1].translate(complement)  #return the reversed and translated sequence

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
                # tags[i]=tags[i][1:]
                tags=tags[i][1:]
                seq = seq[:al.qstart]
                quals = quals[:al.qstart]

            elif tags[i] == "ASR":
                # tags[i]=tags[i][1:]
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


    #return the fastq string, appending tags to read name (rname:tags)
    #NOTE: typical fastq format is rname<\t>comment
    #(BWA can add these FASTQ comments to the output bam as tags, but other aligners cannot.)
    return "@"+name+":"+tags+"\n"+seq+"\n+\n"+quals+"\n"

def check_pairs(al1, in_bam):
    """Determines if a read from a pair is a UU, UR, or RU."""
    #default return values are False
    anc, fq = False, False
    #if this read has been marked as a clipper,

    #get mate mapq
    mate_mapq = al1.opt('MQ')

    #if both are unique:
    if al1.mapq > 0 and mate_mapq > 0:
        al1.setTag("TY","UU")
        fq = fastq_str(al1)
        anc = al1

    #if this al is not unique,
    elif al1.mapq == 0 and mate_mapq > 0:
        #realign this al
        al1.setTag("TY","RU")
        anc = False
        fq = fastq_str(al1)
    
    #if other al is not unique,
    elif al1.mapq > 0 and mate_mapq == 0:
        #realign other al, not this one (unless clipped).
        al1.setTag("TY","UR")
        fq = False
        anc = al1

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

    if disc_pair:
        tag = al.opt("TY")

    #if side=="L":
    #Clipped on L side
    if cigar[0][0] == 4 and cigar[0][1] >= clip_len:
        #if opposite is not clipped more than max opposite clip len, write for realignment
        if cigar[-1][0] != 4 or (cigar[-1][0] == 4 and cigar[-1][1] <= max_opp_clip):
            if disc_pair:
                al.setTag("TY",tag+",ASL")
            else:
                al.setTag("TY","ASL")
            return al, True

    #elif side=="R":
    #Clipped on R side
    elif cigar[-1][0] == 4 and cigar[-1][1] >= clip_len:
    #if opposite is not clipped more than max opposite clip len, write for realignment
        if cigar[0][0] != 4 or (cigar[0][0] == 4 and cigar[0][1] <= max_opp_clip):
            if disc_pair:
                al.setTag("TY",tag+",ASR")
            else:
                al.setTag("TY","ASR")
            return al, True

    return al, False

def check_polyA(seq, 
            polyA_ssw):

    hit_f = polyA_ssw.align(seq, min_score=10, min_len=12)
    hit_r = polyA_ssw.align(reverse_complement(seq), min_score=10,min_len=12)
    
    cutoff = 0.8
    f_percent = 0
    r_percent = 0
    if hit_f:
        f_len = (hit_f.query_end-hit_f.query_begin)+1
        f_percent = hit_f.score/float(f_len*4)

    if hit_r:
        r_len = (hit_r.query_end-hit_r.query_begin)+1
        r_percent = hit_r.score/float(r_len*4)

    if (f_percent >= cutoff) and (r_percent >= cutoff):
        if f_percent > r_percent:
            return hit_f.cigar, "+"
        else:
            return hit_r.cigar, "-"
    elif f_percent >= cutoff:
        return hit_f.cigar, "+"
    elif r_percent >= cutoff:
        return hit_r.cigar, "-"
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
            excludes[chrom].addi(start, stop, 1) #intervaltree.addi(start, stop, data)

    #could probably speed this up using a mapping function instead
    for variant in variants:
        #get interval tree for var chrom, and query with the position.
        if len(excludes[variant.chrom][variant.pos]) == 0:
            filtered.append(variant)
    return filtered

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
bamgroupreads.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: Extract candidates for MEI re-alignment")
    parser.add_argument('-i', '--input', metavar='BAM', required=False, help='Input BAM file [stdin]')
    parser.add_argument('-a', '--anchors', metavar='SAM', required=True, help='Output anchors SAM')
    parser.add_argument('-f', '--fastq', metavar='FASTQ', required=True, help='Output realign FASTQ')
    parser.add_argument('-c', '--clip', metavar='LEN', required=True, type=int, help='Minimum clip length')
    parser.add_argument('-oc', '--opclip', metavar='LEN', required=False, type=int, help='Max opposite clip length')
    parser.add_argument('-s', '--single', required=False, action='store_true', help='Input single-ended')
    parser.add_argument('-S', required=False, action='store_true', help='Input is SAM format')    


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

    extract_candidates(args.input, args.S, args.anchors, args.fastq, args.clip, args.single, args.opclip)

if __name__ == "__main__":
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise

