#!/bin/sh

#Check if correct command line input
if [ $# -ne 4 ]; then
    echo "Usage: bash pipeline.sh <results prefix> <BAM file path> <MEI reference path> <output directory name>"
    exit 1
fi

#Get command line input
PREFIX=$1
INPUT_BAM=$2
MEI_REF=$3
OUTPUT_DIR_NAME=$4
MOBSTER=$5


#Set relative dirs
WORKING_DIR=$(readlink -f $PWD)
SCRIPTS_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
FILE_LABEL=$(basename $INPUT_BAM ".bam")
FILE_DIR=$(dirname $INPUT_BAM)

#Mosaik Files
IA="/gscmnt/gc2719/halllab/users/rsmith/transposons/Mobster/Mobster-0.1.6/mobiome/54_mobiles_inclHERVK"
ANN="/gscmnt/gc2719/halllab/users/rsmith/transposons/Mobster/Mobster-0.1.6/MOSAIK"
REPMASK="/gscmnt/gc2719/halllab/users/rsmith/git/meistro/repmask/nchr.SLOP90.repmask_hg19_Alu.L1.SVA.ERV.bed"

##Set output directory
RESULTS_DIR=${WORKING_DIR}/${OUTPUT_DIR_NAME}
OUTPUT=$RESULTS_DIR/$PREFIX

##Check if output dir exists and ask to overwrite
if [ -d ${RESULTS_DIR} ]; then
    echo "Output directory exists: overwriting."
else
    mkdir ${RESULTS_DIR}
fi

###############
### BWA MEM ###
###############

# sambamba view -t 4 -f bam -l 0 $INPUT_BAM \
#     | python ${SCRIPTS_DIR}/extract_candidates.py -a >(samtools view -b - > ${OUTPUT}.anchors.bam) -f - -c 25 -oc 7 \
#         | bwa mem -t 8 -k 11 -M -C $MEI_REF /dev/stdin \
#             | samtools view -F 4 -b -u - | samtools sort -n - ${OUTPUT}.realigned

# samtools merge - ${OUTPUT}.anchors.bam ${OUTPUT}.realigned.bam | samtools sort -n - ${OUTPUT}.merged

# python ./filter_merged.py -i ${OUTPUT}.merged.bam -o /dev/stdout | samtools view -b - > ${OUTPUT}.filtered.bam


###############
### MOSAIK ###
###############

#extract candidate reads for realingment to MEI library
sambamba view -t 4 -f bam -l 0 $INPUT_BAM | python ${SCRIPTS_DIR}/extract_candidates.py -a >(samtools view -b - > ${OUTPUT}.anchors.bam) -f - -c 20 -oc 7 -R mosaik > ${OUTPUT}.candidates.fq

#realing the candidates to the MEI library
MosaikBuild -q ${OUTPUT}.candidates.fq -st illumina -out ${OUTPUT}.dat -quiet && \
MosaikAligner -in ${OUTPUT}.dat -out ${OUTPUT}.realigned -ia ${IA}.dat -hs 9 -mmp 0.1 -act 20 -p 8 -j ${IA}_hs9 -annpe ${ANN}/2.1.26.pe.100.0065.ann -annse ${ANN}/2.1.26.se.100.005.ann -quiet

samtools view -F 4 ${OUTPUT}.realigned.bam -h | vawk '{ if ($0 !~ /^@/) { split($1,a,":"); $1 = a[1]; $0=$0"\tTY:Z:"a[2] } print $0}' | samtools view -b - > ${OUTPUT}.hits.bam

#merge the hit reads with the anchors
samtools merge - ${OUTPUT}.anchors.bam ${OUTPUT}.hits.bam | samtools sort -n - ${OUTPUT}.merged

#get the refnames from the MEI library
samtools view ${OUTPUT}.realigned.bam  -H | grep "^@SQ" | cut -f 2 | sed -e 's/SN://g' > ${OUTPUT}.mei_refnames.txt

#filter the merged bam for anchor-mei pairs or splitters.
python ./filter_merged.py -i ${OUTPUT}.merged.bam -o /dev/stdout -m ${OUTPUT}.mei_refnames.txt | samtools view -b - > ${OUTPUT}.filtered.bam

#################

# bamToBed -cigar -i ${OUTPUT}.filtered.bam | paste - <(samtools view ${OUTPUT}.filtered.bam \
#     | vawk '{ for(i = 12; i <= NF; i++) { if($i ~ /^ME/) {print $i;} } }' ) \
#         | bedtools sort | bedtools cluster -d 350 \
#             | bedtools intersect -v -a - -b $REPMASK > ${OUTPUT}.intersect_clusters.bed

#with intersection against repmask features:
# bamToBed -cigar -i test.filtered.bam | paste - <(samtools view test.filtered.bam \
#     | vawk '{ for(i = 12; i <= NF; i++) { if($i ~ /^ME/) {print $i;} } }' ) \
#         | bedtools sort | bedtools cluster -d 350 \
#             | bedtools intersect -v -a - -b ../repmask/nchr.SLOP90.repmask_hg19_Alu.L1.SVA.ERV.bed > intersect_clusters.bed

#print only TY tag (and fields 1 and 6) from BAM
#awk '{for (i=1;i<=NF;i++) {if ($i ~/^TY:Z/) print $1"\t"$6"\t"$i;}}'