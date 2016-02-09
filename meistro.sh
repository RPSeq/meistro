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


#Set relative dirs
WORKING_DIR=$(readlink -f $PWD)
SCRIPTS_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
FILE_LABEL=$(basename $INPUT_BAM ".bam")
FILE_DIR=$(dirname $INPUT_BAM)

#Mosaik Files
IA="/gscmnt/gc2719/halllab/users/rsmith/transposons/Mobster/Mobster-0.1.6/mobiome/54_mobiles_inclHERVK"
#IA="assembly/human_youngTE_revisedPolyA"
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

##############
### MOSAIK ###
###############

## If MOSIAK reference hasnt been built:
#MosaikBuild -fr ${IA}.fa -oa ${IA}.dat
#MosaikJump -ia ${IA}.dat -hs 9 -out ${IA}_hs9

#extract candidate reads for realingment to MEI library
sambamba view -t 4 -f bam -l 0 $INPUT_BAM | python ${SCRIPTS_DIR}/extract_candidates.py -a >(samtools view -b - > ${OUTPUT}.anchors.bam) -f - -c 20 -oc 10 > ${OUTPUT}.candidates.fq 2> /dev/null

#realing the candidates to the MEI library
MosaikBuild -q ${OUTPUT}.candidates.fq -st illumina -out ${OUTPUT}.dat -quiet && \
MosaikAligner -in ${OUTPUT}.dat -out ${OUTPUT}.realigned -ia ${IA}.dat -hs 9 -mmp 0.1 -act 20 -p 8 -j ${IA}_hs9 -annpe ${ANN}/2.1.26.pe.100.0065.ann -annse ${ANN}/2.1.26.se.100.005.ann -quiet

samtools view -F 4 ${OUTPUT}.realigned.bam -h | vawk '{ if ($0 !~ /^@/) { split($1,a,":"); $1 = a[1]; $0=$0"\tTY:Z:"a[2] } print $0}' | samtools view -b - > ${OUTPUT}.hits.bam

#merge the hit reads with the anchors
samtools merge - ${OUTPUT}.anchors.bam ${OUTPUT}.hits.bam | samtools sort -n - ${OUTPUT}.merged

#get the refnames from the MEI library
samtools view ${OUTPUT}.realigned.bam  -H | grep "^@SQ" | cut -f 2 | sed -e 's/SN://g' > ${OUTPUT}.mei_refnames.txt

#filter the merged bam for anchor-mei pairs or splitters.
python ./filter_merged.py -i ${OUTPUT}.merged.bam -o /dev/stdout -m ${OUTPUT}.mei_refnames.txt | samtools view -b - | samtools sort - ${OUTPUT}.filtered

#get clusters with bedtools cluster######################
samtools view ${OUTPUT}.filtered.bam | paste - \
    <(bamToBed -i ${OUTPUT}.filtered.bam \
            | bedtools cluster -d 300 \
                | cut -f 7) \
    | cat <(samtools view ${OUTPUT}.filtered.bam -H) - \
    | vawk '{if($0 !~ /^@/) $NF="CL:i:"$NF; print $0}' \
    | samtools view -b - > ${OUTPUT}.clusters.bam
##########################################################

bamToBed -cigar -i ${OUTPUT}.filtered.bam | paste - <(samtools view ${OUTPUT}.filtered.bam \
    | vawk '{ for(i = 12; i <= NF; i++) { if($i ~ /^RA/) {print $i;} } }' ) \
        | bedtools sort | bedtools cluster -d 300 > ${OUTPUT}.clusters.bed

bedtools intersect -v -a ${OUTPUT}.clusters.bed -b repmask/nchr.SLOP90.repmask_hg19_Alu.L1.SVA.ERV.bed > ${OUTPUT}.intersect_clusters.bed

cat ${OUTPUT}.clusters.bed | python get_clusters.py > ${OUTPUT}.filtered_clusters.bed

cat ${OUTPUT}.intersect_clusters.bed | python get_clusters.py > ${OUTPUT}.filtered_intersect_clusters.bed