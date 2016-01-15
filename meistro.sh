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

sambamba view -t 4 -f bam -l 0 $INPUT_BAM | head -n 800000 | python ${SCRIPTS_DIR}/extract_candidates.py -a >(samtools view -b - > ${OUTPUT}.anchors.bam) -f - -c 25 -oc 7 > ${OUTPUT}.candidates.fq

MosaikBuild -q ${OUTPUT}.candidates.fq -st illumina -out ${OUTPUT}.dat -quiet && \
MosaikAligner -in ${OUTPUT}.dat -out ${OUTPUT}.realigned -ia ${IA}.dat -hs 9 -mmp 0.1 -act 20 -p 8 -j ${IA}_hs9 -annpe ${ANN}/2.1.26.pe.100.0065.ann -annse ${ANN}/2.1.26.se.100.005.ann -quiet

cat <(samtools view -H ${OUTPUT}.realigned.bam) <( zjoin -a <(samtools view -F 4 ${OUTPUT}.realigned.bam) -b <(grep "^@ERR" ${OUTPUT}.candidates.fq | sed -e 's/ /\t/g' -e 's/@//g') | cut -f 1-15,17 ) \
    | samtools view -b - > ${OUTPUT}.hits.bam

samtools merge - ${OUTPUT}.anchors.bam ${OUTPUT}.hits.bam | samtools sort -n - ${OUTPUT}.merged

python ./filter_merged.py -i ${OUTPUT}.merged.bam -o /dev/stdout | samtools view -b - > ${OUTPUT}.filtered.bam