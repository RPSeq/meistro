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

##Set output directory
RESULTS_DIR=${WORKING_DIR}/${OUTPUT_DIR_NAME}
OUTPUT=$RESULTS_DIR/$PREFIX

##Check if output dir exists and ask to overwrite
if [ -d ${RESULTS_DIR} ]; then
    echo "Output directory exists: overwriting."
else
    mkdir ${RESULTS_DIR}
fi


#   Note- These three commands run simultaneously. 
#   1. Extract MEI realignment candidates from input bam file
sambamba view -t 4 -f bam -l 0 $INPUT_BAM \
    | python ${SCRIPTS_DIR}/extract_candidates.py -a >(samtools view -b - > ${OUTPUT}.anchors.bam) -f - -c 25 -oc 7 \
        | bwa mem -t 8 -k 11 -M -C $MEI_REF /dev/stdin \
            | samtools view -F 4 -b -u - | samtools sort -n - ${OUTPUT}.realigned


samtools merge - ${OUTPUT}.anchors.bam ${OUTPUT}.realigned.bam | samtools sort -n - ${OUTPUT}.merged

python ./filter_merged.py -i ${OUTPUT}.merged.bam -o /dev/stdout | samtools view -b - > ${OUTPUT}.filtered.bam