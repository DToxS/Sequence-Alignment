#!/usr/bin/env bash
# This script aligns conventional mRNA sequence files.

# Set top directory.
DATASET_DIR="${HOME}/LINCSData/Datasets/Difference/LINCS.Dataset/Coen-Paper/Conv-GEO-Depot"

# Set executable program.
PROG_DIR_PATH="${DATASET_DIR}/Programs/Sequence-Alignment"
PROG_FILE_NAME="Align-RNAseq-Reads.GEO.sh"
PROG_FILE_PATH="${PROG_DIR_PATH}/${PROG_FILE_NAME}"

# Set input directories and files.
SEQ_DIR="${DATASET_DIR}/Seqs"
ALIGN_DIR="${DATASET_DIR}/Aligns"
REF_DIR="${DATASET_DIR}/References/UCSC/hg38"
STAR_VER="$(STAR --version | cut -d _ -f 2)"
REF_INDEX_DIR="${REF_DIR}/STAR-Index-${STAR_VER}-Conv"

# Set input arguments.
SEQ_FILE_SUFFIX="fastq.gz"
THREAD_NUMBER="16"

# Run the program.
echo "${PROG_FILE_PATH} ${SEQ_DIR} ${ALIGN_DIR} ${REF_INDEX_DIR} ${SEQ_FILE_SUFFIX} ${THREAD_NUMBER}"
"${PROG_FILE_PATH}" ${SEQ_DIR} ${ALIGN_DIR} ${REF_INDEX_DIR} "${SEQ_FILE_SUFFIX}" "${THREAD_NUMBER}"

# Exit with error code.
exit $?
