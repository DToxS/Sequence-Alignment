#!/usr/bin/env bash
# Build the STAR index of human reference genome library for conventional
# mRNA sequencing data.

# Initialize error code.
EXIT_CODE=0

# Set the top dataset directory.
DATASET_DIR="$HOME/LINCSData/Datasets/Difference/LINCS.Dataset/Coen-Paper/Conv-GEO-Depot"

# Set executable program.
PROG_DIR_PATH="${DATASET_DIR}/Programs/Sequence-Alignment"
PROG_FILE_NAME="Build-Genome-Index.GEO.sh"
PROG_FILE_PATH="${PROG_DIR_PATH}/${PROG_FILE_NAME}"

# Set the directories and files of reference genome library..
HG_REF_DIR="${DATASET_DIR}/References/UCSC/hg38"
SEQ_FILE_DIR="${HG_REF_DIR}/Sequence"
SEQ_FILE_NAME="genome.hg38.fa"
SEQ_FILE_PATH="${SEQ_FILE_DIR}/${SEQ_FILE_NAME}"
ANNOT_FILE_DIR="${HG_REF_DIR}/Annotation"
ANNOT_FILE_NAME="RefSeq.hg38.gtf"
ANNOT_FILE_PATH="${ANNOT_FILE_DIR}/${ANNOT_FILE_NAME}"
STAR_VER="$(STAR --version | cut -d _ -f 2)"
REF_INDEX_DIR_NAME="STAR-Index-${STAR_VER}-Conv"
REF_INDEX_DIR_PATH="${HG_REF_DIR}/${REF_INDEX_DIR_NAME}"

# Set the input arguments.
# Conventional mRNA sequence data files have a read length of 100.
READ_LENGTH="100"
THREAD_NUMBER="16"

# Create the index directory as needed.
if [ ! -d "${REF_INDEX_DIR_PATH}" ]; then
	echo "Creating ${REF_INDEX_DIR_PATH} ..."
	mkdir "${REF_INDEX_DIR_PATH}"
	EXIT_CODE=$?
	if [ ${EXIT_CODE} -ne 0 ]; then
		echo "ERROR: cannot create the index directory ${REF_INDEX_DIR_PATH}!" 1>&2
	fi
else
	echo "ERROR: the index directory ${REF_INDEX_DIR_PATH} already exists!" 1>&2
	EXIT_CODE=1
fi

# Build the index data of reference genome used by STAR.
if [ ${EXIT_CODE} -eq 0 ]; then
	echo "Building the STAR index of human reference genome library for conventional mRNA-seq data ..."
	"${PROG_FILE_PATH}" ${REF_INDEX_DIR_PATH} ${SEQ_FILE_PATH} "${ANNOT_FILE_PATH}" "${READ_LENGTH}" "${THREAD_NUMBER}"
fi

# Exit with error code.
exit ${EXIT_CODE}
