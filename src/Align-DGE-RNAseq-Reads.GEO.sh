#!/usr/bin/env bash
# This script aligns 3'-DGE mRNA sequence files.

# Initialize error code.
EXIT_CODE=0

# Set top directory.
DATASET_DIR="${HOME}/LINCSData/Datasets/Difference/LINCS.Dataset/Coen-Paper/DGE-GEO-Depot"

# Set executable program.
PROG_DIR_PATH="${DATASET_DIR}/Programs/Sequence-Alignment"
PROG_FILE_NAME="Align-RNAseq-Reads.GEO.sh"
PROG_FILE_PATH="${PROG_DIR_PATH}/${PROG_FILE_NAME}"

# Set the directory of genome reference library.
REF_DIR="${DATASET_DIR}/References/UCSC/hg38"
STAR_VER="$(STAR --version | cut -d _ -f 2)"
REF_INDEX_DIR="${REF_DIR}/STAR-Index-${STAR_VER}-DGE"

# Set input arguments.
SEQ_FILE_SUFFIX="fastq.gz"
THREAD_NUMBER="16"

# Set the names of dataset series.
DATASET_SERIES_NAMES=("Dataset-20150409" "Dataset-20150503" "Dataset-20151120" "Dataset-20161108")

# Process all datasets.
for DATASET_SERIES_NAME in "${DATASET_SERIES_NAMES[@]}"; do
    # Set the directory of current dataset series.
	DATASET_SERIES_DIR="${DATASET_DIR}/${DATASET_SERIES_NAME}"
	if [ -d "${DATASET_SERIES_DIR}" ]; then
		# Set input directories and files.
		SEQ_DIR="${DATASET_SERIES_DIR}/Seqs"
		ALIGN_DIR="${DATASET_SERIES_DIR}/Aligns"
		# Run the program.
		echo "${PROG_FILE_PATH} ${SEQ_DIR} ${ALIGN_DIR} ${REF_INDEX_DIR} ${SEQ_FILE_SUFFIX} ${THREAD_NUMBER}"
		"${PROG_FILE_PATH}" ${SEQ_DIR} ${ALIGN_DIR} ${REF_INDEX_DIR} "${SEQ_FILE_SUFFIX}" "${THREAD_NUMBER}"
		EXIT_CODE=$?
		if [ ${EXIT_CODE} -ne 0 ]; then
			echo "An error occurred in aligning the dataset ${DATASET_SERIES_DIR}!" 1>&2
			break
		fi
	else
		echo "WARNING: the dataset directory ${DATASET_SERIES_DIR} is not found!" 1>&2
	fi
done

# Exit with error code.
exit ${EXIT_CODE}
