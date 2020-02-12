#!/usr/bin/env bash
# Build the index of reference genome used by STAR.

help_msg ()
{
	echo "Usage: "$1" [Read Length] [number of threads]" 1>&2
	exit 1
}

cpus ()
{
	local MAX_CPUS=4
	local CPUS=`getconf _NPROCESSORS_ONLN 2>/dev/null`
	[ -z "$CPUS" ] && CPUS=`getconf NPROCESSORS_ONLN`
	[ -z "$CPUS" ] && CPUS=`ksh93 -c 'getconf NPROCESSORS_ONLN'`
	[ -z "$CPUS" ] && CPUS=1
	[ "$CPUS" -gt "$MAX_CPUS" ] && CPUS=$MAX_CPUS
	echo "${CPUS}"
}

check_read_length ()
{
	if [ "$1" -lt 1 ]; then
		echo "ERROR: Read length must be greater than zero!" 1>&2
		exit 1
	fi
}

check_thread_number ()
{
	if [ "$1" -lt 1 ]; then
		echo "ERROR: The number of threads must be greater than zero!" 1>&2
		exit 1
	fi
}

check_dir ()
{
	if [ ! -d "$1" ]; then
		echo "ERROR: Directory "$1" is not found!" 1>&2
		exit 1
	fi
}

check_file ()
{
	if [[ ! -f "$1" && ! -h "$1" ]] || [[ ! -r "$1" ]]; then
		echo "ERROR: File "$1" is not found or not accessible!" 1>&2
		exit 1
	fi
}

#
# The main program begins here
#

PROG_NAME="$(basename "$0")"
if [ $# -lt 1 ]; then
	READ_LENGTH="100"
	THREAD_NUMBER="$(cpus)"
elif [ $# -lt 2 ]; then
	READ_LENGTH="$1"
	THREAD_NUMBER="$(cpus)"
elif [ $# -lt 3 ]; then
	READ_LENGTH="$1"
	THREAD_NUMBER="$2"
else
	help_msg "${PROG_NAME}"
fi
check_read_length "${READ_LENGTH}"
check_thread_number "${THREAD_NUMBER}"

# Top dataset directory and its sub-directories.
DATASET_DIR="$HOME/LINCSData/Datasets/Difference/LINCS.Dataset/Coen-Paper/Conv-GEO-Depot"
check_dir "${DATASET_DIR}"
REF_DIR="${DATASET_DIR}/References/UCSC/hg38"
check_dir "${REF_DIR}"
SEQ_DIR="${REF_DIR}/Sequence"
check_dir "${SEQ_DIR}"
SEQ_FILE="${SEQ_DIR}/genome.hg38.fa"
check_file "${SEQ_FILE}"
ANNOT_DIR="${REF_DIR}/Annotation"
check_dir "${ANNOT_DIR}"
ANNOT_FILE="${ANNOT_DIR}/RefSeq.hg38.gtf"
check_file "${ANNOT_FILE}"
STAR_VER="$(STAR --version | cut -d _ -f 2)"
INDEX_DIR="${REF_DIR}/STAR-Index-${STAR_VER}-Test"

# Create index directory as needed.
if [ ! -d "${INDEX_DIR}" ]; then
	echo "Creating ${INDEX_DIR}"
	mkdir -p "${INDEX_DIR}"
fi

# Build the index data of reference genome used by STAR.
echo "Building the STAR index of reference genome..."
STAR \
	--runMode genomeGenerate \
	--genomeDir "${INDEX_DIR}" \
	--genomeFastaFiles "${SEQ_FILE}" \
	--sjdbGTFfile "${ANNOT_FILE}" \
	--sjdbOverhang "${READ_LENGTH}" \
	--runThreadN "${THREAD_NUMBER}"

exit 0
