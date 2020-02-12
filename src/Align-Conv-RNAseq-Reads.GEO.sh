#!/usr/bin/env bash
# This script invokes STAR to align mRNA sequences (in fastq.gz format) to
# UCSC's reference genome library (downloaded from UCSC's Genome Browser).

help_msg ()
{
	echo "Usage: "$1" [mRNA sequence stranded-ness] [number of threads]" 1>&2
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

#
# The main program begins here
#

PROG_NAME="$(basename "$0")"
if [ $# -lt 1 ]; then
	THREAD_NUMBER="$(cpus)"
elif [ $# -lt 2 ]; then
	THREAD_NUMBER="$1"
else
	help_msg "${PROG_NAME}"
fi
check_thread_number "${THREAD_NUMBER}"

# Top dataset directory and its sub-directories.
DATASET_DIR="$HOME/LINCSData/Datasets/Difference/LINCS.Dataset/Coen-Paper/Conv-GEO-Depot"
check_dir "${DATASET_DIR}"
SEQ_DIR="${DATASET_DIR}/Seqs"
check_dir "${SEQ_DIR}"
ALIGN_DIR="${DATASET_DIR}/Aligns"
check_dir "${ALIGN_DIR}"
REF_DIR="${DATASET_DIR}/References/UCSC/hg38"
check_dir "${REF_DIR}"
STAR_VER="$(STAR --version | cut -d _ -f 2)"
INDEX_DIR="${REF_DIR}/STAR-Index-${STAR_VER}"
check_dir "${INDEX_DIR}"

# Get a list of all sequence files.
SEQ_FILE_SUFFIX="fastq.gz"
readarray -t -d $'\0' SEQ_FILES < <(find -L "${SEQ_DIR}" -maxdepth 1 -type f -name "*\.${SEQ_FILE_SUFFIX}" -print0)

# Generate a list of alignment file names.
ALIGN_FILE_SUFFIX="bam"
ALIGN_FILES=()
for SEQ_FILE in "${SEQ_FILES[@]}"; do
	SEQ_FILE_NAME="$(basename "${SEQ_FILE}")"
	ALIGN_FILE_NAME="$(echo ${SEQ_FILE_NAME} | sed "s/${SEQ_FILE_SUFFIX}$/${ALIGN_FILE_SUFFIX}/g")"
	ALIGN_FILE="${ALIGN_DIR}/${ALIGN_FILE_NAME}"
	ALIGN_FILES+=("${ALIGN_FILE}")
done


SEQ_FILE_SUFFIX="fastq.gz"
UNZIP_CMD="gzip -dc"
STAR_ALIGN_FILE="Aligned.out.bam"

# Start aligning paired-end sequence to reference library using STAR.
((IDX=0))
for SEQ_FILE in "${SEQ_FILES[@]}"; do
	echo "STAR is aligning ${SEQ_FILE} to ${INDEX_DIR} ..."
	STAR \
		--runMode alignReads \
		--runThreadN "${THREAD_NUMBER}" \
		--readFilesIn "${SEQ_FILE}" \
		--readFilesCommand "${UNZIP_CMD}" \
		--quantMode TranscriptomeSAM GeneCounts \
		--genomeDir "${INDEX_DIR}" \
		--outFileNamePrefix "${ALIGN_DIR}/" \
		--outSAMtype BAM Unsorted \
		--outFilterMultimapNmax 10 \
		--outFilterMismatchNmax 10 \
		--outMultimapperOrder Old_2.4 \
		--outSAMunmapped Within \
		--outSAMattributes NH HI NM MD AS nM jI jM
	mv "${ALIGN_DIR}/${STAR_ALIGN_FILE}" "${ALIGN_FILES[${IDX}]}"
	((IDX+=1))
done

exit $?
