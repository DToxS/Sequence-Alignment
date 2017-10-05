#!/usr/bin/env bash

check_program ()
{
	if [ -z "$(which $1)" ]; then
		echo "ERROR: "$1" is not found!" 1>&2
		exit 1
	fi
}

# 1 Parameters

# 1.1 Global

TOP_DIR="$HOME/LINCSData"

# 1.2 Dataset
DATA_DIR="${TOP_DIR}/Datasets/Alignment/LINCS.Dataset/LINCS.Dataset.Gene.LINCS.20170710"
SEQ_DIR="${DATA_DIR}/Seqs"
ALIGN_DIR="${DATA_DIR}/Aligns"
COUNT_DIR="${DATA_DIR}/Counts"
SAMPLE_ID="DGE-Plate-6"
WELL_NUMBERS=()
for LETTER in {A..P}; do
	WELL_NUMBERS=(${WELL_NUMBERS[@]} $(seq 1 24 | sed "s/\([[:digit:]]\+\)/$LETTER\1/g"))
done

# 1.3 Reference
REF_DIR="${TOP_DIR}/References/Broad_UMI"
SPECIES_DIR="${REF_DIR}/Human_RefSeq"
REF_SEQ_FILE="${SPECIES_DIR}/refMrna_ERCC_polyAstrip.hg19.fa"
SYM2REF_FILE="${SPECIES_DIR}/refGene.hg19.sym2ref.dat"
ERCC_SEQ_FILE="${REF_DIR}/ERCC92.fa"
BARCODE_DIR="${REF_DIR}/barcodes_trugrade_384_set1"

# 1.4 Program
BWA_PROG="bwa"
PYTHON_PROG="python"
PROG_DIR="${TOP_DIR}/Programs/Broad-DGE"
MERGE_COUNT_PROG="${PROG_DIR}/merge_and_count.py"
BWA_ALN_SEED_LENGTH=24
BWA_SAM_MAX_ALIGNS_FOR_XA_TAG=20
THREAD_NUMBER="4"

check_program "${BWA_PROG}"
check_program "${PYTHON_PROG}"

# 2 Computation

for WELL_NUMBER in ${WELL_NUMBERS[@]}; do
	echo "Sample" "${WELL_NUMBER}"
	BARCODE_FILE="${BARCODE_DIR}/barcodes_trugrade_384_set1.${WELL_NUMBER}.dat"
	SAMPLE_ALIGN_DIR="${ALIGN_DIR}/${WELL_NUMBER}"
	if [ ! -d "${SAMPLE_ALIGN_DIR}" ]; then
		mkdir -p "${SAMPLE_ALIGN_DIR}"
	fi
	SEQ_FILE="${SEQ_DIR}/${SAMPLE_ID}.${WELL_NUMBER}.fastq"
	SAM_FILE="${SAMPLE_ALIGN_DIR}/${SAMPLE_ID}.${WELL_NUMBER}.sam"
	# 2.1 Alignment
	# Align sequence fragments to reference genome library.
	"${BWA_PROG}" aln -l "${BWA_ALN_SEED_LENGTH}" -t "${THREAD_NUMBER}" "${REF_SEQ_FILE}" "${SEQ_FILE}" | "${BWA_PROG}" samse -n "${BWA_SAM_MAX_ALIGNS_FOR_XA_TAG}" "${REF_SEQ_FILE}" - "${SEQ_FILE}" > "${SAM_FILE}"
	# 2.2 Counting
	# Count the number of sequence alignments for reference genes.
	SAMPLE_COUNT_DIR="${COUNT_DIR}/${WELL_NUMBER}"
	if [ ! -d "${SAMPLE_COUNT_DIR}" ]; then
		mkdir -p "${SAMPLE_COUNT_DIR}"
	fi
	"${PYTHON_PROG}" "${MERGE_COUNT_PROG}" "${SAMPLE_ID}" "${SYM2REF_FILE}" "${ERCC_SEQ_FILE}" "${BARCODE_FILE}" "${SAMPLE_ALIGN_DIR}/" "${SAMPLE_COUNT_DIR}/" False
done

exit 0
