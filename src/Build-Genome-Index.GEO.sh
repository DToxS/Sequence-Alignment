#!/usr/bin/env bash
# Build STAR index for reference genome library.

#
# Function definitions
#

help_msg ()
{
	echo "Usage: "$1" [Ref Index Dir] [Genome Seq File] [Genome Annot File] [Read Length] [Max Threads]" 1>&2
	echo "       [Ref Index Dir] is a directory of the output STAR index files for reference genome library (Value: directory)" 1>&2
	echo "       [Genome Seq File] is a genome sequence file in FASTA format for reference genome library (Value: file)" 1>&2
	echo "       [Genome Annot File] is a genome annotation file in GTF format for reference genome library (Value: [no file], file)" 1>&2
	echo "       [Read Length] is the read length of a sequence fragment (Value: 1, ..., [100], ...)" 1>&2
	echo "       [Max Threads] is the maximum number of threads to use (Value: [1], 2, ...)" 1>&2
	return 1
}

check_file ()
{
	if [ ! -f "$1" ]; then
		echo "ERROR: "$1" is not found!" 1>&2
		return 1
	else
		return 0
	fi
}

create_dir ()
{
	if [ ! -d "$1" ]; then
		echo "Creating $1"
		mkdir -p "$1"
		return $?
	fi
	return 0
}

check_read_length ()
{
	if [ "$1" -lt 1 ]; then
		echo "ERROR: read length must be greater than zero!" 1>&2
		return 1
	else
		return 0
	fi
}

check_n_threads ()
{
	if [ "$1" -lt 1 ]; then
		echo "ERROR: the maximum number of threads must be greater than 0!" 1>&2
		return 1
	fi
}

#
# Main program
#

# Initialize error code.
EXIT_CODE=0

# Obtain program path.
PROG_PATH="$0"
PROG_DIR="$(dirname "${PROG_PATH}")"
PROG_NAME="$(basename "${PROG_PATH}")"

# Specify the range of the number of input arguments.
N_ARGS_MIN=2
N_ARGS_MAX=5
if [ $# -ge "${N_ARGS_MIN}" ] && [ $# -le "${N_ARGS_MAX}" ]; then
	ANNOT_FILE_PATH=""
	READ_LENGTH="100"
	N_THREADS_MAX="1"
	if [ ${EXIT_CODE} == 0 ] && [ $# -ge $((N_ARGS_MIN)) ] && [ "${N_ARGS_MIN}" -ge 1 ]; then
		REF_INDEX_DIR="${@:((N_ARGS_MIN-1)):1}"
		SEQ_FILE_PATH="${@:((N_ARGS_MIN)):1}"
		check_file "${SEQ_FILE_PATH}" && create_dir "${REF_INDEX_DIR}"
		EXIT_CODE=$?
	fi
	if [ ${EXIT_CODE} == 0 ] && [ $# -ge $((N_ARGS_MIN+1)) ]; then
		ANNOT_FILE_PATH="${@:((N_ARGS_MIN+1)):1}"
		if [ "${ANNOT_FILE_PATH}" != "" ]; then
			check_file "${ANNOT_FILE_PATH}"
			EXIT_CODE=$?
		fi
	fi
	if [ ${EXIT_CODE} == 0 ] && [ $# -ge $((N_ARGS_MIN+2)) ]; then
		READ_LENGTH="${@:((N_ARGS_MIN+2)):1}"
		check_read_length "${READ_LENGTH}"
		EXIT_CODE=$?
	fi
	if [ ${EXIT_CODE} == 0 ] && [ $# -ge $((N_ARGS_MIN+3)) ]; then
		N_THREADS_MAX="${@:((N_ARGS_MIN+3)):1}"
		check_n_threads "${N_THREADS_MAX}"
		EXIT_CODE=$?
	fi
else
	help_msg "$(basename "$0")"
	EXIT_CODE=$?
fi

# Set STAR program.
STAR_PROG_NAME="STAR"
# Set the name of the program for obtaining the number of processors.
GET_CPUS_PROG_NAME="get-cpus.sh"

# Get the path of STAR program.
if [ ${EXIT_CODE} == 0 ]; then
	STAR_PROG_PATH="$(which "${STAR_PROG_NAME}")"
	EXIT_CODE=$?
	if [ "${EXIT_CODE}" != 0 ]; then
		echo "${STAR_PROG_NAME} is not found!" 1>&2
	fi
fi

# Set overhang length.
if [ ${EXIT_CODE} == 0 ]; then
	OVERHANG_LENGTH=$((READ_LENGTH-1))
	EXIT_CODE=$?
fi

# Check the path of the program for obtaining the number of processors.
if [ ${EXIT_CODE} == 0 ]; then
	GET_CPUS_PROG_PATH="$(which "${GET_CPUS_PROG_NAME}")"
	EXIT_CODE=$?
	if [ "${EXIT_CODE}" != 0 ]; then
		echo "${GET_CPUS_PROG_NAME} is not found!" 1>&2
	fi
fi

# Set the number of processors.
if [ ${EXIT_CODE} == 0 ]; then
	N_THREADS="$("${GET_CPUS_PROG_PATH}" "${N_THREADS_MAX}")"
	EXIT_CODE=$?
fi

# Build STAR index.
if [ ${EXIT_CODE} == 0 ]; then
	# Start building STAR index for reference genome library.
	echo "Building ${STAR_PROG_NAME} index of reference genome..."
	"${STAR_PROG_PATH}" \
		--runMode genomeGenerate \
		--genomeDir "${REF_INDEX_DIR}" \
		--genomeFastaFiles "${SEQ_FILE_PATH}" \
		$([ "${ANNOT_FILE_PATH}" != "" ] && echo "--sjdbGTFfile ${ANNOT_FILE_PATH} --sjdbOverhang ${OVERHANG_LENGTH}" || :) \
		--runThreadN "${N_THREADS}"
	EXIT_CODE=$?
fi

# Print a warning message.
if [ ${EXIT_CODE} != 0 ]; then
	echo "Abort executing program due to error!" 1>&2
fi

# Exit with error code.
exit ${EXIT_CODE}
