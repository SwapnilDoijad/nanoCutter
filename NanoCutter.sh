#!/bin/bash
###############################################################################
# Usage/help
show_usage() {
	cat <<'USAGE'
Usage: NanoCutter.sh [OPTIONS]

Options:
	-h, --help               Show this help message and exit
	-c, --current-dir DIR    Set the current working directory (default: pwd)
	-i, --input-dir DIR      Set the input directory containing FASTA/FASTQ files
	-o, --output-dir DIR     Set the output directory for results

Examples:
	./NanoCutter.sh -i test_data/fasta -o test_results

USAGE
}

while [[ $# -gt 0 ]]; do
	case "$1" in
		-h|--help)
			show_usage
			exit 0
			;;
		-c|--current-dir)
			if [[ -n "$2" && ${2:0:1} != "-" ]]; then
				current_dir="$2"
				shift 2
			else
				echo "Error: $1 requires a value" >&2
				exit 1
			fi
			;;
		--current-dir=*)
			current_dir="${1#*=}"
			shift
			;;
		-i|--input-dir)
			if [[ -n "$2" && ${2:0:1} != "-" ]]; then
				input_dir="$2"
				shift 2
			else
				echo "Error: $1 requires a value" >&2
				exit 1
			fi
			;;
		--input-dir=*)
			input_dir="${1#*=}"
			shift
			;;
		-o|--output-dir)
			if [[ -n "$2" && ${2:0:1} != "-" ]]; then
				output_dir="$2"
				shift 2
			else
				echo "Error: $1 requires a value" >&2
				exit 1
			fi
			;;
		--output-dir=*)
			output_dir="${1#*=}"
			shift
			;;
		--)
			shift
			break
			;;
		*)
			echo "Unknown option: $1" >&2
			show_usage
			exit 1
			;;
	esac
done

# Basic validation
if [[ ! -d "$input_dir" ]]; then
	echo "Error: input directory '$input_dir' does not exist." >&2
	exit 1
fi

mkdir -p $output_dir/fastq > /dev/null 2>&1

ls "$input_dir" > list.txt
list=list.txt
###############################################################################
	> $output_dir/read_counts.txt
	source venv/bin/activate

	for i in $(cat $list); do
		echo $i 

		mkdir -p $output_dir/$i > /dev/null 2>&1

		python3.9 scripts/Collapser.py \
		$input_dir/$i \
		$i \
		--out-dir $output_dir/$i \
		--dotplot \
		--rc-repeats \
		--direct-repeats \
		--bin-size 300 \
		--max-gap 1000 \
		--min-pct-id 90 \
		--min-length 2000 \
		--jobs 1 \
		--workers 10 

		python3.9 scripts/cut_sequences.py \
		--input-type fastq \
		--input-file $input_dir/$i \
		--coords-dir $output_dir/$i \
		--out-dir $output_dir/fastq \
		--threads 40

		python3.9 scripts/count.py \
		--input-file $input_dir/$i \
		--input-type fastq \
		--threads 40 \
		>> $output_dir/read_counts.txt

	done

	deactivate
###############################################################################
