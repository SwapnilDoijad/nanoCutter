#!/bin/bash
###############################################################################
## inputs
	current_dir=$(pwd)
	input_dir=test_data/fasta
	# input_dir=test_data/fastq
	output_dir=$current_dir/test_results

	ls $input_dir > list.txt
	list=list.txt
###############################################################################

	source venv/bin/activate

	mkdir -p $output_dir/fastq > /dev/null 2>&1
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
		--coords-dir $output_dir/raw_files/$i \
		--out-dir $output_dir/fastq \
		--threads 40

		python3.9 scripts/count.py \
		--input-file $input_dir/$i \
		--input-type fastq \
		--threads 40

	done

	deactivate
###############################################################################
