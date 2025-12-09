#!/usr/bin/env python3
"""Small utility to count sequences/reads in FASTA/FASTQ files (supports .gz).

Usage examples:
  python3 scripts/count.py --input-file reads.fastq.gz --input-type fastq --out-dir results/
  python3 scripts/count.py --input-file seqs.fasta --out-dir .

Writes a file <basename>.count.txt in --out-dir and prints the count to stdout.
"""

import argparse
import gzip
import os
from Bio import SeqIO
import shutil
import subprocess


def _open_text(path):
	"""Open a file path for text reading; handles .gz automatically."""
	if path.lower().endswith(".gz"):
		return gzip.open(path, "rt")
	return open(path, "r")


def guess_format_from_extension(path):
	"""Return 'fastq' or 'fasta' if recognized, otherwise None."""
	p = path.lower()
	if p.endswith(".gz"):
		p = p[:-3]
	if p.endswith((".fastq", ".fq")):
		return "fastq"
	if p.endswith((".fasta", ".fa", ".fna")):
		return "fasta"
	return None


def count_sequences(path, fmt=None):
	"""Count sequences/reads using Bio.SeqIO.parse.

	Parameters
	- path: filesystem path (can be .gz)
	- fmt: optional Biopython format string ('fastq' or 'fasta'); if None, guessed from extension

	Returns: integer count
	"""
	if fmt is None:
		fmt = guess_format_from_extension(path)
		if fmt is None:
			raise ValueError("Can't guess format from extension; pass --input-type")
	with _open_text(path) as handle:
		return sum(1 for _ in SeqIO.parse(handle, fmt))


def _strip_known_extensions(name):
	"""Return base name without common sequence extensions (handles .gz combos)."""
	lower = name.lower()
	for ext in (".fastq.gz", ".fq.gz", ".fastq", ".fq", ".fasta", ".fa", ".fna"):
		if lower.endswith(ext):
			return name[: len(name) - len(ext)]
	return os.path.splitext(name)[0]


def main():
	parser = argparse.ArgumentParser(
		description="Count sequences/reads in FASTA/FASTQ files. Supports .gz files."
	)
	parser.add_argument(
		"--input-type",
		choices=["fastq", "fasta", "auto"],
		default="auto",
		help="format of input file; 'auto' will guess from filename",
	)
	parser.add_argument("--input-file", required=True, help="input file path (can be .gz)")
	parser.add_argument("--threads", type=int, default=1, help="number of threads to use for fast (external) counting")

	args = parser.parse_args()

	fmt = None if args.input_type == "auto" else args.input_type
	threads = max(1, int(args.threads))

	# Fast paths using external tools when possible
	is_gz = args.input_file.lower().endswith('.gz')
	base = _strip_known_extensions(os.path.basename(args.input_file))

	cnt = None
	# For FASTQ: count lines and divide by 4. For FASTA: count headers starting with '>'.
	if fmt == 'fastq' or (fmt is None and guess_format_from_extension(args.input_file) == 'fastq'):
		try:
			if is_gz:
				# Prefer pigz for threaded decompression
				pigz = shutil.which('pigz')
				gzipexe = shutil.which('gzip')
				if threads > 1 and pigz:
					p1 = subprocess.Popen([pigz, '-dc', '-p', str(threads), args.input_file], stdout=subprocess.PIPE)
					p2 = subprocess.Popen(['wc', '-l'], stdin=p1.stdout, stdout=subprocess.PIPE)
					p1.stdout.close()
					out = p2.communicate()[0].decode().strip()
					lines = int(out.split()[0])
					cnt = lines // 4
				elif gzipexe:
					p1 = subprocess.Popen([gzipexe, '-dc', args.input_file], stdout=subprocess.PIPE)
					p2 = subprocess.Popen(['wc', '-l'], stdin=p1.stdout, stdout=subprocess.PIPE)
					p1.stdout.close()
					out = p2.communicate()[0].decode().strip()
					lines = int(out.split()[0])
					cnt = lines // 4
			else:
				# uncompressed: fast wc -l
				wc = shutil.which('wc')
				if wc:
					out = subprocess.check_output([wc, '-l', args.input_file]).decode().strip()
					lines = int(out.split()[0])
					cnt = lines // 4
		except Exception:
			cnt = None

	elif fmt == 'fasta' or (fmt is None and guess_format_from_extension(args.input_file) == 'fasta'):
		try:
			if is_gz:
				pigz = shutil.which('pigz')
				gzipexe = shutil.which('gzip')
				grep = shutil.which('grep')
				if threads > 1 and pigz and grep:
					p1 = subprocess.Popen([pigz, '-dc', '-p', str(threads), args.input_file], stdout=subprocess.PIPE)
					p2 = subprocess.Popen([grep, '-c', '^>'], stdin=p1.stdout, stdout=subprocess.PIPE)
					p1.stdout.close()
					out = p2.communicate()[0].decode().strip()
					cnt = int(out)
				elif gzipexe and grep:
					p1 = subprocess.Popen([gzipexe, '-dc', args.input_file], stdout=subprocess.PIPE)
					p2 = subprocess.Popen([grep, '-c', '^>'], stdin=p1.stdout, stdout=subprocess.PIPE)
					p1.stdout.close()
					out = p2.communicate()[0].decode().strip()
					cnt = int(out)
			else:
				grep = shutil.which('grep')
				if grep:
					out = subprocess.check_output([grep, '-c', '^>', args.input_file]).decode().strip()
					cnt = int(out)
		except Exception:
			cnt = None

	# Fallback to SeqIO parsing if no fast path succeeded
	if cnt is None:
		try:
			cnt = count_sequences(args.input_file, fmt)
		except Exception as e:
			parser.error(str(e))

	# Print only the file basename without sequence suffix and the total count
	print(f"{base}\t{cnt}")


if __name__ == "__main__":
	main()