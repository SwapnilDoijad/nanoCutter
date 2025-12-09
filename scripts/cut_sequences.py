#!/usr/bin/env python3
"""
cut_sequences.py

Utility to cut a sequence (FASTA/FASTQ) using a coordinates TSV (1-based inclusive).
Writes the fragments corresponding to intervals into a single combined output file by default.

Usage examples:
    # Default: write combined fragments into a single multi-record output file
    python3 scripts/cut_sequences.py --input-file data_1_116_O/0202f909-770f-433c-a08c-3b4b649a9988.fasta \
            --coords results/0202f909-770f-433c-a08c-3b4b649a9988/0202f909-770f-433c-a08c-3b4b649a9988.repeats.cut_coordinates.tsv \
            --out-dir results/0202f909-770f-433c-a08c-3b4b649a9988

    # The script writes a single combined fragments file by default
    python3 scripts/cut_sequences.py --input-file input.fasta --input-type fasta --coords cuts.tsv

This script prefers to load the sequence into memory (fast for reads up to many MBs).
For very large genomes or many random-access slices, prefer `samtools faidx` or `pyfaidx`.
"""
import argparse
from pathlib import Path
import gzip
import multiprocessing
import tempfile
import shutil
import os
try:
    from Bio import SeqIO
    _HAS_BIO = True
except Exception:
    SeqIO = None
    _HAS_BIO = False


def _read_first_record_no_bio(path):
    """Minimal FASTA reader to extract the first record id and sequence.

    This is a lightweight fallback for environments without Biopython.
    It assumes the FASTA is text and returns (id, seq) for the first record.
    """
    seq_lines = []
    seq_header = None
    with open(path, 'r') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            if line.startswith('>'):
                if seq_header is not None:
                    break
                # header; keep the full header line (without leading '>')
                seq_header = line[1:].strip()
                continue
            if seq_header is None:
                # skip leading non-header lines
                continue
            seq_lines.append(line.strip())
    if seq_header is None:
        return None, ''
    return seq_header, ''.join(seq_lines)


def _read_first_record_fastq_no_bio(path):
    """Minimal FASTQ reader to extract the first record id, sequence and quality list.

    Returns (id, seq, qual_list) where qual_list is a list of ints (phred scores).
    This is a lightweight fallback for environments without Biopython.
    """
    with open(path, 'r') as fh:
        # read lines until we find a non-empty header starting with '@'
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            if line.startswith('@'):
                hdr = line[1:].strip()
                break
        else:
            return None, '', []

        # read sequence, plus line, quality line
        seq = ''
        qual = ''
        try:
            # sequence line(s) are usually a single line; we read next non-empty
            seq = fh.readline().rstrip('\n').strip()
            plus = fh.readline().rstrip('\n')
            qual = fh.readline().rstrip('\n').strip()
        except Exception:
            return None, '', []

    # convert quality string to list of ints (phred+33)
    qual_list = [ord(c) - 33 for c in qual]
    return hdr, seq, qual_list


def parse_coords(path):
    """Read a tab-delimited coords file with header '#start\tend' and return
    a sorted list of non-overlapping 1-based inclusive intervals.

    NOTE: do not merge adjacent/contiguous intervals. If two intervals are
    directly adjacent (cur_b + 1 == a) we keep them separate so downstream
    fragment writers will produce distinct fragments rather than a merged
    full-length slice.
    """
    intervals = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                a = int(parts[0])
                b = int(parts[1])
            except Exception:
                continue
            if a > b:
                a, b = b, a
            intervals.append((a, b))
    if not intervals:
        return []
    intervals.sort(key=lambda x: x[0])
    merged = []
    cur_a, cur_b = intervals[0]
    for a, b in intervals[1:]:
        # merge only if intervals overlap; do NOT merge adjacent/contiguous
        # intervals because users may want explicit fragment boundaries.
        if a <= cur_b:
            cur_b = max(cur_b, b)
        else:
            merged.append((cur_a, cur_b))
            cur_a, cur_b = a, b
    merged.append((cur_a, cur_b))
    return merged


def write_fragments(seq_id, seq, intervals, out_dir, prefix, quals=None, is_fastq=False):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for idx, (a, b) in enumerate(intervals, start=1):
        # convert to 0-based slice
        s0 = a - 1
        s1 = b
        sub = seq[s0:s1]
        if is_fastq:
            # prefix can be tuple (base, ext) or base string
            if isinstance(prefix, tuple):
                base, ext = prefix
            else:
                base, ext = prefix, '.fastq'
            fname = out_dir / f"{base}.part{idx}{ext}"
            # slice quality list
            sub_quals = None
            if quals is not None:
                sub_quals = quals[s0:s1]
            # convert to quality string
            qstr = ''.join(chr(q + 33) for q in sub_quals) if sub_quals is not None else ('~' * len(sub))
            with open(fname, 'w') as fh:
                fh.write(f"@{seq_id}_{a}_{b}\n")
                fh.write(str(sub) + "\n")
                fh.write("+\n")
                fh.write(qstr + "\n")
        else:
            if isinstance(prefix, tuple):
                base, ext = prefix
            else:
                base, ext = prefix, '.fasta'
            fname = out_dir / f"{base}.part{idx}{ext}"
            with open(fname, 'w') as fh:
                # preserve the original header and append a coordinate suffix
                fh.write(f">{seq_id}_{a}_{b}\n")
                fh.write(str(sub) + "\n")
        written.append(str(fname))
    return written


def write_fragments_combined(seq_id, seq, intervals, out_dir, prefix, quals=None, is_fastq=False):
    """Write all fragments into a single multi-record FASTA.

    Header format: >{seq_id}_{start}_{end}
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    # prefix expected to already include the '.collapsed' suffix; prefix is used
    # without adding extra '.fragments' text so final files look like
    # {inputname}.collapsed.fasta or {inputname}.collapsed.fastq(.gz)
    # 'prefix' here is the base filename (including .collapsed) and the
    # function caller should supply the right extension via prefix_ext.
    # Determine extension from prefix if provided (prefix may be a tuple),
    # but for simplicity we'll accept prefix as a tuple (base, ext) or a
    # plain string (base) and expect callers to manage extensions.
    if isinstance(prefix, tuple):
        base, ext = prefix
    else:
        # caller provided base string; we don't know ext here, default to .fasta
        base, ext = prefix, '.fasta'

    fname = out_dir / f"{base}{ext}"
    if ext.endswith('.fastq') or ext.endswith('.fq'):
        with open(fname, 'w') as fh:
            for (a, b) in intervals:
                s0 = a - 1
                s1 = b
                sub = seq[s0:s1]
                sub_quals = quals[s0:s1] if quals is not None else None
                qstr = ''.join(chr(q + 33) for q in sub_quals) if sub_quals is not None else ('~' * len(sub))
                fh.write(f"@{seq_id}_{a}_{b}\n")
                fh.write(str(sub) + "\n")
                fh.write("+\n")
                fh.write(qstr + "\n")
    else:
        with open(fname, 'w') as fh:
            for (a, b) in intervals:
                s0 = a - 1
                s1 = b
                sub = seq[s0:s1]
                header = f">{seq_id}_{a}_{b}"
                fh.write(header + "\n")
                fh.write(str(sub) + "\n")
    return str(fname)


def _write_record_handle(fh, seq_id, seq, quals=None, is_fastq=False):
    """Write a single FASTA/FASTQ record to an open text handle.

    fh: text-mode file-like object
    quals: list of ints or None
    """
    if is_fastq:
        qstr = ''.join(chr(q + 33) for q in quals) if quals is not None else ('~' * len(seq))
        fh.write(f"@{seq_id}\n")
        fh.write(str(seq) + "\n")
        fh.write("+\n")
        fh.write(qstr + "\n")
    else:
        fh.write(f">{seq_id}\n")
        fh.write(str(seq) + "\n")


def _write_fragments_handle(fh, seq_id, seq, intervals, quals=None, is_fastq=False):
    """Write multiple fragment records (for a single input record) to fh."""
    for (a, b) in intervals:
        s0 = a - 1
        s1 = b
        sub = seq[s0:s1]
        if is_fastq:
            sub_quals = quals[s0:s1] if quals is not None else None
            qstr = ''.join(chr(q + 33) for q in sub_quals) if sub_quals is not None else ('~' * len(sub))
            fh.write(f"@{seq_id}_{a}_{b}\n")
            fh.write(str(sub) + "\n")
            fh.write("+\n")
            fh.write(qstr + "\n")
        else:
            fh.write(f">{seq_id}_{a}_{b}\n")
            fh.write(str(sub) + "\n")


def _process_record(args):
    """Worker that writes fragments (or the original record) for a single record to a temp file.

    Returns a tuple (idx, tmp_path_str).
    """
    idx, rid, seq, quals, coords_dir, is_fastq, tmp_dir = args
    coords_dir = Path(coords_dir)
    tmp_dir = Path(tmp_dir)
    # use an index-prefixed filename so merge ordering is trivial
    safe_rid = rid.replace('/', '_')
    tmp_path = tmp_dir / f"{idx:06d}_{safe_rid}.part"
    with open(tmp_path, 'w') as fh:
        # find coords file for this record (match *{rid}*cut_coordinates.tsv)
        matches = list(coords_dir.glob(f"*{rid}*cut_coordinates.tsv"))
        if not matches:
            # write record unchanged
            _write_record_handle(fh, rid, seq, quals=quals, is_fastq=is_fastq)
        else:
            coords_path = matches[0]
            ivals = parse_coords(coords_path)
            if not ivals:
                _write_record_handle(fh, rid, seq, quals=quals, is_fastq=is_fastq)
            else:
                _write_fragments_handle(fh, rid, seq, ivals, quals=quals, is_fastq=is_fastq)
    return idx, str(tmp_path)


def write_removed(seq_id, seq, intervals, out_dir, prefix, quals=None, is_fastq=False):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    parts = []
    qual_parts = []
    last = 0
    n = len(seq)
    for (a, b) in intervals:
        s0 = a - 1
        if last < s0:
            parts.append(seq[last:s0])
            if quals is not None:
                qual_parts.append(quals[last:s0])
        last = b
    if last < n:
        parts.append(seq[last:])
        if quals is not None:
            qual_parts.append(quals[last:])
    newseq = ''.join(parts)
    # prefix may be tuple (base, ext) or base string
    if isinstance(prefix, tuple):
        base, ext = prefix
    else:
        base, ext = prefix, '.fasta'
    fname = out_dir / f"{base}.cut_removed{ext}"
    if is_fastq:
        # write FASTQ (ext may already be .fastq or .fastq.gz)
        fname = out_dir / f"{base}.cut_removed{ext}"
        # flatten qual_parts to list
        newqual = []
        for qpart in qual_parts:
            newqual.extend(qpart)
        qstr = ''.join(chr(q + 33) for q in newqual) if newqual else ('~' * len(newseq))
        with open(fname, 'w') as fh:
            fh.write(f"@{seq_id}.cut_removed len={len(newseq)}\n")
            fh.write(str(newseq) + "\n")
            fh.write("+\n")
            fh.write(qstr + "\n")
    else:
        with open(fname, 'w') as fh:
            fh.write(f">{seq_id}.cut_removed len={len(newseq)}\n")
            fh.write(str(newseq) + "\n")
    return str(fname)


def main():
    ap = argparse.ArgumentParser(description="Cut a FASTA by coords TSV")
    ap.add_argument('--input-type', choices=['fasta', 'fastq'], default='fasta', help='Input file type: fasta or fastq; determines output format')
    ap.add_argument('--input-file', required=True, help='Input FASTA/FASTQ file (may be multi-record, may be gzipped)')
    ap.add_argument('--coords-dir', required=True, help='Directory containing *cut_coordinates.tsv files')
    ap.add_argument('--out-dir', default='.', help='Output directory')
    ap.add_argument('-j', '--threads', type=int, default=1, help='Number of worker processes to use (default: 1)')
    # output naming: derived from input filename (see derive_output_name below)
    args = ap.parse_args()

    input_path = Path(args.input_file)
    if not input_path.exists() or not input_path.is_file():
        raise SystemExit(f"Input file not found: {input_path}")

    coords_dir = Path(args.coords_dir)
    if not coords_dir.exists() or not coords_dir.is_dir():
        raise SystemExit(f"Coords directory not found: {coords_dir}")

    # coords files are handled per-record by searching --coords-dir; no global
    # coords parse here.

    # For file-based processing we require Biopython
    if not _HAS_BIO:
        raise SystemExit('Biopython is required; please install biopython')

    # decide output extension: preserve gzippedness of input file
    in_name = input_path.name
    lower = in_name.lower()
    gz_in = lower.endswith('.gz')
    if args.input_type == 'fastq':
        out_ext = '.fastq.gz' if gz_in else '.fastq'
    else:
        out_ext = '.fasta.gz' if gz_in else '.fasta'

    # derive output path from input filename
    if lower.endswith('.fastq.gz'):
        base = in_name[:-len('.fastq.gz')]
    elif lower.endswith('.fasta.gz'):
        base = in_name[:-len('.fasta.gz')]
    else:
        base = input_path.stem
    file_base = f"{base}.collapsed"
    out_path = Path(args.out_dir) / f"{file_base}{out_ext}"

    # open input handle (gzip if needed) and collect records into memory
    if gz_in:
        in_fh = gzip.open(input_path, 'rt')
    else:
        in_fh = open(input_path, 'r')

    is_fastq = (args.input_type == 'fastq')
    fmt = 'fastq' if is_fastq else 'fasta'

    # collect records (we store plain strings/lists so they are picklable for workers)
    records = []
    try:
        for i, rec in enumerate(SeqIO.parse(in_fh, fmt)):
            rid = rec.id
            seq = str(rec.seq)
            quals = rec.letter_annotations.get('phred_quality') if is_fastq else None
            records.append((i, rid, seq, quals))
    finally:
        in_fh.close()

    # If user requested single-threaded behavior, keep original streaming-style writing
    if args.threads <= 1:
        # open output handle (gzip if needed)
        if out_ext.endswith('.gz'):
            out_fh = gzip.open(out_path, 'wt')
        else:
            out_fh = open(out_path, 'w')
        try:
            for (i, rid, seq, quals) in records:
                matches = list(coords_dir.glob(f"*{rid}*cut_coordinates.tsv"))
                if not matches:
                    _write_record_handle(out_fh, rid, seq, quals=quals, is_fastq=is_fastq)
                    continue
                coords_path = matches[0]
                ivals = parse_coords(coords_path)
                if not ivals:
                    _write_record_handle(out_fh, rid, seq, quals=quals, is_fastq=is_fastq)
                else:
                    _write_fragments_handle(out_fh, rid, seq, ivals, quals=quals, is_fastq=is_fastq)
        finally:
            out_fh.close()
    else:
        # Parallel mode: spawn worker pool to write per-record temp files, then merge them.
        tmp_dir = Path(args.out_dir) / f".parallel_tmp_{os.getpid()}"
        tmp_dir.mkdir(parents=True, exist_ok=True)

        worker_args = []
        for (i, rid, seq, quals) in records:
            worker_args.append((i, rid, seq, quals, str(coords_dir), is_fastq, str(tmp_dir)))

        pool = multiprocessing.Pool(processes=args.threads)
        try:
            results = pool.map(_process_record, worker_args)
        finally:
            pool.close()
            pool.join()

        # results: list of (idx, tmp_path). Merge in index order into final output
        results.sort(key=lambda x: x[0])
        if out_ext.endswith('.gz'):
            out_fh = gzip.open(out_path, 'wt')
        else:
            out_fh = open(out_path, 'w')
        try:
            for idx, tmp_path in results:
                with open(tmp_path, 'r') as fh:
                    shutil.copyfileobj(fh, out_fh)
        finally:
            out_fh.close()

        # cleanup temp files
        for idx, tmp_path in results:
            try:
                os.remove(tmp_path)
            except Exception:
                pass
        try:
            tmp_dir.rmdir()
        except Exception:
            pass

    print('Wrote combined fragments:')
    print('  ', out_path)

if __name__ == '__main__':
    main()
