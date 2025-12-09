#!/usr/bin/env python3
"""
Collapser.py — repeat finder (direct + inverted), optional dotplot,
RC-junction cut, and monomerization (automatic periods ≥ --min-period).

This file contains a compact implementation of several steps used when
analysing long DNA reads / sequences for internal repeats and periodicity.

High-level stages (executed in main):
    1) Read the first sequence from the supplied FASTA file
    2) Collect exact k-mer seeds (direct + reverse-complement) used to find repeats
    3) Optionally produce a self-dotplot PNG (requires matplotlib)
    4) Find long repeats using seed chaining and report counts

Optional/auxiliary steps that are kept but run later in the pipeline:
    • Optionally detect and write the coordinates of a dominant RC (reverse-complement) junction
    • Optionally monomerize the sequence by detecting periodicity and writing cut positions

Outputs (note: several outputs can be suppressed per user settings):
        • <prefix>.png                    — self-dotplot (if --dotplot and matplotlib installed)
        • <prefix>.rc_cut_coords.txt      — RC junction cut (1-based inclusive) if --rc-repeats
        • <prefix>.monomers.breaks.txt    — 1-based cut positions (union of all periods used)

Recommended for noisy reads (example):
    --k 9 --bin-size 150-300 --max-gap 2000-3000 --min-length 1000

"""

import argparse
from collections import defaultdict, Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from Bio import SeqIO
from Bio.Align import PairwiseAligner
# Use non-interactive backend for matplotlib in headless environments
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except Exception:
    plt = None
import os
import gzip
import concurrent.futures
import multiprocessing as mp
import sys

# Prefer the 'fork' start method on POSIX to avoid spawn-related
# pickling/import issues when worker callables are defined in `__main__`.
if sys.platform != 'win32':
    try:
        mp.set_start_method('fork')
    except RuntimeError:
        # start method already set; ignore
        pass

# per-worker aligner (initialized in child processes when using ProcessPoolExecutor)
_WORKER_ALIGNER = None

def _init_worker_aligner():
    """Initializer for worker processes: create a PairwiseAligner once per process."""
    global _WORKER_ALIGNER
    try:
        a = PairwiseAligner()
        a.mode = 'global'
        a.match_score = 1.0
        a.mismatch_score = 0.0
        a.open_gap_score = -1.0
        a.extend_gap_score = -0.5
        _WORKER_ALIGNER = a
    except Exception:
        _WORKER_ALIGNER = None


def compute_alignment_stats_worker(seqA: str, seqB: str):
    """Compute alignment-derived stats for two sequences (global align).

    Returns: (matches, mismatches, bases_aln, bases_not_aln, gaps, gap_opens, pct_id, pct_id_ungapped)
    """
    # Prefer a per-process preinitialized aligner to avoid constructing a new
    # PairwiseAligner object for every alignment (expensive). When using a
    # ProcessPoolExecutor with initializer `_init_worker_aligner` the module
    # global `_WORKER_ALIGNER` will be set in the child; fall back to creating one.
    aligner = _WORKER_ALIGNER
    if aligner is None:
        aligner = PairwiseAligner()
        aligner.mode = 'global'
        aligner.match_score = 1.0
        aligner.mismatch_score = 0.0
        aligner.open_gap_score = -1.0
        aligner.extend_gap_score = -0.5

    try:
        alns = aligner.align(seqA, seqB)
        aln = alns[0]
        a_blocks, b_blocks = aln.aligned
    except Exception:
        a_blocks = [(0, len(seqA))]
        b_blocks = [(0, len(seqB))]

    last_a = last_b = 0
    a_al = []
    b_al = []
    for (a_start, a_end), (b_start, b_end) in zip(a_blocks, b_blocks):
        if a_start > last_a:
            a_al.append(seqA[last_a:a_start])
            b_al.append('-' * (a_start - last_a))
        if b_start > last_b:
            a_al.append('-' * (b_start - last_b))
            b_al.append(seqB[last_b:b_start])
        a_al.append(seqA[a_start:a_end])
        b_al.append(seqB[b_start:b_end])
        last_a = a_end
        last_b = b_end
    if last_a < len(seqA):
        a_al.append(seqA[last_a:])
        b_al.append('-' * (len(seqA) - last_a))
    if last_b < len(seqB):
        a_al.append('-' * (len(seqB) - last_b))
        b_al.append(seqB[last_b:])

    a_al_s = ''.join(a_al)
    b_al_s = ''.join(b_al)
    aln_len = len(a_al_s)
    matches = mismatches = gaps = gap_opens = 0
    pct_id = pct_id_ungapped = 0.0
    if aln_len > 0:
        matches = sum(1 for x, y in zip(a_al_s, b_al_s) if x == y and x != '-')
        mismatches = sum(1 for x, y in zip(a_al_s, b_al_s) if x != y and x != '-' and y != '-')
        gaps = a_al_s.count('-') + b_al_s.count('-')
        gap_opens = sum(1 for i in range(aln_len) if a_al_s[i] == '-' and (i == 0 or a_al_s[i-1] != '-'))
        gap_opens += sum(1 for i in range(aln_len) if b_al_s[i] == '-' and (i == 0 or b_al_s[i-1] != '-'))
        pct_id = (matches / aln_len * 100.0) if aln_len > 0 else 0.0
        pct_id_ungapped = (matches / (matches + mismatches) * 100.0) if (matches + mismatches) > 0 else 0.0

    bases_aln = matches + mismatches
    bases_not_aln = gaps
    return (matches, mismatches, bases_aln, bases_not_aln, gaps, gap_opens, pct_id, pct_id_ungapped)


# Note: worker initializer `_init_worker_aligner` is defined above and sets
# the module-global `_WORKER_ALIGNER`. We avoid duplicating another
# initializer or `_GLOBAL_ALIGNER` to keep behavior explicit and simple.

# ---------- basic io ----------
_rc_map = str.maketrans("ACGTNacgtn", "TGCANtgcan")

def revcomp(s: str) -> str:
    """Return the reverse-complement of DNA sequence `s`.

    Example: revcomp('ACG') -> 'CGT'
    The translation map handles upper- and lower-case bases and 'U'->'T'
    is enforced later during sequence reading.
    """
    return s.translate(_rc_map)[::-1]

def read_first_fasta_sequence(path: str) -> str:
    for rec in SeqIO.parse(path, "fasta"):
        return str(rec.seq).upper().replace("U", "T")
    return ""

def read_first_fasta_header(path: str) -> Optional[str]:
    for rec in SeqIO.parse(path, "fasta"):
        return rec.id
    return None


def iter_reads_file(path: str):
    """Yield (rec_id, sequence) for all records in FASTA/FASTQ (supports .gz).

    Sequences are returned as upper-case strings with U->T normalization.
    """
    p = Path(path)
    # compressed (gz) files: open with gzip in text mode
    if str(path).lower().endswith('.gz'):
        # choose format heuristically from filename
        fname = str(path).lower()
        fmt = 'fastq' if (fname.endswith('.fastq.gz') or fname.endswith('.fq.gz')) else 'fasta'
        with gzip.open(path, 'rt') as fh:
            for rec in SeqIO.parse(fh, fmt):
                yield rec.id, str(rec.seq).upper().replace('U', 'T')
        return

    # uncompressed files: choose format by suffix
    ext = p.suffix.lower()
    fmt = 'fastq' if ext in ('.fastq', '.fq') else 'fasta'
    for rec in SeqIO.parse(path, fmt):
        yield rec.id, str(rec.seq).upper().replace('U', 'T')

# ---------- seeds & dotplot ----------
@dataclass
class Seed:
    i: int
    j: int
    ori: int  # +1 direct, -1 reverse-comp

def iter_kmers(s: str, k: int, step: int = 1):
    for i in range(0, len(s) - k + 1, step):
        yield i, s[i:i+k]
def collect_seeds(s: str, k: int, step: int = 1, include_self: bool = True) -> List[Seed]:
    """Collect exact k-mer self-matches (direct + RC).

    Returns a list of Seed(i,j,ori) where ori==+1 for direct matches and
    ori==-1 for reverse-complement matches.

    Parameters
    ----------
    s, k, step, include_self
    """
    idx: Dict[str, List[int]] = defaultdict(list)
    for i, kmer in iter_kmers(s, k, step):
        idx[kmer].append(i)
    seeds: List[Seed] = []
    # direct
    for i, kmer in iter_kmers(s, k, step):
        for j in idx.get(kmer, []):
            if include_self or i != j:
                seeds.append(Seed(i, j, +1))
    # reverse-comp
    for i, kmer in iter_kmers(s, k, step):
        for j in idx.get(revcomp(kmer), []):
            if i != j:
                seeds.append(Seed(i, j, -1))
    return seeds

def write_dotplot_png(seeds: List[Seed], n: int, out_png: str, title: str):
    if plt is None:
        return
    xs_f = [s.i for s in seeds if s.ori == +1]
    ys_f = [s.j for s in seeds if s.ori == +1]
    xs_r = [s.i for s in seeds if s.ori == -1]
    ys_r = [s.j for s in seeds if s.ori == -1]
    fig = plt.figure(figsize=(7.2, 7.2), dpi=160)
    ax = fig.add_subplot(111)
    # draw a light-grey grid behind the scatter points
    # ensure grid lines are drawn below plotted artists
    ax.set_axisbelow(True)
    ax.grid(color='lightgrey', linestyle='-', linewidth=0.5, alpha=0.7)

    # plot points above the grid (use zorder > grid default)
    ax.scatter(xs_f, ys_f, s=2, alpha=0.35, label="Direct", zorder=3)
    ax.scatter(xs_r, ys_r, s=2, alpha=0.35, label="Reverse-comp", zorder=3)
    ax.set_xlim(0, n); ax.set_ylim(0, n)
    ax.set_xlabel(title); ax.set_ylabel(title)
    ax.set_title(f"Self-dotplot {title}")
    ax.legend(loc="upper left", markerscale=3)
    fig.tight_layout(); fig.savefig(out_png); plt.close(fig)

# ---------- repeat finding (bin + chain) ----------
@dataclass
class Repeat:
    start1: int; end1: int; start2: int; end2: int; ori: str; length: int

def _segments_from_bin(points: List[Tuple[int,int]], k: int, max_gap: int) -> List[Tuple[int,int,int,int]]:
    """Chain sorted (i,j) into segments allowing gaps; return half-open coords."""
    if not points: return []
    points.sort(key=lambda x: x[0])
    segs: List[Tuple[int,int,int,int]] = []
    imin = imax = points[0][0]; jmin = jmax = points[0][1]
    prev_i, prev_j = points[0]
    for i, j in points[1:]:
        if (i - prev_i) <= max_gap and abs(j - prev_j) <= max_gap:
            imax = max(imax, i); jmin = min(jmin, j); jmax = max(jmax, j)
        else:
            segs.append((imin, imax + k, jmin, jmax + k))
            imin = imax = i; jmin = jmax = j
        prev_i, prev_j = i, j
    segs.append((imin, imax + k, jmin, jmax + k))
    return segs

def find_long_repeats(s: str, k: int, step: int, min_length: int = 1000, bin_size: int = 300, max_gap: int = 3000) -> List[Repeat]:
    """Find repeats (direct '+', inverted '-') >= min_length using exact seeds + chaining.

    Steps:
      1) index k-mers sampled by `step` -> position lists
      2) generate direct and RC seed pairs (i,j)
      3) bin pairs to coarse diagonals/anti-diagonals
      4) chain nearby points inside each bin into segments
      5) keep segments whose overlap length >= min_length and merge them
    """
    idx: Dict[str, List[int]] = defaultdict(list)
    for i in range(0, len(s) - k + 1, step):
        idx[s[i:i+k]].append(i)

    plus: List[Tuple[int,int]] = []
    minus: List[Tuple[int,int]] = []
    for i in range(0, len(s) - k + 1, step):
        kmer = s[i:i+k]
        for j in idx.get(kmer, []):
            if i != j:
                plus.append((i, j))
        for j in idx.get(revcomp(kmer), []):
            if i != j:
                minus.append((i, j))

    bins_plus: Dict[int, List[Tuple[int,int]]] = defaultdict(list)
    for i, j in plus:
        bins_plus[(i - j) // bin_size].append((i, j))

    bins_minus: Dict[int, List[Tuple[int,int]]] = defaultdict(list)
    for i, j in minus:
        bins_minus[(i + j) // bin_size].append((i, j))

    reps: List[Repeat] = []
    for pts in bins_plus.values():
        for a,b,c,d in _segments_from_bin(pts, k, max_gap):
            L = min(b-a, d-c)
            if L >= min_length:
                reps.append(Repeat(a,b,c,d,'+',L))
    for pts in bins_minus.values():
        for a,b,c,d in _segments_from_bin(pts, k, max_gap):
            L = min(b-a, d-c)
            if L >= min_length:
                reps.append(Repeat(a,b,c,d,'-',L))

    reps.sort(key=lambda r: (r.ori, r.start1, r.start2))
    merged: List[Repeat] = []
    for r in reps:
        if merged and r.ori == merged[-1].ori and r.start1 <= merged[-1].end1 and r.start2 <= merged[-1].end2:
            m = merged[-1]
            m.end1 = max(m.end1, r.end1)
            m.end2 = max(m.end2, r.end2)
            m.length = min(m.end1 - m.start1, m.end2 - m.start2)
        else:
            merged.append(r)
    return merged

# ---------- RC junction (1-based) ----------
def rc_repeats_coords(s: str, k: int, step: int, min_seeds: int = 1) -> tuple[int,int,int,int]:
    """Detect dominant reverse-complement (RC) anti-diagonal.

    If a sequence contains a strong RC-junction (for example a circular
    molecule concatenated in opposite orientation), then many k-mers will
    match to reverse-complemented positions and these pairs concentrate
    on an anti-diagonal in the self dotplot. This function finds the
    dominant anti-diagonal by simple binning and returns a 1-based
    inclusive coordinate split (start1,end1,start2,end2) describing the
    best place to cut the sequence into two parts (left/right).

    If no RC pairs are detected, the whole sequence range is returned.
    """
    n = len(s)
    idx = defaultdict(list)
    for i in range(0, n - k + 1, step):
        idx[s[i:i+k]].append(i)
    rc_pairs = []
    for i in range(0, n - k + 1, step):
        for j in idx.get(revcomp(s[i:i+k]), []):
            if i != j: rc_pairs.append((i, j))
    if not rc_pairs: return (1, n, 1, n)
    bin_size = max(5, k)
    hist = defaultdict(int)
    for i, j in rc_pairs: hist[(i + j)//bin_size] += 1
    peak = max(hist, key=hist.get)
    # require a minimum number of RC k-mer seed pairs in the peak to accept a split
    peak_count = hist.get(peak, 0)
    if peak_count < max(1, min_seeds):
        return (1, n, 1, n)
    cand = [(i,j) for i,j in rc_pairs if abs(((i+j)//bin_size)-peak) <= 1] or rc_pairs
    max_i = max(i for i,_ in cand); min_j = min(j for _,j in cand)
    t = (max_i + k + min_j) // 2  # 0-based split point
    return (1, t, t+1, n)

# ---------- Period detection & monomerization ----------
def estimate_all_periods_from_seeds(
    seeds: List[Seed],
    bin_w: int = 200,
    min_period: int = 1000,
    max_period: Optional[int] = None,
    harmonic_tol: float = 0.10,
    min_hits: int = 1000,
) -> List[int]:
    """Return all strong direct-repeat periods ≥ min_period, filtering near-harmonics.

    Compute |i-j| for direct seeds, bin them by bin_w, and return peak periods
    that have at least min_hits counts and are not harmonics of already
    selected periods (within harmonic_tol).
    """
    diffs = [abs(s.i - s.j) for s in seeds if s.ori == +1 and s.i != s.j]
    if not diffs:
        return []
    if max_period is None:
        max_period = max(diffs)
    binned = [d // bin_w for d in diffs if min_period <= d <= max_period]
    if not binned:
        return []
    counts = Counter(binned)
    counts.pop(0, None)
    peaks = sorted(((b*bin_w, c) for b,c in counts.items() if c >= min_hits), key=lambda x: x[0])
    selected: List[int] = []
    for p,_ in peaks:
        keep = True
        for q in selected:
            r = p / q
            nearest = round(r)
            if nearest >= 2 and abs(r - nearest) <= harmonic_tol:
                keep = False
                break
        if keep:
            selected.append(p)
    return selected
# Monomerization helpers removed per user request (simplified CLI)

# ---------- writing ----------
def write_seeds_txt(out_prefix: str, seeds: List[Seed], fasta_path: Optional[str] = None, read_id: Optional[str] = None):
    """Write k-mer match coordinates (direct and RC) to a TSV file.

    The requested filename is based on the input FASTA stem: <stem>.kmerMatchCordinates.tsv
    If `fasta_path` is not provided, fall back to using <out_prefix>.kmerMatchCordinates.tsv.

    Output columns (tab-separated, 1-based inclusive):
      i	j	ori	kmer
    where `ori` is '+' for direct matches and '-' for reverse-complement matches.
    """
    if not seeds:
        return None
    # determine filename: prefer fasta stem if provided
    # prefer to place the TSV in the same directory as out_prefix (this follows --out-dir)
    out_dir = Path(out_prefix).parent if out_prefix else Path('.')
    out_dir.mkdir(parents=True, exist_ok=True)
    if fasta_path:
        try:
            # preserve the original filename (including extensions) so
            # output looks like <input.fastq>.kmerMatchCordinates.tsv
            fname = Path(fasta_path).name
            if read_id:
                path = str(out_dir / (f"{fname}.{read_id}.kmerMatchCordinates.tsv"))
            else:
                path = str(out_dir / (fname + ".kmerMatchCordinates.tsv"))
        except Exception:
            path = str(Path(out_prefix) / (Path(out_prefix).name + ".kmerMatchCordinates.tsv"))
    else:
        path = str(Path(out_prefix) / (Path(out_prefix).name + ".kmerMatchCordinates.tsv"))

    # write TSV with 1-based coordinates
    with open(path, "w") as fh:
        fh.write("#i\tj\tori\n")
        # Note: we don't have kmer string available here by default; reconstruct from seed positions
        # Attempt to infer k from the minimum positive difference between i and j where ori==+1
        # Fallback: if unable to infer, leave kmer field empty.
        # Try to infer k by looking at repeated seed lengths from nearby seeds with same i and ori
        for s in seeds:
            i1 = s.i + 1
            j1 = s.j + 1
            ori = '+' if s.ori == +1 else '-'
            fh.write(f"{i1}\t{j1}\t{ori}\t\n")
    return path

def write_repeats(out_prefix: str, repeats: List[Repeat], s: Optional[str] = None, min_pct_id: float = 0.0, jobs: int = 1):
    # Write repeats TSVs (filtered + unfiltered) with alignment stats when
    # sequence `s` is provided. Returns (filtered_path, unfiltered_path).
    # Always create the filtered and unfiltered files (headers at minimum)
    # so downstream pipelines don't fail when there are no repeats.
    path_filtered = str(out_prefix) + ".repeats.filtered.tsv"
    path_unfiltered = str(out_prefix) + ".repeats.unfiltered.tsv"
    header = "#start1\tend1\tstart2\tend2\tori\tlength"
    stats_header = "\tmatches\tmismatches\tbases_aln\tbases_not_aln\tgaps\tgap_opens\tpct_id\tpct_id_ungapped\n"
    filtered_lines = [header + stats_header]
    unfiltered_lines = [header + stats_header]

    # Precompute reverse-complement of the full sequence to avoid per-repeat revcomp
    n = len(s) if s is not None else 0
    s_rc = revcomp(s) if s is not None else None

    # If there are no repeats, write header-only files and return paths
    path_filtered = str(out_prefix) + ".repeats.filtered.tsv"
    path_unfiltered = str(out_prefix) + ".repeats.unfiltered.tsv"
    legacy_path = str(out_prefix) + ".repeats.chained.tsv"
    if not repeats:
        # No repeats found: do not write header-only repeat files to avoid
        # cluttering output with empty placeholders.
        return (None, None)

    all_entries = []
    # collect tasks for alignment/stat computation when sequence `s` is provided
    tasks_r = []
    tasks_seqA = []
    tasks_seqB = []
    for r in repeats:
        if r.start1 == r.start2 and r.end1 == r.end2:
            continue
        a = r.start1 + 1
        b = r.end1
        c = r.start2 + 1
        d = r.end2
        L = r.length

        if s is None:
            line = f"{a}\t{b}\t{c}\t{d}\t{r.ori}\t{L}\n"
            unfiltered_lines.append(line)
            all_entries.append((r.start1, r.end1, L, line))
            if 0.0 >= min_pct_id:
                filtered_lines.append(line)
            continue

        seqA = s[r.start1:r.start1 + L]
        if r.ori == '-':
            # reverse-complement of s[start2:start2+L] == s_rc[n-(start2+L) : n-start2]
            seqB = s_rc[n - (r.start2 + L) : n - r.start2]
        else:
            seqB = s[r.start2:r.start2 + L]
        # queue task for later stat computation
        tasks_r.append(r)
        tasks_seqA.append(seqA)
        tasks_seqB.append(seqB)

    # compute stats either serially or in parallel
    # number of worker processes passed in by caller (CLI --jobs)
    jobs = max(1, int(jobs or 1))
    if tasks_r:
        if jobs > 1:
            # choose a reasonable chunksize to reduce IPC overhead; ensure at least 1
            try:
                chunksize = max(1, len(tasks_seqA) // (jobs * 4))
            except Exception:
                chunksize = 16
            # Use a 'fork' start method to avoid spawn-related pickling issues
            # when submitting top-level callables defined in this script.
            # This is safe on POSIX systems (Linux) used in the cluster.
            mp_ctx = mp.get_context('fork')
            with concurrent.futures.ProcessPoolExecutor(max_workers=jobs, initializer=_init_worker_aligner, mp_context=mp_ctx) as exc:
                for r, (matches, mismatches, bases_aln, bases_not_aln, gaps, gap_opens, pct_id, pct_id_ungapped) in zip(tasks_r, exc.map(compute_alignment_stats_worker, tasks_seqA, tasks_seqB, chunksize)):
                    a = r.start1 + 1
                    b = r.end1
                    c = r.start2 + 1
                    d = r.end2
                    L = r.length
                    line = f"{a}\t{b}\t{c}\t{d}\t{r.ori}\t{L}\t{matches}\t{mismatches}\t{bases_aln}\t{bases_not_aln}\t{gaps}\t{gap_opens}\t{pct_id:.2f}\t{pct_id_ungapped:.2f}\n"
                    unfiltered_lines.append(line)
                    all_entries.append((r.start1, r.end1, L, line))
                    if pct_id >= min_pct_id:
                        filtered_lines.append(line)
        else:
            for r, seqA, seqB in zip(tasks_r, tasks_seqA, tasks_seqB):
                matches, mismatches, bases_aln, bases_not_aln, gaps, gap_opens, pct_id, pct_id_ungapped = compute_alignment_stats_worker(seqA, seqB)
                a = r.start1 + 1
                b = r.end1
                c = r.start2 + 1
                d = r.end2
                L = r.length
                line = f"{a}\t{b}\t{c}\t{d}\t{r.ori}\t{L}\t{matches}\t{mismatches}\t{bases_aln}\t{bases_not_aln}\t{gaps}\t{gap_opens}\t{pct_id:.2f}\t{pct_id_ungapped:.2f}\n"
                unfiltered_lines.append(line)
                all_entries.append((r.start1, r.end1, L, line))
                if pct_id >= min_pct_id:
                    filtered_lines.append(line)

    wrote_any = False
    written_filtered = False
    written_unfiltered = False
    if len(filtered_lines) > 1:
        try:
            with open(path_filtered, "w") as fh:
                fh.writelines(filtered_lines)
            written_filtered = True
            wrote_any = True
        except Exception:
            written_filtered = False
    if len(unfiltered_lines) > 1:
        try:
            with open(path_unfiltered, "w") as fu:
                fu.writelines(unfiltered_lines)
            written_unfiltered = True
            wrote_any = True
        except Exception:
            written_unfiltered = False

    # build legacy aggregated chained-fragment file: select dominating
    # (longest) non-overlapping left-blocks greedily so legacy file contains
    # one or more representative chained fragments depending on cutoffs
    legacy_path = str(out_prefix) + ".repeats.chained.tsv"
    legacy_lines = [header + stats_header]
    try:
        # sort entries by length descending and pick non-overlapping on left block
        all_entries_sorted = sorted(all_entries, key=lambda x: x[2], reverse=True)
        chosen_intervals: List[Tuple[int,int]] = []
        for start1, end1, length, line in all_entries_sorted:
            a1 = start1 + 1
            b1 = end1
            # check overlap with already chosen left blocks
            overlap = False
            for ca, cb in chosen_intervals:
                if not (b1 < ca or a1 > cb):
                    overlap = True
                    break
            if overlap:
                continue
            chosen_intervals.append((a1, b1))
            legacy_lines.append(line)
        # if nothing chosen, still write header-only file
        with open(legacy_path, "w") as lf:
            lf.writelines(legacy_lines)
    except Exception:
        # fallback: ensure a header-only file exists
        try:
            with open(legacy_path, "w") as lf:
                lf.writelines([header + stats_header])
        except Exception:
            pass

    return (path_filtered if written_filtered else None, path_unfiltered if written_unfiltered else None)


def write_trimmed_repeats(
    out_prefix: str,
    repeats: List[Repeat],
    s: str,
    min_pct_id: float = 0.0,
    suffix: str = '',
    ori_filter: Optional[str] = None,
):
    """Write tight (trimmed) repeat coordinates and alignment stats.

    For each Repeat r (0-based, half-open stored), the true aligned block
    length is r.length. We write 1-based inclusive coordinates along with
    computed alignment statistics using Bio.Align.PairwiseAligner:

    start1 end1 start2 end2 ori length matches mismatches gaps gap_opens pct_id pct_id_ungapped
    """
    # optionally filter repeats by orientation ('+' or '-')
    # Always create filtered/unfiltered files (header at minimum) so downstream
    # pipelines won't fail when there are no trimmed repeats.
    path = str(out_prefix) + f".repeats{suffix}.filtered.tsv"
    unfiltered_path = str(out_prefix) + f".repeats{suffix}.unfiltered.tsv"
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    # scoring tuned for counting matches/mismatches rather than scoring
    aligner.match_score = 1.0
    aligner.mismatch_score = 0.0
    aligner.open_gap_score = -1.0
    aligner.extend_gap_score = -0.5

    # Prepare buffered output lists
    header = "#start1\tend1\tstart2\tend2\tori\tlength\tmatches\tmismatches\tbases_aln\tbases_not_aln\tgaps\tgap_opens\tpct_id\tpct_id_ungapped\n"
    unfiltered_lines = [header]
    filtered_lines = [header]
    # track covered positions on left block (1-based intervals) to avoid overlaps
    covered_until_left = 0
    # process repeats in order of left block start (apply orientation filter earlier)
    # precompute reverse-complement once for efficiency
    n = len(s)
    s_rc = revcomp(s)
    if ori_filter in ('+', '-'):
        repeats_sorted = sorted([r for r in repeats if r.ori == ori_filter], key=lambda x: (x.start1, x.start2))
    else:
        repeats_sorted = sorted(repeats, key=lambda x: (x.start1, x.start2))
    for r in repeats_sorted:
        # skip exact self-matches where the two blocks are identical coordinates
        if r.start1 == r.start2 and r.end1 == r.end2:
            continue
        # if this repeat's left block is fully <= covered region, skip
        left_a = r.start1 + 1
        left_b = r.start1 + r.length
        if left_b <= covered_until_left:
            continue
        # if partially overlapping, trim the left block to start just after covered_until_left
        trim = 0
        if left_a <= covered_until_left < left_b:
            trim = covered_until_left - left_a + 1
        # compute new coordinates after trimming `trim` bases from left start
        a0 = r.start1 + trim
        L = r.length - trim
        if L <= 0:
            continue
        c0 = r.start2 + trim
        # now seqA/seqB reflect the possibly-trimmed blocks
        # seqA starts at a0 and has length L
        seqA = s[a0:a0+L]
        # for inverted repeats (ori == '-') the matched block should be
        # reverse-complemented before alignment so that aligner compares
        # sequences in the same orientation
        if r.ori == '-':
            # use precomputed reverse-complement slice
            seqB = s_rc[n - (c0 + L) : n - c0]
        else:
            seqB = s[c0:c0+L]
        matches = mismatches = gaps = gap_opens = 0
        pct_id = pct_id_ungapped = 0.0

        try:
            alns = aligner.align(seqA, seqB)
            aln = alns[0]
            aligned = aln.aligned  # tuple (a_blocks, b_blocks)
            a_blocks, b_blocks = aligned
        except Exception:
            # fallback: treat as ungapped comparison of trimmed length
            a_blocks = [(0, len(seqA))]
            b_blocks = [(0, len(seqB))]

        # rebuild aligned strings with gaps from aligned blocks
        last_a = last_b = 0
        a_al = []
        b_al = []
        for (a_start, a_end), (b_start, b_end) in zip(a_blocks, b_blocks):
            if a_start > last_a:
                a_al.append(seqA[last_a:a_start])
                b_al.append('-' * (a_start - last_a))
            if b_start > last_b:
                a_al.append('-' * (b_start - last_b))
                b_al.append(seqB[last_b:b_start])
            a_al.append(seqA[a_start:a_end])
            b_al.append(seqB[b_start:b_end])
            last_a = a_end
            last_b = b_end
        # tails
        if last_a < len(seqA):
            a_al.append(seqA[last_a:])
            b_al.append('-' * (len(seqA) - last_a))
        if last_b < len(seqB):
            a_al.append('-' * (len(seqB) - last_b))
            b_al.append(seqB[last_b:])

        a_al_s = ''.join(a_al)
        b_al_s = ''.join(b_al)
        aln_len = len(a_al_s)
        if aln_len > 0:
            matches = sum(1 for x, y in zip(a_al_s, b_al_s) if x == y and x != '-')
            mismatches = sum(1 for x, y in zip(a_al_s, b_al_s) if x != y and x != '-' and y != '-')
            gaps = a_al_s.count('-') + b_al_s.count('-')
            # count gap opens in both sequences
            gap_opens = sum(1 for i in range(aln_len) if a_al_s[i] == '-' and (i == 0 or a_al_s[i-1] != '-'))
            gap_opens += sum(1 for i in range(aln_len) if b_al_s[i] == '-' and (i == 0 or b_al_s[i-1] != '-'))
            pct_id = (matches / aln_len * 100.0) if aln_len > 0 else 0.0
            pct_id_ungapped = (matches / (matches + mismatches) * 100.0) if (matches + mismatches) > 0 else 0.0

        a = a0 + 1
        b = a0 + L
        c = c0 + 1
        d = c0 + L
        bases_aln = matches + mismatches
        bases_not_aln = gaps
        line = f"{a}\t{b}\t{c}\t{d}\t{r.ori}\t{L}\t{matches}\t{mismatches}\t{bases_aln}\t{bases_not_aln}\t{gaps}\t{gap_opens}\t{pct_id:.2f}\t{pct_id_ungapped:.2f}\n"
        # buffer unfiltered and filtered lines
        unfiltered_lines.append(line)
        if pct_id >= min_pct_id:
            filtered_lines.append(line)
            # update covered region to the end of this written left block
            covered_until_left = b

    wrote_any = False
    written_filtered = False
    written_unfiltered = False
    # write filtered file only if we have data rows beyond the header
    try:
        if len(filtered_lines) > 1:
            with open(path, "w") as fh:
                fh.writelines(filtered_lines)
            written_filtered = True
            wrote_any = True
        else:
            written_filtered = False
    except Exception:
        written_filtered = False

    # write unfiltered file only if we have data rows beyond the header
    try:
        if len(unfiltered_lines) > 1:
            with open(unfiltered_path, "w") as fu:
                fu.writelines(unfiltered_lines)
            written_unfiltered = True
            wrote_any = True
        else:
            written_unfiltered = False
    except Exception:
        written_unfiltered = False

    if not wrote_any:
        # nothing useful written
        return (None, None)
    return (path if written_filtered else None, unfiltered_path if written_unfiltered else None)


def _read_chained_intervals(path: str, want_ori: Optional[str] = None) -> List[Tuple[int,int]]:
    """Read a chained repeats TSV and return a list of 1-based inclusive tight intervals.

    The chained TSV contains bounding coordinates (start1,end1,start2,end2,ori,length,...).
    For cutting we prefer the tight aligned block: start..(start+length-1) for each side.
    This function extracts `start1`, `start2` and `length` when available and returns
    both left and right tight intervals. If parsing fails it falls back to using the
    provided start/end columns.

    - path: path to a chained TSV (header line starting with '#')
    - want_ori: if '+' or '-', only include rows matching that orientation; if None include all
    """
    intervals: List[Tuple[int,int]] = []
    try:
        with open(path, "r") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                # need at least start1 and start2 and ideally length
                if len(parts) < 4:
                    continue
                ori = None
                try:
                    if len(parts) >= 6:
                        ori = parts[4]
                except Exception:
                    ori = None
                if want_ori is not None and ori != want_ori:
                    continue

                try:
                    s1 = int(parts[0])
                    e1 = int(parts[1])
                    s2 = int(parts[2])
                    e2 = int(parts[3])
                except Exception:
                    # malformed coordinates
                    continue

                # prefer to use reported length (6th column) to compute tight blocks
                length = None
                if len(parts) >= 6:
                    try:
                        length = int(parts[5])
                    except Exception:
                        length = None

                if length and length > 0 and ori != '-':
                    # For normal chained rows prefer the reported tight `length`
                    # (start..start+length-1). However, for RC-derived chained
                    # entries the `length` field may represent the left-block
                    # extent (split point) and can span a very large region
                    # (covering most of the sequence). Using the bounding
                    # start/end columns for RC rows is safer for producing
                    # cut coordinates that highlight internal gaps rather
                    # than covering the whole sequence.
                    intervals.append((s1, s1 + length - 1))
                    intervals.append((s2, s2 + length - 1))
                else:
                    # fallback: use bounding start/end columns
                    intervals.append((s1, e1))
                    intervals.append((s2, e2))
    except Exception:
        return []
    return intervals


def _merge_intervals(intervals: List[Tuple[int,int]]) -> List[Tuple[int,int]]:
    """Merge overlapping or adjacent 1-based inclusive intervals.

    Returns a sorted list of non-overlapping intervals.
    """
    if not intervals:
        return []
    intervals_sorted = sorted(intervals, key=lambda x: x[0])
    merged: List[Tuple[int,int]] = []
    cur_a, cur_b = intervals_sorted[0]
    for a, b in intervals_sorted[1:]:
        # Only merge when intervals overlap. Do NOT merge adjacent
        # intervals (where a == cur_b + 1) because adjacent tight
        # blocks (e.g., left/right parts of an RC split) should remain
        # separate for cutting.
        if a <= cur_b:
            cur_b = max(cur_b, b)
        else:
            merged.append((cur_a, cur_b))
            cur_a, cur_b = a, b
    merged.append((cur_a, cur_b))
    return merged


def write_repeats_cut_coordinates(out_prefix: str, seq_len: Optional[int] = None) -> Optional[str]:
    """Create <out_prefix>.repeats.cut_coordinates.tsv by reading chained TSVs.

    This function looks for two possible chained files:
      - <out_prefix>.repeats.chained.tsv      (use '+' hits)
      - <out_prefix>.repeats.RC.chained.tsv   (use '-' hits)

    It collects intervals from the requested orientations, merges overlaps,
    and writes a TSV with header `#start\tend` (1-based inclusive). Returns
    the path written or None on failure.
    """
    p_plus = str(out_prefix) + ".repeats.chained.tsv"
    p_minus = str(out_prefix) + ".repeats.RC.chained.tsv"
    # collect intervals separately for '+' and '-' chained files
    plus_intervals: List[Tuple[int,int]] = []
    minus_intervals: List[Tuple[int,int]] = []
    if Path(p_plus).exists():
        plus_intervals = _read_chained_intervals(p_plus, want_ori='+') or []
    if Path(p_minus).exists():
        minus_intervals = _read_chained_intervals(p_minus, want_ori='-') or []

    # Prefer using direct-repeat (plus) intervals when available. RC-derived
    # intervals tend to be less precise (can span large split regions) and can
    # accidentally merge separate plus-intervals; using plus intervals when
    # present preserves the explicit left/gap/right structure expected by
    # downstream cut logic.
    if plus_intervals:
        intervals = plus_intervals
    else:
        intervals = minus_intervals

    # if no chained files found, nothing to do: do not write header-only file
    out_cut = str(out_prefix) + ".repeats.cut_coordinates.tsv"
    try:
        merged = _merge_intervals(intervals)
        if not merged:
            # No chained repeat intervals found -> skip creating a file
            return None

        with open(out_cut, "w") as fh:
            fh.write("#start\tend\n")
            # If seq_len is provided, include leading/trailing complementary
            # regions so the output covers 1..seq_len when possible.
            if seq_len is not None:
                # leading gap
                first_a, first_b = merged[0]
                if first_a > 1:
                    fh.write(f"1\t{first_a - 1}\n")
                # write merged intervals and internal gaps
                for i, (a, b) in enumerate(merged):
                    fh.write(f"{a}\t{b}\n")
                    if i + 1 < len(merged):
                        na, nb = merged[i+1]
                        gap_start = b + 1
                        gap_end = na - 1
                        if gap_start <= gap_end:
                            fh.write(f"{gap_start}\t{gap_end}\n")
                # trailing gap
                last_a, last_b = merged[-1]
                if last_b < seq_len:
                    fh.write(f"{last_b + 1}\t{seq_len}\n")
            else:
                # original behavior: write merged intervals + internal gaps
                for i, (a, b) in enumerate(merged):
                    fh.write(f"{a}\t{b}\n")
                    if i + 1 < len(merged):
                        na, nb = merged[i+1]
                        gap_start = b + 1
                        gap_end = na - 1
                        if gap_start <= gap_end:
                            fh.write(f"{gap_start}\t{gap_end}\n")
        return out_cut
    except Exception:
        return None


def process_single_read(rec_id: str, s: str, out_prefix_path_str: str, fasta_path: str, params: dict):
    """Process a single read: collect seeds, find repeats, write outputs.

    Returns: tuple(rec_id, wrote_any: bool, messages: List[str])
    """
    messages = []
    if not s:
        messages.append(f"Skipping empty read {rec_id}")
        return (rec_id, False, messages)
    title = rec_id
    n = len(s)
    messages.append(f"Processing read '{title}' length={n} bp from {fasta_path}")

    # Collect seeds
    step = 1
    k = int(params.get('k', 13))
    messages.append(f"  Collecting seeds with k={k}, step={step} ...")
    seeds = collect_seeds(s, k, step, include_self=True)
    messages.append(f"  Collected {len(seeds)} seeds (direct + RC).")

    # Find repeats
    reps = find_long_repeats(s, k, step, min_length=int(params.get('min_length', 2000)),
                             bin_size=int(params.get('bin_size', 300)), max_gap=int(params.get('max_gap', 1000)))
    plus_count = sum(1 for r in reps if r.ori == '+')
    minus_count = sum(1 for r in reps if r.ori == '-')
    messages.append(f"  Detected repeats for {rec_id}: +={plus_count}, -={minus_count} (total={len(reps)})")

    if not reps:
        messages.append(f"  No repeats to write for {rec_id}; skipping outputs.")
        return (rec_id, False, messages)

    safe_id = rec_id.replace('/', '_').replace('\\', '_')
    out_pref_read = str(out_prefix_path_str) + f".{safe_id}"

    # write seeds TSV
    try:
        _ = write_seeds_txt(out_pref_read, seeds, fasta_path=fasta_path, read_id=safe_id)
    except Exception:
        pass

    # optional dotplot
    if params.get('dotplot', False):
        try:
            write_dotplot_png(seeds, n, out_pref_read + ".png", title)
        except Exception:
            pass

    # Avoid nested multiprocessing oversubscription: if we're running multiple
    # per-read worker processes, force per-repeat alignment jobs to 1 unless
    # the user explicitly requested otherwise and is aware.
    local_jobs = int(params.get('jobs', 1))
    if int(params.get('workers', 1)) > 1 and local_jobs > 1:
        local_jobs = 1

    # write repeats
    try:
        repeats_res = write_repeats(out_pref_read, reps, s, min_pct_id=float(params.get('min_pct_id', 0.0)), jobs=local_jobs)
        if repeats_res:
            try:
                filtered_rep, unfiltered_rep = repeats_res
                if filtered_rep:
                    messages.append(f"  Wrote repeats (filtered) for {rec_id}: {filtered_rep}")
                if unfiltered_rep:
                    messages.append(f"  Wrote repeats (unfiltered) for {rec_id}: {unfiltered_rep}")
            except Exception:
                messages.append(f"  Wrote repeats summary for {rec_id}: {repeats_res}")
    except Exception:
        messages.append(f"  Failed to write repeats for {rec_id}")

    # Optionally detect RC junction and write RC split per-read
    if params.get('rc_repeats', False):
        try:
            min_seeds = max(1, int(params.get('min_length', 2000)) // max(1, k))
            s1, e1, s2, e2 = rc_repeats_coords(s, k, step, min_seeds=min_seeds)
            path_filtered = out_pref_read + ".repeats.RC.filtered.tsv"
            path_unfiltered = out_pref_read + ".repeats.RC.unfiltered.tsv"
            if not (s1 == 1 and e1 == n and s2 == 1 and e2 == n):
                pf = Path(path_filtered)
                pu = Path(path_unfiltered)
                write_header_filtered = (not pf.exists()) or (pf.exists() and pf.stat().st_size == 0)
                write_header_unfiltered = (not pu.exists()) or (pu.exists() and pu.stat().st_size == 0)
                left_len = e1 - s1 + 1
                seqA = s[s1-1:e1]
                seqB = revcomp(s[s2-1:e2])
                try:
                    aligner_rc = PairwiseAligner()
                    aligner_rc.mode = 'global'
                    aligner_rc.match_score = 1.0
                    aligner_rc.mismatch_score = 0.0
                    aligner_rc.open_gap_score = -1.0
                    aligner_rc.extend_gap_score = -0.5
                    alns = aligner_rc.align(seqA, seqB)
                    aln = alns[0]
                    a_blocks, b_blocks = aln.aligned
                except Exception:
                    a_blocks = [(0, len(seqA))]
                    b_blocks = [(0, len(seqB))]

                last_a = last_b = 0
                a_al = []
                b_al = []
                for (a_start, a_end), (b_start, b_end) in zip(a_blocks, b_blocks):
                    if a_start > last_a:
                        a_al.append(seqA[last_a:a_start])
                        b_al.append('-' * (a_start - last_a))
                    if b_start > last_b:
                        a_al.append('-' * (b_start - last_b))
                        b_al.append(seqB[last_b:b_start])
                    a_al.append(seqA[a_start:a_end])
                    b_al.append(seqB[b_start:b_end])
                    last_a = a_end
                    last_b = b_end
                if last_a < len(seqA):
                    a_al.append(seqA[last_a:])
                    b_al.append('-' * (len(seqA) - last_a))
                if last_b < len(seqB):
                    a_al.append('-' * (len(seqB) - last_b))
                    b_al.append(seqB[last_b:])

                a_al_s = ''.join(a_al)
                b_al_s = ''.join(b_al)
                aln_len = len(a_al_s)
                matches = mismatches = gaps = gap_opens = 0
                pct_id = pct_id_ungapped = 0.0
                if aln_len > 0:
                    matches = sum(1 for x, y in zip(a_al_s, b_al_s) if x == y and x != '-')
                    mismatches = sum(1 for x, y in zip(a_al_s, b_al_s) if x != y and x != '-' and y != '-')
                    gaps = a_al_s.count('-') + b_al_s.count('-')
                    gap_opens = sum(1 for i in range(aln_len) if a_al_s[i] == '-' and (i == 0 or a_al_s[i-1] != '-'))
                    gap_opens += sum(1 for i in range(aln_len) if b_al_s[i] == '-' and (i == 0 or b_al_s[i-1] != '-'))
                    pct_id = (matches / aln_len * 100.0) if aln_len > 0 else 0.0
                    pct_id_ungapped = (matches / (matches + mismatches) * 100.0) if (matches + mismatches) > 0 else 0.0

                line = f"{s1}\t{e1}\t{s2}\t{e2}\t-\t{left_len}\t{matches}\t{mismatches}\t{gaps}\t{gap_opens}\t{pct_id:.2f}\t{pct_id_ungapped:.2f}\n"
                with open(path_unfiltered, "a") as fh_un:
                    if write_header_unfiltered:
                        fh_un.write("#start1\tend1\tstart2\tend2\tori\tlength\tmatches\tmismatches\tgaps\tgap_opens\tpct_id\tpct_id_ungapped\n")
                    fh_un.write(line)
                try:
                    pct_thr = float(params.get('min_pct_id', 0.0))
                except Exception:
                    pct_thr = 0.0
                if pct_id >= pct_thr:
                    with open(path_filtered, "a") as fh_f:
                        if write_header_filtered:
                            fh_f.write("#start1\tend1\tstart2\tend2\tori\tlength\tmatches\tmismatches\tgaps\tgap_opens\tpct_id\tpct_id_ungapped\n")
                        fh_f.write(line)
                # write legacy RC chained file per-read
                legacy_rc_path = out_pref_read + ".repeats.RC.chained.tsv"
                try:
                    legacy_p = Path(legacy_rc_path)
                    write_legacy_header = (not legacy_p.exists()) or (legacy_p.exists() and legacy_p.stat().st_size == 0)
                    with open(legacy_rc_path, "a") as lf:
                        if write_legacy_header:
                            lf.write("#start1\tend1\tstart2\tend2\tori\tlength\tmatches\tmismatches\tgaps\tgap_opens\tpct_id\tpct_id_ungapped\n")
                        lf.write(line)
                except Exception:
                    pass
                messages.append(f"  [revComplRepeats] Suggested RC split for {rec_id}: {s1} {e1} {s2} {e2}")
            else:
                messages.append(f"  [revComplRepeats] No dominant RC junction detected for {rec_id}; no RC files written.")
        except Exception:
            messages.append(f"  Failed RC detection for {rec_id}")

    # Optionally write trimmed repeats
    if params.get('direct_repeats', False) or params.get('write_trimmed_repeats', False):
        try:
            res = write_trimmed_repeats(out_pref_read, reps, s, min_pct_id=float(params.get('min_pct_id', 0.0)), suffix='')
            if res:
                try:
                    filtered, unfiltered = res
                    messages.append(f"  Wrote trimmed repeats (filtered) for {rec_id}: {filtered}")
                    messages.append(f"  Wrote trimmed repeats (unfiltered) for {rec_id}: {unfiltered}")
                except Exception:
                    messages.append(f"  Wrote trimmed repeats for {rec_id}: {res}")
        except Exception:
            messages.append(f"  Failed writing trimmed repeats for {rec_id}")
        try:
            res_rc = write_trimmed_repeats(out_pref_read, reps, s, min_pct_id=float(params.get('min_pct_id', 0.0)), suffix='.RC', ori_filter='-')
            if res_rc:
                try:
                    filtered_rc, unfiltered_rc = res_rc
                    messages.append(f"  Wrote RC trimmed repeats (filtered) for {rec_id}: {filtered_rc}")
                    messages.append(f"  Wrote RC trimmed repeats (unfiltered) for {rec_id}: {unfiltered_rc}")
                except Exception:
                    messages.append(f"  Wrote RC trimmed repeats for {rec_id}: {res_rc}")
        except Exception:
            messages.append(f"  Failed writing RC trimmed repeats for {rec_id}")

    # build per-read cut coordinates TSV
    try:
        cut_path = write_repeats_cut_coordinates(out_pref_read, seq_len=n)
        if cut_path:
            messages.append(f"  Wrote repeats cut coordinates for {rec_id}: {cut_path}")
        else:
            messages.append(f"  No repeats cut coordinates written for {rec_id}")
    except Exception:
        messages.append(f"  Failed writing cut coordinates for {rec_id}")

    return (rec_id, True, messages)

# ---------- CLI ----------
def main():
    ap = argparse.ArgumentParser(
        description="Repeat finder (direct + inverted) with optional dotplot and auto monomerization",
        formatter_class=argparse.RawTextHelpFormatter
    )
    ap.add_argument("fasta"); ap.add_argument("out_prefix")
    ap.add_argument("--out-dir", type=str, default=None)
    
    ap.add_argument("--k", type=int, default=13, help="k-mer size for exact seeds (direct + RC)")
    ap.add_argument("--min-length", type=int, default=2000, help="minimum length for k-mer chains to report as repeats")
    ap.add_argument("--bin-size", type=int, default=300, 
        help=("coarse bin size in bp used to group k-mer seed pairs into "
              "diagonals/anti-diagonals (used for coarse binning before "
              "chaining). Smaller values = finer resolution (higher "
              "sensitivity but more noise); larger values = coarser "
              "grouping (more tolerant to offsets/indels). Typical: "
              "150-300 for noisy long reads."),
    )
    ap.add_argument("--max-gap", type=int, default=1000,
        help=("maximum allowed gap in bp between successive k-mer seed points "
              "when chaining points into segments. The chaining step compares "
              "delta_i and delta_j to this threshold (points with both deltas "
              "<= --max-gap are chained). Smaller values make chaining stricter "
              "(require denser seeds); larger values allow sparser or offset matches "
              "(tolerant to indels/noise). Typical for noisy long reads: 2000-3000."),
    )
    ap.add_argument("--dotplot", action="store_true")
    ap.add_argument("--jobs", type=int, default=1,
                    help="number of worker processes to use for alignment stats (default: 1 = serial)")
    ap.add_argument("--workers", type=int, default=1,
                    help="number of parallel worker processes to use for per-read processing (default: 1)")

    ap.add_argument("--rc-repeats", action="store_true",
                    help="detect dominant RC junction and write 1-based cut coords")

    # keep only the trimmed-repeat related flags the user requested
    ap.add_argument("--direct-repeats", action="store_true",
                    help="write tight trimmed repeat coordinates to <out_prefix>.repeats.trimmed.filtered.tsv")
    ap.add_argument("--min-pct-id", type=float, default=0.0,
                    help="minimum percent identity (pct_id) required to write a trimmed repeat line")

    args = ap.parse_args()
    out_prefix_path = Path(args.out_prefix)
    if args.out_dir:
        Path(args.out_dir).mkdir(parents=True, exist_ok=True)
        out_prefix_path = Path(args.out_dir) / out_prefix_path.name
    out_prefix = str(out_prefix_path)

    # Process all records in the input (supports FASTA/FASTQ and .gz)
    any_written = False
    params = {
        'k': args.k,
        'min_length': args.min_length,
        'bin_size': args.bin_size,
        'max_gap': args.max_gap,
        'dotplot': args.dotplot,
        'jobs': args.jobs,
        'rc_repeats': args.rc_repeats,
        'min_pct_id': args.min_pct_id,
        'direct_repeats': args.direct_repeats,
        'write_trimmed_repeats': getattr(args, 'write_trimmed_repeats', False),
        'workers': args.workers,
    }

    # Use a ProcessPoolExecutor to process reads in parallel when requested.
    if args.workers and int(args.workers) > 1:
        max_workers = int(args.workers)
        # Create executor with 'fork' context to avoid pickling __main__-scoped
        # callables when the default start method is 'spawn' in some environments.
        mp_ctx = mp.get_context('fork')
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers, mp_context=mp_ctx) as exc:
            futures = set()
            # Submit tasks as we stream reads to avoid loading entire file into memory
            for rec_id, s in iter_reads_file(args.fasta):
                # bound outstanding futures to 2 * workers
                while len(futures) >= max_workers * 2:
                    done, futures = concurrent.futures.wait(futures, return_when=concurrent.futures.FIRST_COMPLETED)
                    for fut in done:
                        try:
                            rid, wrote, msgs = fut.result()
                            for m in msgs:
                                print(m)
                            any_written = any_written or wrote
                        except Exception as e:
                            print(f"Worker failed: {e}")
                futures.add(exc.submit(process_single_read, rec_id, s, str(out_prefix_path), args.fasta, params))

            # drain remaining futures
            for fut in concurrent.futures.as_completed(futures):
                try:
                    rid, wrote, msgs = fut.result()
                    for m in msgs:
                        print(m)
                    any_written = any_written or wrote
                except Exception as e:
                    print(f"Worker failed: {e}")
    else:
        # Sequential fallback
        for rec_id, s in iter_reads_file(args.fasta):
            rid, wrote, msgs = process_single_read(rec_id, s, str(out_prefix_path), args.fasta, params)
            for m in msgs:
                print(m)
            any_written = any_written or wrote

    if not any_written:
        print("No reads with repeats found; no per-read outputs written.")

if __name__ == "__main__":
    main()
