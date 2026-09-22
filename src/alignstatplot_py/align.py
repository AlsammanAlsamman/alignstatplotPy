"""Lightweight multiple sequence alignment.

The R package shells out to ClustalW / Clustal Omega / MUSCLE. To keep the
Python port dependency-light and fully self-contained (no external binary
required), this module implements a simple progressive aligner: every
sequence is globally pairwise-aligned to a running consensus profile using
Biopython's ``PairwiseAligner``. This is *not* a drop-in replacement for a
production MSA tool on large/divergent datasets, but it produces a valid
gapped alignment of consistent width and is sufficient for the statistics
and visualisations this package builds on.

If you already have an alignment (e.g. from MAFFT/ClustalW/MUSCLE), skip
this step entirely and load it directly with :func:`io_utils.read_fasta`
or Biopython's ``AlignIO`` -- every downstream function only needs an
``{name: aligned_sequence}`` mapping of equal-length strings.
"""
from __future__ import annotations

from collections import OrderedDict

from Bio.Align import PairwiseAligner

GAP = "-"


def _aligner() -> PairwiseAligner:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5
    return aligner


def is_aligned(seqs: dict[str, str]) -> bool:
    lengths = {len(s) for s in seqs.values()}
    return len(lengths) == 1


def align_sequences(seqs: dict[str, str], method: str = "auto") -> "OrderedDict[str, str]":
    """Align sequences, returning an equal-length gapped alignment.

    ``method`` is accepted for API/CLI parity with the R package
    (``"ClustalW"``, ``"ClustalOmega"``, ``"Muscle"``, ``"auto"``) but this
    implementation always uses the same built-in progressive aligner; the
    parameter is kept so scripts/config files written for either tool stay
    valid.
    """
    if is_aligned(seqs):
        return OrderedDict(seqs)

    names = list(seqs.keys())
    ordered = sorted(names, key=lambda n: len(seqs[n]), reverse=True)
    aligner = _aligner()

    profile: "OrderedDict[str, str]" = OrderedDict()
    ref_name = ordered[0]
    profile[ref_name] = seqs[ref_name]

    for name in ordered[1:]:
        query = seqs[name]
        reference = profile[ref_name]
        alignment = aligner.align(reference.replace(GAP, ""), query)[0]
        aligned_ref, aligned_query = str(alignment[0]), str(alignment[1])

        # Merge any new gap columns the pairwise alignment introduced in the
        # reference into every sequence aligned so far.
        if GAP in aligned_ref:
            insert_positions = [i for i, c in enumerate(aligned_ref) if c == GAP]
            for existing_name, existing_seq in profile.items():
                profile[existing_name] = _insert_gaps(existing_seq, insert_positions)
        profile[name] = aligned_query

    max_len = max(len(s) for s in profile.values())
    for name, seq in profile.items():
        if len(seq) < max_len:
            profile[name] = seq + GAP * (max_len - len(seq))

    return OrderedDict((name, profile[name]) for name in names)


def _insert_gaps(seq: str, positions: list[int]) -> str:
    # positions are indices into the *reference-with-gaps* coordinate space;
    # since we always align against the current (already-gapped) reference,
    # inserting at those same indices keeps every sequence in lock-step.
    out = list(seq)
    for pos in sorted(positions):
        pos = min(pos, len(out))
        out.insert(pos, GAP)
    return "".join(out)
