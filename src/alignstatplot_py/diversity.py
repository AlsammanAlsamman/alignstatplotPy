"""Population-genetics summary statistics computed directly from an
aligned sequence table (mirrors R's ``nucDiversityStats``,
``positionConservation``, and ``regionStats``)."""
from __future__ import annotations

import math
from itertools import combinations

import numpy as np
import pandas as pd

GAP = "-"
TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
VALID_BASES = {"A", "C", "T", "G"}


def diversity_stats(table: pd.DataFrame) -> pd.DataFrame:
    """Segregating sites, nucleotide diversity (pi), Watterson's theta, Ts/Tv."""
    n = len(table)
    length = table.shape[1]

    segregating = 0
    ts, tv = 0, 0
    pairwise_diffs = []
    seqs = table.values

    pair_count = n * (n - 1) // 2 if n > 1 else 0
    diff_sum = 0
    compared_sum = 0

    for col in range(length):
        column = seqs[:, col]
        valid = [b for b in column if b in VALID_BASES]
        if len(set(valid)) > 1:
            segregating += 1
        for i, j in combinations(range(n), 2):
            a, b = column[i], column[j]
            if a in VALID_BASES and b in VALID_BASES:
                compared_sum += 1
                if a != b:
                    diff_sum += 1
                    if (a, b) in TRANSITIONS:
                        ts += 1
                    else:
                        tv += 1

    nucleotide_diversity = diff_sum / compared_sum if compared_sum else 0.0
    harmonic = sum(1 / i for i in range(1, n)) if n > 1 else 1
    watterson_theta = segregating / harmonic if harmonic else 0.0
    ts_tv_ratio = ts / tv if tv else float("nan")

    return pd.DataFrame(
        [
            {
                "NumSequences": n,
                "AlignmentLength": length,
                "SegregatingSites": segregating,
                "NucleotideDiversity": nucleotide_diversity,
                "WattersonsTheta": watterson_theta,
                "TsTvRatio": ts_tv_ratio,
            }
        ]
    )


def position_conservation(table: pd.DataFrame) -> pd.DataFrame:
    """Majority allele, conservation fraction, Shannon entropy, and
    missingness at every alignment position."""
    n = len(table)
    rows = []
    for pos, col in enumerate(table.columns, start=1):
        column = table[col]
        valid = column[column != GAP]
        missing = (len(column) - len(valid)) / len(column) if len(column) else 0.0
        if len(valid) == 0:
            rows.append(
                {"Position": pos, "MajorityAllele": GAP, "Conservation": 0.0, "Entropy": 0.0, "Missing": missing}
            )
            continue
        counts = valid.value_counts()
        majority = counts.index[0]
        conservation = counts.iloc[0] / len(valid)
        freqs = counts / len(valid)
        entropy = -sum(p * math.log2(p) for p in freqs if p > 0)
        rows.append(
            {
                "Position": pos,
                "MajorityAllele": majority,
                "Conservation": conservation,
                "Entropy": entropy,
                "Missing": missing,
            }
        )
    return pd.DataFrame(rows)


def region_stats(position_stats: pd.DataFrame, anno: pd.DataFrame) -> pd.DataFrame:
    """Aggregate :func:`position_conservation` output within each annotated
    region (e.g. exon/intron)."""
    rows = []
    for _, region in anno.iterrows():
        mask = (position_stats["Position"] >= region["Start"]) & (position_stats["Position"] <= region["End"])
        subset = position_stats[mask]
        if subset.empty:
            continue
        rows.append(
            {
                "GeneName": region.get("GeneName", ""),
                "Region": region.get("Region", ""),
                "Start": region["Start"],
                "End": region["End"],
                "NumPositions": len(subset),
                "MeanConservation": subset["Conservation"].mean(),
                "MeanEntropy": subset["Entropy"].mean(),
                "MeanMissing": subset["Missing"].mean(),
                "VariantRate": (subset["Conservation"] < 1.0).mean(),
            }
        )
    return pd.DataFrame(rows)
