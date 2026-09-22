"""Alignment-to-table conversion, filtering, and per-column nucleotide frequency."""
from __future__ import annotations

import numpy as np
import pandas as pd

NUCLEOTIDES = ["A", "C", "T", "G", "N"]
GAP = "-"


def alignment_to_table(aligned: dict[str, str]) -> pd.DataFrame:
    """Turn an ``{name: seq}`` alignment into a sequences x positions table.

    Rows are sequence names, columns are ``N1..Nk`` (1-based alignment
    position), values are the single-character base at that position.
    """
    names = list(aligned.keys())
    length = len(next(iter(aligned.values())))
    data = np.array([list(seq) for seq in aligned.values()])
    columns = [f"N{i}" for i in range(1, length + 1)]
    return pd.DataFrame(data, index=names, columns=columns)


def nuc_table_filter(table: pd.DataFrame, max_miss_per: float = 0.2, remove_mono: bool = True) -> pd.DataFrame:
    """Drop alignment columns that are too gappy and/or monomorphic."""
    keep = []
    n = len(table)
    for col in table.columns:
        vals = table[col]
        miss = (vals == GAP).sum() / n if n else 0
        if miss > max_miss_per:
            continue
        if remove_mono:
            non_gap = vals[vals != GAP]
            if non_gap.nunique() <= 1:
                continue
        keep.append(col)
    return table[keep].copy()


def nuc_frequency(table: pd.DataFrame) -> pd.DataFrame:
    """Per-column proportion of A/C/T/G/N (and implicitly gap) bases."""
    n = len(table)
    rows = {}
    for base in NUCLEOTIDES:
        rows[base] = (table == base).sum(axis=0) / n if n else 0
    return pd.DataFrame(rows).T
