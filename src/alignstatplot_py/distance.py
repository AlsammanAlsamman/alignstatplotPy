"""Pairwise sequence distance matrix (mirrors R's ``getDistanceMatrixTabel``)."""
from __future__ import annotations

import numpy as np
import pandas as pd

VALID_BASES = {"A", "C", "T", "G"}


def distance_matrix(table: pd.DataFrame) -> pd.DataFrame:
    """Pairwise p-distance (fraction of differing, non-gap positions)."""
    names = list(table.index)
    n = len(names)
    values = table.values
    dist = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            a, b = values[i], values[j]
            valid = np.isin(a, list(VALID_BASES)) & np.isin(b, list(VALID_BASES))
            compared = valid.sum()
            if compared == 0:
                d = 0.0
            else:
                diffs = (a[valid] != b[valid]).sum()
                d = diffs / compared
            dist[i, j] = dist[j, i] = d
    return pd.DataFrame(dist, index=names, columns=names)
