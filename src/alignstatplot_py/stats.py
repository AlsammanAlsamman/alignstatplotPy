"""Per-sequence alignment statistics (mirrors R's ``AlignmentStatsPerSeq``)."""
from __future__ import annotations

import pandas as pd


def per_sequence_stats(aligned: dict[str, str]) -> pd.DataFrame:
    """Per-sequence base counts, gap %, and GC % of an alignment."""
    rows = []
    for i, (name, seq) in enumerate(aligned.items(), start=1):
        seq = seq.upper()
        length = len(seq)
        a, t, c, g = seq.count("A"), seq.count("T"), seq.count("C"), seq.count("G")
        gap = seq.count("-")
        gc = c + g
        rows.append(
            {
                "SeqNo": i,
                "Sequence.Name": name,
                "Sequence.Length": length,
                "A": a,
                "T": t,
                "C": c,
                "G": g,
                "Gap": gap,
                "Gap.Percentage": round(100 * gap / length, 2) if length else 0.0,
                "GC": gc,
                "GC.Percentage": round(100 * gc / length, 2) if length else 0.0,
            }
        )
    return pd.DataFrame(rows)
