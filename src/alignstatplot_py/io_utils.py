"""FASTA / annotation / alignment file reading utilities."""
from __future__ import annotations

from collections import OrderedDict
from pathlib import Path

import pandas as pd


def read_fasta(path: str | Path) -> "OrderedDict[str, str]":
    """Read a FASTA file into an ordered {name: sequence} mapping.

    Sequence names are the FASTA header up to the first whitespace,
    matching how the R package's ``getSeqInfo`` names sequences.
    """
    seqs: "OrderedDict[str, str]" = OrderedDict()
    name = None
    chunks: list[str] = []
    with open(path, "r") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        seqs[name] = "".join(chunks).upper()
    if not seqs:
        raise ValueError(f"No sequences found in FASTA file: {path}")
    return seqs


def read_annotation(path: str | Path) -> pd.DataFrame:
    """Read a gene annotation table.

    Expected whitespace/tab separated columns (matching the R package's
    example annotation format)::

        GeneName  Region  Start  End  Direction
        gene1     intron  1      14   forward
    """
    df = pd.read_csv(path, sep=r"\s+", header=None, engine="python")
    ncols = df.shape[1]
    names = ["GeneName", "Region", "Start", "End", "Direction"][:ncols]
    df.columns = names
    if "Start" in df.columns:
        df["Start"] = df["Start"].astype(int)
    if "End" in df.columns:
        df["End"] = df["End"].astype(int)
    return df


def write_fasta(seqs: dict[str, str], path: str | Path) -> None:
    with open(path, "w") as fh:
        for name, seq in seqs.items():
            fh.write(f">{name}\n{seq}\n")
