"""Neighbour-joining phylogenetic tree construction (mirrors R's ``getTree``,
built on ``ape``'s NJ algorithm there; here on Biopython's)."""
from __future__ import annotations

import pandas as pd
from Bio.Phylo.BaseTree import Tree
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor


def build_tree(dist_df: pd.DataFrame) -> Tree:
    names = list(dist_df.index)
    matrix = [list(dist_df.iloc[i, : i + 1]) for i in range(len(names))]
    dm = DistanceMatrix(names=names, matrix=matrix)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)
    tree.root_at_midpoint()
    return tree


def tree_tip_order(tree: Tree) -> list[str]:
    return [t.name for t in tree.get_terminals()]
