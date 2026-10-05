#!/usr/bin/env python3
"""Clustered heatmaps of expression: all DEGs, and the N most variable genes.

Input: a table whose first column is the gene name and whose other columns are
normalised expression values per sample. Values are log2(x + 1) transformed, then each
gene is scaled to 0-1 across samples (row-wise), so the colours show relative
expression per gene, not absolute levels. Samples are not clustered (column order kept).

Example
-------
python scripts/04_heatmap.py --input data/18mo_counts_heatmap.xlsx --label 18m --top-n 50
"""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from common import read_table

CMAP = "RdBu_r"  # red = high, blue = low


def heatmap(df, title, out, show_names):
    n = len(df)
    g = sns.clustermap(df, cmap=CMAP, standard_scale=0, col_cluster=False, row_cluster=True,
                       figsize=(10, max(8, 0.28 * n) if show_names else 12),
                       yticklabels=show_names, cbar_pos=(0.02, 0.8, 0.03, 0.15))
    g.ax_heatmap.set_title(title, pad=20)
    if show_names:
        g.ax_heatmap.tick_params(axis="y", labelsize=7)
    g.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(g.fig)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True)
    p.add_argument("--sheet", default=None)
    p.add_argument("--label", default="", help="suffix for output file names, e.g. 18m")
    p.add_argument("--top-n", type=int, default=50)
    p.add_argument("--out-dir", default="results")
    a = p.parse_args()

    df = read_table(a.input, a.sheet)
    df = df.set_index(df.columns[0])
    df = df.select_dtypes("number")
    log_df = np.log2(df + 1)

    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    tag = f"_{a.label}" if a.label else ""

    heatmap(log_df, f"All DEGs (n={len(log_df)})", out / f"heatmap_all{tag}.png", show_names=len(log_df) <= 80)

    top = log_df.var(axis=1).nlargest(a.top_n).index
    heatmap(log_df.loc[top], f"{len(top)} most variable genes", out / f"heatmap_top{a.top_n}{tag}.png",
            show_names=True)


if __name__ == "__main__":
    main()
