#!/usr/bin/env python3
"""Volcano plots and up/down DEG counts for P3, P10, P20 and 18 months.

Input: one Excel workbook with one sheet per age and columns `gene`, `log2FC`, `padj`.
DEG definition used everywhere: padj < 0.05 (no fold-change cut-off).

Direction: "up" means log2FC > 0 *as written in the table*. If your contrast was
Wistar vs OXYS rather than OXYS vs Wistar, pass --flip-sign so that "up" means
higher in OXYS.

Example
-------
python scripts/02_deg_summary.py --input data/volcano_data.xlsx --out-dir results \
    [--gene Hebp1] [--flip-sign]
"""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# label used in plots -> sheet name in the workbook
AGES = {"P3": "p3", "P10": "p10", "P20": "p20", "18mo": "m18"}
PADJ = 0.05
COLORS = {"UP": "#d62728", "DOWN": "#1f77b4", "NS": "lightgray"}


def load(path, sheet, flip):
    df = pd.read_excel(path, sheet_name=sheet)
    df["log2FC"] = pd.to_numeric(df["log2FC"], errors="coerce") * (-1 if flip else 1)
    df["padj"] = pd.to_numeric(df["padj"], errors="coerce")
    df = df.dropna(subset=["log2FC", "padj"]).copy()
    df["padj"] = df["padj"].clip(lower=1e-300)  # avoid log10(0)
    df["status"] = "NS"
    df.loc[(df.padj < PADJ) & (df.log2FC > 0), "status"] = "UP"
    df.loc[(df.padj < PADJ) & (df.log2FC < 0), "status"] = "DOWN"
    return df


def volcano_grid(data, out):
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    for ax, (age, df) in zip(axes.flat, data.items()):
        for cat in ("NS", "DOWN", "UP"):
            sub = df[df.status == cat]
            ax.scatter(sub.log2FC, -np.log10(sub.padj), s=10, alpha=0.7, c=COLORS[cat],
                       label=cat if age == "P3" else None)
        ax.axhline(-np.log10(PADJ), ls="--", lw=1, c="gray")
        ax.set(title=age, xlabel="log2FC", ylabel="-log10(padj)")
    axes.flat[0].legend(frameon=False)
    fig.suptitle("Differential expression, OXYS vs Wistar", fontsize=14)
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)


def count_plot(counts, out):
    x = np.arange(len(counts))
    w = 0.35
    fig, ax = plt.subplots(figsize=(8, 5))
    up = ax.bar(x - w / 2, counts["up"], w, label="Up-regulated", color=COLORS["UP"], edgecolor="black")
    down = ax.bar(x + w / 2, counts["down"], w, label="Down-regulated", color=COLORS["DOWN"], edgecolor="black")
    ax.bar_label(up, padding=3)
    ax.bar_label(down, padding=3)
    ax.set_xticks(x, counts.index)
    ax.set(ylabel="Number of genes", title=f"Significant DEGs (padj < {PADJ})")
    ax.grid(axis="y", ls="--", alpha=0.3)
    ax.set_axisbelow(True)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)


def gene_dynamics(data, gene):
    print(f"\n{gene}")
    for age, df in data.items():
        row = df[df["gene"].astype(str).str.upper() == gene.upper()]
        if row.empty:
            print(f"  {age}: not found")
            continue
        r = row.iloc[0]
        print(f"  {age}: {r.status:<4} log2FC={r.log2FC:.3f} padj={r.padj:.3g}")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True)
    p.add_argument("--out-dir", default="results")
    p.add_argument("--flip-sign", action="store_true", help="multiply log2FC by -1")
    p.add_argument("--gene", help="also print status of one gene at every age")
    a = p.parse_args()

    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    data = {age: load(a.input, sheet, a.flip_sign) for age, sheet in AGES.items()}

    counts = pd.DataFrame({age: {"up": (d.status == "UP").sum(), "down": (d.status == "DOWN").sum()}
                           for age, d in data.items()}).T
    counts["total"] = counts.up + counts.down
    print(counts)
    counts.to_csv(out / "deg_counts.csv")

    volcano_grid(data, out / "volcano_grid.png")
    count_plot(counts, out / "deg_counts.png")
    if a.gene:
        gene_dynamics(data, a.gene)


if __name__ == "__main__":
    main()
