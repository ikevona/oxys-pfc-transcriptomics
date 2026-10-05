#!/usr/bin/env python3
"""Overlap of DEG lists across ages: 4-set Venn diagram + a gene x age membership table.

Each input is a table with a `Gene` column. The membership table answers
"which genes are changed at several ages / reappear at 18 months?".

Example
-------
python scripts/03_venn_overlap.py --p3 data/p3_rgd_degs.xlsx --p10 data/p10_rgd_degs.xlsx \
    --p20 data/p20_rgd_degs.xlsx --m18 data/18mo_degs_rgd.xlsx --out-dir results
"""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from venn import venn

from common import norm_symbols, read_table

AGES = ["P3", "P10", "P20", "18mo"]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for flag in ("p3", "p10", "p20", "m18"):
        p.add_argument(f"--{flag}", required=True)
    p.add_argument("--gene-col", default="Gene")
    p.add_argument("--out-dir", default="results")
    a = p.parse_args()

    files = dict(zip(AGES, [a.p3, a.p10, a.p20, a.m18]))
    sets = {age: set(norm_symbols(read_table(f)[a.gene_col])) for age, f in files.items()}

    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(10, 10))
    venn(sets, cmap=["#D6EAF8", "#AED6F1", "#5DADE2", "#1F4E79"])
    plt.title("Overlap of DEGs across ages", fontsize=18, fontweight="bold", color="#16324F")
    plt.savefig(out / "venn.png", dpi=300, bbox_inches="tight", transparent=True)
    plt.close()

    genes = sorted(set().union(*sets.values()))
    member = pd.DataFrame({age: [g in s for g in genes] for age, s in sets.items()}, index=genes)
    member["n_ages"] = member[AGES].sum(axis=1)
    member.index.name = "gene"
    member.sort_values("n_ages", ascending=False).to_csv(out / "deg_membership.csv")

    print(member.n_ages.value_counts().sort_index(ascending=False).to_string())
    early = member[["P3", "P10", "P20"]].any(axis=1)
    print("DEGs at 18mo that were also DEGs at an early age:", int((early & member["18mo"]).sum()))


if __name__ == "__main__":
    main()
