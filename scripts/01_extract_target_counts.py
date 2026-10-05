#!/usr/bin/env python3
"""Pull expression rows for a gene list out of a big table, optionally adding gene names.

Used for every "give me the counts of my DEGs" step of the project (gene symbols for
P3/P10, XLOC ids for P20/18 months). Matching is case-insensitive and ignores spaces.

Examples
--------
# match by gene symbol
python scripts/01_extract_target_counts.py \
    --targets data/genelists_all_ages.xlsx --targets-sheet p10 --targets-col Gene \
    --table data/CORTEX_strain_and_age_2021.xls --table-sheet "P10_Padj<0.05" --table-col "New Gene" \
    --out results/p10_targeted_counts.xlsx

# match by XLOC id and attach gene names from an annotation file
python scripts/01_extract_target_counts.py \
    --targets data/p20_XLOC.xlsx --targets-col XLOC \
    --table data/COUNTS_PFC_EXONS.xlsx --table-col "Unnamed: 0" \
    --annotation data/annotation_xloc.xlsx --annotation-key XLOC \
    --out results/p20_targeted_counts_fin.xlsx
"""
import argparse
from pathlib import Path

import pandas as pd

from common import norm_symbols, read_table


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--targets", required=True, help="file with the gene list")
    p.add_argument("--targets-sheet", default=None)
    p.add_argument("--targets-col", required=True, help="column with gene ids in the gene list")
    p.add_argument("--table", required=True, help="big expression table to filter")
    p.add_argument("--table-sheet", default=None)
    p.add_argument("--table-col", required=True, help="column with gene ids in the big table")
    p.add_argument("--annotation", default=None, help="optional table mapping ids to gene names")
    p.add_argument("--annotation-key", default=None, help="id column in the annotation table")
    p.add_argument("--out", required=True)
    a = p.parse_args()

    targets = read_table(a.targets, a.targets_sheet)
    table = read_table(a.table, a.table_sheet)

    wanted = set(norm_symbols(targets[a.targets_col]))
    keys = table[a.table_col].astype(str).str.strip().str.upper()
    matches = table[keys.isin(wanted)].copy()

    missing = wanted - set(keys)
    print(f"{len(wanted)} target ids, {len(matches)} rows matched, {len(missing)} not found")
    if missing:
        print("  not found (first 10):", sorted(missing)[:10])

    if a.annotation:
        if not a.annotation_key:
            p.error("--annotation needs --annotation-key")
        annot = read_table(a.annotation)
        matches = matches.rename(columns={a.table_col: a.annotation_key})
        for df in (matches, annot):
            df[a.annotation_key] = df[a.annotation_key].astype(str).str.strip()
        matches = matches.merge(annot, on=a.annotation_key, how="left")

    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    matches.to_excel(a.out, index=False)
    print("saved", a.out)


if __name__ == "__main__":
    main()
