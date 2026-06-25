#!/usr/bin/env python3
"""
Build the cell-type proportions matrix that DESP needs, derived from THIS pipeline's scRNA
reference-mapping output (the Module-C proteogenomic integration hook).

Consumes the per-sample `<sample>_celltypes.csv` files emitted by REFERENCE_MAPPING (columns:
barcode index, sample_id, ref_cell_type, ...), computes each scRNA sample's cell-type fractions,
and assembles a `cell_type` x proteomics-sample matrix (the format DESP::DESP expects: proportions
transposed to samples x cell_types downstream).

scRNA sample ids rarely equal the bulk-proteomics sample ids, so an optional --map TSV
(rna_sample<TAB>prot_sample) renames/selects columns. Without a map, scRNA sample_id is used as-is.

Usage:
  prot_ms_proportions.py --celltypes a_celltypes.csv b_celltypes.csv [--map map.tsv] \
      [--celltype_col ref_cell_type] --out f157_celltype_proportions.tsv
"""
from __future__ import annotations

import argparse

import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--celltypes", nargs="+", required=True,
                    help="one or more <sample>_celltypes.csv files from REFERENCE_MAPPING")
    ap.add_argument("--celltype_col", default="ref_cell_type")
    ap.add_argument("--sample_col", default="sample_id")
    ap.add_argument("--map", default=None, help="TSV: rna_sample<TAB>prot_sample")
    ap.add_argument("--out", default="celltype_proportions.tsv")
    a = ap.parse_args()

    frames = []
    for path in a.celltypes:
        df = pd.read_csv(path)
        if a.celltype_col not in df.columns or a.sample_col not in df.columns:
            raise SystemExit(f"{path} missing '{a.celltype_col}' or '{a.sample_col}'")
        frames.append(df[[a.sample_col, a.celltype_col]])
    cells = pd.concat(frames, ignore_index=True).dropna()

    # per-sample cell-type fractions -> cell_type x sample matrix
    counts = cells.groupby([a.sample_col, a.celltype_col]).size().unstack(fill_value=0)
    props = counts.div(counts.sum(axis=1), axis=0)        # rows=sample, cols=cell_type
    mat = props.T                                          # rows=cell_type, cols=sample

    if a.map:
        m = pd.read_csv(a.map, sep="\t", header=None, names=["rna_sample", "prot_sample"], dtype=str)
        rename = dict(zip(m["rna_sample"], m["prot_sample"]))
        keep = [c for c in mat.columns if c in rename]
        if not keep:
            raise SystemExit("No scRNA samples in --celltypes matched --map rna_sample column.")
        mat = mat[keep].rename(columns=rename)

    mat = mat.reset_index().rename(columns={mat.index.name or "index": "cell_type"})
    if "cell_type" not in mat.columns:
        mat = mat.rename(columns={mat.columns[0]: "cell_type"})
    mat.to_csv(a.out, sep="\t", index=False)
    print(f"proportions: {mat.shape[0]} cell types x {mat.shape[1]-1} samples -> {a.out}")


if __name__ == "__main__":
    main()
