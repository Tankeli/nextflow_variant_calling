#!/usr/bin/env python3
"""
Protein/ADT normalization — ported from scripts/13_prot_normalization.ipynb.

Normalizes protein counts with CLR (centered log-ratio) or DSB. Raw counts are stashed in the
`counts` layer first. DSB needs the raw/empty-droplet matrix (table of droplets); pass it with
--raw_h5, otherwise the script falls back to CLR (the notebook does the same on DSB failure).
Checkpoint: prot_02_normalized.h5mu.

Usage: prot_normalize.py --in IN.h5mu --sample S --method clr|dsb [--raw_h5 R] --out OUT.h5mu
"""
import argparse
import os
import warnings

import anndata as ad
import muon as mu
import scanpy as sc

warnings.filterwarnings("ignore")
os.environ["PYTHONWARNINGS"] = "ignore"
sc.settings.verbosity = 0
ad.settings.allow_write_nullable_strings = True


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--in", dest="inp", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--method", default="clr", choices=["clr", "dsb"])
    p.add_argument("--raw_h5", default=None, help="raw_feature_bc_matrix.h5 for DSB background")
    p.add_argument("--out", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    mdata = mu.read(args.inp)
    mdata["prot"].layers["counts"] = mdata["prot"].X.copy()

    method = args.method
    if method == "dsb":
        if not args.raw_h5:
            print("DSB requested but no --raw_h5; falling back to CLR")
            method = "clr"
        else:
            try:
                mdata_raw = mu.read_10x_h5(args.raw_h5)
                mu.prot.pp.dsb(mdata, mdata_raw)
                print("DSB normalization complete")
            except Exception as e:
                print(f"DSB failed ({e}); falling back to CLR")
                method = "clr"

    if method == "clr":
        mu.prot.pp.clr(mdata)
        print("CLR normalization complete")

    mdata.write(args.out)
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
