#!/usr/bin/env python3
"""
Flatten a CloneTracer *_out.pickle into a tidy per-cell clone-assignment CSV.

Picks the selected tree with the lowest final ELBO, then for each cell records the
argmax clone and its posterior probability. Output columns:
    barcode, sample, clone, max_prob, tree
Consumable by the DDE_32 downstream integration + the per-patient clone-overlay figures.
"""
import csv
import pickle
import sys

import numpy as np


def main():
    inp, out = sys.argv[1], sys.argv[2]
    with open(inp, "rb") as f:
        obj = pickle.load(f)

    cp = obj.get("clonal_prob", {})
    bc = obj.get("cell_barcode", [])

    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["barcode", "sample", "clone", "max_prob", "tree"])
        if not cp:
            # Trivial / single-clone result (e.g. healthy control): header only.
            return

        elbo = obj.get("ELBO", {})

        def final_elbo(t):
            e = elbo.get(t)
            return float(np.asarray(e).ravel()[-1]) if e is not None else np.inf

        best = min(cp.keys(), key=final_elbo) if elbo else next(iter(cp))
        P = np.asarray(cp[best])              # cells x clones
        clone = P.argmax(axis=1)
        mprob = P.max(axis=1)
        for i, b in enumerate(bc):
            sample = b.split("__")[0] if "__" in str(b) else ""
            w.writerow([b, sample, int(clone[i]), round(float(mprob[i]), 4), int(best)])


if __name__ == "__main__":
    main()
