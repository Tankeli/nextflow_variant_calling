#!/usr/bin/env python3
"""Build a proportion-controlled barcode list for a souporcell titration run.

Given the full combined (``<sample>__<bc>``) barcode list of a 2-sample mix, draw a fixed TOTAL of
cells split so the minority sample contributes ``--minority-frac`` of them. Total is held constant
across ratios so any accuracy change is due to proportion, not cluster size. Sampling is without
replacement and seeded.

Usage:
  souporcell_subsample_barcodes.py --bclist mix.barcodes.tsv --minority-sample PBMMC_2 \
      --total 4000 --minority-frac 0.05 --seed 0 --out r05/barcodes.txt
"""
import argparse
import os
import random


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bclist", required=True)
    ap.add_argument("--minority-sample", required=True)
    ap.add_argument("--total", type=int, required=True)
    ap.add_argument("--minority-frac", type=float, required=True)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    bcs = [l.strip() for l in open(a.bclist) if l.strip()]
    minp = a.minority_sample + "__"
    minority = [b for b in bcs if b.startswith(minp)]
    majority = [b for b in bcs if not b.startswith(minp)]
    if not minority or not majority:
        raise SystemExit(f"minority/majority empty (min={len(minority)} maj={len(majority)}) "
                         f"— check --minority-sample matches a prefix in {a.bclist}")

    n_min = round(a.total * a.minority_frac)
    n_maj = a.total - n_min
    n_min = min(n_min, len(minority))
    n_maj = min(n_maj, len(majority))

    rng = random.Random(a.seed)
    pick = rng.sample(minority, n_min) + rng.sample(majority, n_maj)
    rng.shuffle(pick)

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w") as f:
        f.write("\n".join(pick) + "\n")
    print(f"{a.out}: {n_min} minority ({a.minority_sample}) + {n_maj} majority = {len(pick)} cells "
          f"(target frac {a.minority_frac}, actual {n_min/len(pick):.4f})")


if __name__ == "__main__":
    main()
