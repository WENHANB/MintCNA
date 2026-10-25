"""Bin-level ground-truth log2 ratios of one simulation setting, for Supplementary Fig. 1.

Total copy number per bin and cell from CNAsim's profiles.tsv (the two haplotype values summed),
converted to log2(CN / 2) with CN = 0 floored at 0.25. Rows follow the bin order of the processed
log2-ratio matrix so the three matrices of the figure line up.

    python reproduce/analysis/simulation/denoising_comparison/build_truth.py
Input : <setting>/profiles.tsv, <setting>/l2r_smoothed_all_cells_100kb.csv
Output: <setting>/FigureS1/data/truth_log2.csv
Run with MINTCNA_DATA set (see reproduce/README.md).
"""
import os
import numpy as np
import pandas as pd

ROOT = os.path.join(os.environ["MINTCNA_DATA"], "Simulation/cnasim/c3_N100_100kb_l1e4_n0.1")
PROFILES = os.path.join(ROOT, "profiles.tsv")
L2R = os.path.join(ROOT, "l2r_smoothed_all_cells_100kb.csv")
OUT = os.path.join(ROOT, "FigureS1", "data", "truth_log2.csv")
CN_FLOOR = 0.25


def main():
    l2r = pd.read_csv(L2R, index_col=0)
    bins, cells = list(l2r.index), list(l2r.columns)
    bin_index = {(lab.split(":")[0], int(lab.split(":")[1].split("-")[0])): i for i, lab in enumerate(bins)}
    cell_idx = {c: i for i, c in enumerate(cells)}

    cn = np.full((len(cells), len(bins)), np.nan, dtype=np.float32)
    with open(PROFILES) as f:
        f.readline()
        for line in f:
            cell, chrom, start, _end, states = line.rstrip("\n").split("\t")
            if cell not in cell_idx:
                continue
            j = bin_index.get((chrom, int(start)))
            if j is None:
                continue
            a, b = states.split(",")
            cn[cell_idx[cell], j] = int(a) + int(b)
    if np.isnan(cn).any():
        print(f"{int(np.isnan(cn).sum())} missing entries filled with 2 (diploid)")
        cn[np.isnan(cn)] = 2.0

    log2r = np.log2(np.maximum(cn, CN_FLOOR) / 2.0)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    pd.DataFrame(log2r.T, index=bins, columns=cells).to_csv(OUT)
    print(f"saved {OUT}  ({len(bins)} bins x {len(cells)} cells)")


if __name__ == "__main__":
    main()
