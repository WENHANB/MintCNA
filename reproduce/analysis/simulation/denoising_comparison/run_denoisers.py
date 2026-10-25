"""Denoise the log2 ratios of one simulation setting with AG-DCAE and with the CAE + MSE baseline,
for Supplementary Fig. 1 (setting c3_N100_100kb_l1e4_n0.1, latent 200, seed 0).

    python reproduce/analysis/simulation/denoising_comparison/run_denoisers.py [agdcae|caemse|both]
Input : <setting>/l2r_smoothed_all_cells_100kb.csv
Output: <setting>/FigureS1/data/agdcae_denoised.csv, caemse_denoised.csv
Run from the root of the MintCNA repository with MINTCNA_DATA set (see reproduce/README.md).
"""
import os
import sys
import numpy as np
import pandas as pd
import torch

sys.path.insert(0, "MintCNA")                                                   # ag_dcae package
sys.path.insert(0, "reproduce/analysis/simulation/denoising_comparison")        # cae_mse
from ag_dcae.train import train_ag_dcae, xs_gen
from cae_mse import train_cae_mse

ROOT = os.path.join(os.environ["MINTCNA_DATA"], "Simulation/cnasim/c3_N100_100kb_l1e4_n0.1")
INPUT = os.path.join(ROOT, "l2r_smoothed_all_cells_100kb.csv")
OUT_DIR = os.path.join(ROOT, "FigureS1", "data")

COMMON = dict(latent_dim=200, kernel_size=5, lr=1e-4, batch_size=16, sigma=0.08,
              min_delta=0.001, test_size=0.2, random_state=30)
HP_AGDCAE = dict(COMMON, epochs=200, patience=30)   # more epochs: the 80-epoch run had not converged on this setting
HP_CAEMSE = dict(COMMON, epochs=80, patience=15)


def infer(model, data, batch_size, device):
    model.eval()
    out = []
    for _, x in xs_gen(data, batch_size, shuffle=False):
        x = torch.from_numpy(x).float().to(device).unsqueeze(1)
        with torch.no_grad():
            _, y = model(x)
        out.append(y.cpu().numpy())
    return np.concatenate(out, axis=0)


def main(which):
    df = pd.read_csv(INPUT, index_col=0)
    data = df.values.T                                   # cells x bins
    device = "cuda" if torch.cuda.is_available() else "cpu"
    os.makedirs(OUT_DIR, exist_ok=True)
    print(f"{data.shape[0]} cells x {data.shape[1]} bins, device {device}", flush=True)

    if which in ("agdcae", "both"):
        torch.manual_seed(0); np.random.seed(0)
        model, _, _ = train_ag_dcae(data, alpha=0.1, edge_lambda=0.2, device=device, **HP_AGDCAE)
        out = infer(model, data, HP_AGDCAE["batch_size"], torch.device(device))
        pd.DataFrame(out.T, index=df.index, columns=df.columns).to_csv(os.path.join(OUT_DIR, "agdcae_denoised.csv"))
        print("saved agdcae_denoised.csv", flush=True)

    if which in ("caemse", "both"):
        torch.manual_seed(0); np.random.seed(0)
        model, _, _ = train_cae_mse(data, device=device, **HP_CAEMSE)
        out = infer(model, data, HP_CAEMSE["batch_size"], torch.device(device))
        pd.DataFrame(out.T, index=df.index, columns=df.columns).to_csv(os.path.join(OUT_DIR, "caemse_denoised.csv"))
        print("saved caemse_denoised.csv", flush=True)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "both")
