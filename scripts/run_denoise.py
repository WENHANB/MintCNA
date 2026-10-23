#!/usr/bin/env python
"""CLI entry point for AG-DCAE denoising.

Usage:
    python scripts/run_denoise.py --input normalized_dna.csv --output denoised_dna.csv
    python scripts/run_denoise.py --input normalized_dna.csv --output denoised_dna.csv --config config.yaml
"""

import argparse
import yaml
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "MintCNA"))
from ag_dcae.denoise import denoise


def main():
    parser = argparse.ArgumentParser(description="AG-DCAE Denoising for MintCNA")
    parser.add_argument("--input", required=True, help="Input CSV: normalized log2-ratio matrix (bins x cells)")
    parser.add_argument("--output", required=True, help="Output CSV: denoised matrix (bins x cells)")
    parser.add_argument("--latent", default=None, help="Output CSV: latent embeddings (cells x latent_dim)")
    parser.add_argument("--loss_plot", default=None, help="Output PDF: training loss curve")
    parser.add_argument("--config", default=None, help="YAML config file (overrides defaults)")

    # Allow individual parameter overrides
    parser.add_argument("--latent_dim", type=int, default=None)
    parser.add_argument("--kernel_size", type=int, default=None)
    parser.add_argument("--lr", type=float, default=None)
    parser.add_argument("--epochs", type=int, default=None)
    parser.add_argument("--batch_size", type=int, default=None)
    parser.add_argument("--alpha", type=float, default=None)
    parser.add_argument("--edge_lambda", type=float, default=None)
    parser.add_argument("--sigma", type=float, default=None)
    parser.add_argument("--patience", type=int, default=None)
    parser.add_argument("--device", type=str, default=None)

    args = parser.parse_args()

    # Start with defaults
    params = dict(
        latent_dim=250,
        kernel_size=5,
        lr=0.0001,
        epochs=100,
        batch_size=16,
        alpha=0.1,
        edge_lambda=0.2,
        sigma=0.08,
        patience=15,
        min_delta=0.001,
        test_size=0.2,
        random_state=30,
        device=None,
    )

    # Override from config file
    if args.config is not None:
        with open(args.config) as f:
            cfg = yaml.safe_load(f)
        if "denoise" in cfg:
            dc = cfg["denoise"]
            param_map = {
                "latent_dim": "latent_dim",
                "kernel_size": "kernel_size",
                "learning_rate": "lr",
                "epochs": "epochs",
                "batch_size": "batch_size",
                "alpha": "alpha",
                "edge_lambda": "edge_lambda",
                "sigma": "sigma",
                "patience": "patience",
                "min_delta": "min_delta",
                "test_size": "test_size",
                "random_state": "random_state",
                "device": "device",
            }
            for cfg_key, param_key in param_map.items():
                if cfg_key in dc and dc[cfg_key] is not None:
                    params[param_key] = dc[cfg_key]

    # Override from CLI flags (highest priority)
    cli_overrides = {
        "latent_dim": args.latent_dim,
        "kernel_size": args.kernel_size,
        "lr": args.lr,
        "epochs": args.epochs,
        "batch_size": args.batch_size,
        "alpha": args.alpha,
        "edge_lambda": args.edge_lambda,
        "sigma": args.sigma,
        "patience": args.patience,
        "device": args.device,
    }
    for key, val in cli_overrides.items():
        if val is not None:
            params[key] = val

    denoise(
        input_csv=args.input,
        output_csv=args.output,
        latent_csv=args.latent,
        loss_plot=args.loss_plot,
        **params,
    )


if __name__ == "__main__":
    main()
