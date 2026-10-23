import os
import math

import numpy as np
import pandas as pd
import torch
import matplotlib.pyplot as plt

from .model import AE
from .train import xs_gen, train_ag_dcae


def denoise(
    input_csv,
    output_csv,
    latent_csv=None,
    loss_plot=None,
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
):
    """Run AG-DCAE denoising: train model and produce denoised output.

    Reads a CSV of log2-ratio profiles (bins x cells), trains the AG-DCAE,
    and writes the denoised matrix back as CSV.

    Args:
        input_csv: Path to input CSV (bins x cells, first column = bin names).
        output_csv: Path to write denoised CSV (bins x cells).
        latent_csv: Optional path to write latent embeddings (cells x latent_dim).
        loss_plot: Optional path to save training loss curve plot.
        latent_dim: Latent space dimensionality.
        kernel_size: Conv1d kernel size.
        lr: Learning rate.
        epochs: Maximum training epochs.
        batch_size: Samples per batch.
        alpha: CBAM residual gate strength.
        edge_lambda: Edge-preserving loss weight.
        sigma: Gaussian noise std for training.
        patience: Early stopping patience.
        min_delta: Early stopping min improvement.
        test_size: Train/test split fraction.
        random_state: Random seed.
        device: Torch device string ('cpu' or 'cuda').
    """
    # Load data: CSV is bins x cells
    data_df = pd.read_csv(input_csv, index_col=0)
    data = data_df.values.T  # Transpose to cells x bins for model

    if device is None:
        device_str = "cuda" if torch.cuda.is_available() else "cpu"
    else:
        device_str = device

    # Train model
    model, train_losses, test_losses = train_ag_dcae(
        data,
        latent_dim=latent_dim,
        kernel_size=kernel_size,
        lr=lr,
        epochs=epochs,
        batch_size=batch_size,
        alpha=alpha,
        edge_lambda=edge_lambda,
        sigma=sigma,
        patience=patience,
        min_delta=min_delta,
        test_size=test_size,
        random_state=random_state,
        device=device_str,
    )

    # Inference on full dataset
    torch_device = torch.device(device_str)
    model.eval()
    zs = []
    x_denoised = []

    for _, x in xs_gen(data, batch_size, shuffle=False):
        x = torch.from_numpy(x).float().to(torch_device).unsqueeze(1)
        with torch.no_grad():
            z, y = model(x)
            zs.append(z.cpu().numpy())
            x_denoised.append(y.cpu().numpy())

    Z = np.concatenate(zs, axis=0)
    x_denoised = np.concatenate(x_denoised, axis=0)

    # Convert back to bins x cells DataFrame with original names
    x_denoised_df = pd.DataFrame(x_denoised.T, index=data_df.index, columns=data_df.columns)

    # Save outputs
    os.makedirs(os.path.dirname(os.path.abspath(output_csv)), exist_ok=True)
    x_denoised_df.to_csv(output_csv, index=True)
    print(f"Denoised matrix saved: {output_csv} (shape {x_denoised_df.shape})")

    if latent_csv is not None:
        z_df = pd.DataFrame(Z, index=data_df.columns)
        os.makedirs(os.path.dirname(os.path.abspath(latent_csv)), exist_ok=True)
        z_df.to_csv(latent_csv, index=True)
        print(f"Latent embeddings saved: {latent_csv} (shape {z_df.shape})")

    # Save loss plot
    if loss_plot is not None:
        os.makedirs(os.path.dirname(os.path.abspath(loss_plot)), exist_ok=True)
        epoch_list = list(range(1, len(train_losses) + 1))
        plt.figure(figsize=(8, 6))
        plt.plot(epoch_list, train_losses, label="Training Loss", marker="o", linestyle="--", color="b", markersize=4)
        plt.plot(epoch_list, test_losses, label="Test Loss", marker="s", linestyle="-", color="r", markersize=4)
        plt.xlabel("Epochs")
        plt.ylabel("Loss")
        plt.title("AG-DCAE Training and Test Loss")
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.savefig(loss_plot, format="pdf", bbox_inches="tight")
        plt.close()
        print(f"Loss curve saved: {loss_plot}")

    # Sanity check
    cell = data_df.columns[0]
    corr = np.corrcoef(data_df[cell].values, x_denoised_df[cell].values)[0, 1]
    print(f"Sanity check — correlation for {cell}: {corr:.4f}")
