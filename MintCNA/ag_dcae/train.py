import math
import time
from copy import deepcopy

import numpy as np
import torch
from sklearn.model_selection import train_test_split

from .model import AE
from .loss import edge_preserving_loss


def xs_gen(X, batch_size, shuffle=True):
    """Generator that yields batches of data.

    Args:
        X: Numpy array of shape (n_samples, n_features).
        batch_size: Number of samples per batch.
        shuffle: Whether to shuffle indices each epoch.

    Yields:
        (step_index, batch_array) tuples.
    """
    n = X.shape[0]
    idx = np.random.permutation(n) if shuffle else np.arange(n)
    steps = math.ceil(n / batch_size)
    for i in range(steps):
        sl = idx[i * batch_size : (i + 1) * batch_size]
        yield i, X[sl]


def train_ag_dcae(
    data,
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
    """Train the AG-DCAE model.

    Args:
        data: Numpy array of shape (n_cells, n_bins) — log2-ratio profiles.
        latent_dim: Latent space dimensionality.
        kernel_size: Conv1d kernel size.
        lr: Learning rate for Adam optimizer.
        epochs: Maximum training epochs.
        batch_size: Samples per batch.
        alpha: CBAM residual gate strength.
        edge_lambda: Edge-preserving loss weight.
        sigma: Gaussian noise std injected during training.
        patience: Early stopping patience (epochs without improvement).
        min_delta: Minimum improvement to reset patience counter.
        test_size: Fraction of data for validation.
        random_state: Random seed (train/test split, weight init, shuffling, noise).
        device: Torch device. Auto-detected if None.

    Returns:
        model: Trained AE model (best weights restored).
        train_losses: List of per-epoch training losses.
        test_losses: List of per-epoch test losses.
    """
    if device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    else:
        device = torch.device(device)

    # Seed everything the training touches -- weight init (torch), the per-epoch
    # shuffle and noise injection (numpy / torch) -- so a config reproduces bit-for-bit.
    np.random.seed(random_state)
    torch.manual_seed(random_state)
    if device.type == "cuda":
        torch.cuda.manual_seed_all(random_state)

    in_dim = data.shape[1]
    model = AE(in_dim=in_dim, z_dim=latent_dim, k_size=kernel_size, alpha=alpha, device=str(device)).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)

    train_data, test_data = train_test_split(data, test_size=test_size, random_state=random_state)

    train_losses = []
    test_losses = []
    best_loss = float("inf")
    bad_epochs = 0
    best_state = None

    for epoch in range(epochs):
        start_time = time.time()

        # Training phase
        model.train()
        total_train_loss = 0
        for _, x in xs_gen(train_data, batch_size, shuffle=True):
            x = torch.from_numpy(x).float().to(device).unsqueeze(1)
            x_noisy = x + sigma * torch.randn_like(x)
            _, y = model(x_noisy)
            loss = edge_preserving_loss(y, x.squeeze(1), lam=edge_lambda)
            total_train_loss += loss.item()
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

        # Validation phase
        model.eval()
        total_test_loss = 0
        with torch.no_grad():
            for _, x in xs_gen(test_data, batch_size, shuffle=False):
                x = torch.from_numpy(x).float().to(device).unsqueeze(1)
                _, y = model(x)
                loss = edge_preserving_loss(y, x.squeeze(1), lam=edge_lambda)
                total_test_loss += loss.item()

        train_losses.append(total_train_loss)
        test_losses.append(total_test_loss)
        elapsed = time.time() - start_time
        print(f"epoch: {epoch+1} | train loss: {total_train_loss:.4f} | test loss: {total_test_loss:.4f} | time: {elapsed:.2f}s")

        # Early stopping
        if total_test_loss < best_loss - min_delta:
            best_loss = total_test_loss
            bad_epochs = 0
            best_state = deepcopy(model.state_dict())
        else:
            bad_epochs += 1
            if bad_epochs >= patience:
                print(f"Early stopping at epoch {epoch+1}: no improvement in {patience} epochs (best test loss={best_loss:.4f})")
                break

    if best_state is not None:
        model.load_state_dict(best_state)
        print("Best model restored (lowest test loss)")

    return model, train_losses, test_losses
