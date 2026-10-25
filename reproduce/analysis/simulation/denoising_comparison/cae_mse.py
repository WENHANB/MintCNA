"""CAE + MSE baseline for Supplementary Fig. 1: the AG-DCAE encoder/decoder without the CBAM attention
module, trained with a plain mean-squared-error reconstruction loss (no edge-preserving term)."""
import math
import time
from copy import deepcopy

import numpy as np
import torch
import torch.nn as nn
from sklearn.model_selection import train_test_split

from ag_dcae.model import initialize_weights   # MintCNA/ag_dcae (run from the repository root)


class CAE_MSE(nn.Module):
    def __init__(self, in_dim, z_dim, k_size):
        super().__init__()
        d = in_dim
        for _ in range(3):
            d = int(np.floor((d - k_size) / 1 + 1))
        self.encoder = nn.Sequential(
            nn.Conv1d(1, 128, k_size, stride=1), nn.LeakyReLU(),
            nn.Conv1d(128, 64, k_size, stride=1), nn.LeakyReLU(),
            nn.Conv1d(64, 32, k_size, stride=1), nn.LeakyReLU(),
        )
        self.fc1 = nn.Linear(d * 32, z_dim)
        self.fc2 = nn.Linear(z_dim, d * 32)
        self.decoder = nn.Sequential(
            nn.LeakyReLU(),
            nn.ConvTranspose1d(32, 64, k_size, stride=1), nn.LeakyReLU(),
            nn.ConvTranspose1d(64, 128, k_size, stride=1), nn.LeakyReLU(),
            nn.ConvTranspose1d(128, 1, k_size, stride=1),
        )
        initialize_weights(self)

    def encode(self, x):
        h = self.encoder(x)
        return self.fc1(h.view(h.size(0), -1))

    def decode(self, z):
        z = self.fc2(z)
        z = z.view(z.size(0), 32, z.size(1) // 32)
        return self.decoder(z).squeeze(1)

    def forward(self, x):
        z = self.encode(x)
        return z, self.decode(z)


def _xs_gen(X, batch_size, shuffle=True):
    n = X.shape[0]
    idx = np.random.permutation(n) if shuffle else np.arange(n)
    for i in range(math.ceil(n / batch_size)):
        yield i, X[idx[i * batch_size:(i + 1) * batch_size]]


def train_cae_mse(data, latent_dim=200, kernel_size=5, lr=1e-4, epochs=100, batch_size=16, sigma=0.08,
                  patience=15, min_delta=0.001, test_size=0.2, random_state=30, device=None):
    """Denoising training (Gaussian input noise of s.d. sigma), early stopping on the held-out loss."""
    device = torch.device(device) if device is not None else torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = CAE_MSE(in_dim=data.shape[1], z_dim=latent_dim, k_size=kernel_size).to(device)
    optim = torch.optim.Adam(model.parameters(), lr=lr)
    mse = nn.MSELoss()
    train_data, test_data = train_test_split(data, test_size=test_size, random_state=random_state)

    train_losses, test_losses = [], []
    best_loss, bad, best_state = float("inf"), 0, None
    for ep in range(epochs):
        t0 = time.time()
        model.train()
        tr = 0.0
        for _, x in _xs_gen(train_data, batch_size, shuffle=True):
            x = torch.from_numpy(x).float().to(device).unsqueeze(1)
            _, y = model(x + sigma * torch.randn_like(x))
            loss = mse(y, x.squeeze(1))
            optim.zero_grad()
            loss.backward()
            optim.step()
            tr += loss.item()
        model.eval()
        te = 0.0
        with torch.no_grad():
            for _, x in _xs_gen(test_data, batch_size, shuffle=False):
                x = torch.from_numpy(x).float().to(device).unsqueeze(1)
                _, y = model(x)
                te += mse(y, x.squeeze(1)).item()
        train_losses.append(tr)
        test_losses.append(te)
        print(f"[CAE+MSE] ep {ep+1}: train={tr:.4f}  test={te:.4f}  ({time.time()-t0:.1f}s)")
        if te < best_loss - min_delta:
            best_loss, bad, best_state = te, 0, deepcopy(model.state_dict())
        else:
            bad += 1
            if bad >= patience:
                print(f"[CAE+MSE] early stop at ep {ep+1} (best test={best_loss:.4f})")
                break
    if best_state is not None:
        model.load_state_dict(best_state)
    return model, train_losses, test_losses
