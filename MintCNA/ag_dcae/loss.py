import torch.nn as nn

criterion = nn.MSELoss()


def edge_preserving_loss(y, x, lam=0.2):
    """Combined reconstruction + edge-preserving loss.

    Preserves sharp CNA boundary transitions during denoising by penalizing
    differences in first derivatives between reconstruction and ground truth.

    Args:
        y: Reconstructed signal, shape (B, L).
        x: Ground truth signal, shape (B, L).
        lam: Weight for the edge-preserving term.

    Returns:
        Loss = MSE(y, x) + lam * MSE(dy, dx)
    """
    rec = criterion(y, x)
    dy = y[:, 1:] - y[:, :-1]
    dx = x[:, 1:] - x[:, :-1]
    edge = criterion(dy, dx)
    return rec + lam * edge
