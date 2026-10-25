import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F


def initialize_weights(net):
    """Xavier-style weight initialization for Conv1d, ConvTranspose1d, and Linear layers."""
    for m in net.modules():
        if isinstance(m, (nn.Conv1d, nn.ConvTranspose1d)):
            m.weight.data.normal_(0, 0.02)
            if m.bias is not None:
                m.bias.data.zero_()
        elif isinstance(m, nn.Linear):
            m.weight.data.normal_(0, 0.02)
            if m.bias is not None:
                m.bias.data.zero_()


class CAM1D(nn.Module):
    """Channel Attention Module (1D).

    Applies attention along the channel dimension using fully connected layers
    on max-pooled and avg-pooled features.
    """

    def __init__(self, channels, reduction_ratio=16):
        super().__init__()
        self.linear = nn.Sequential(
            nn.Linear(channels, channels // reduction_ratio, bias=True),
            nn.ReLU(inplace=True),
            nn.Linear(channels // reduction_ratio, channels, bias=True),
        )

    def forward(self, x):
        # x: (B, C, L)
        b, c, l = x.size()
        max_pool = F.adaptive_max_pool1d(x, output_size=1).view(b, c)
        avg_pool = F.adaptive_avg_pool1d(x, output_size=1).view(b, c)
        attention = torch.sigmoid(self.linear(max_pool) + self.linear(avg_pool)).unsqueeze(-1)
        return x * attention


class SAM1D(nn.Module):
    """Spatial Attention Module (1D).

    Applies attention along the spatial dimension using Conv1d on
    max-pooled and avg-pooled channels.
    """

    def __init__(self, bias=False):
        super().__init__()
        self.conv = nn.Conv1d(2, 1, kernel_size=7, stride=1, padding=3, bias=bias)

    def forward(self, x):
        max_pool = torch.max(x, dim=1, keepdim=True)[0]
        avg_pool = torch.mean(x, dim=1, keepdim=True)
        concat = torch.cat((max_pool, avg_pool), dim=1)
        attention = self.conv(concat)
        return torch.sigmoid(attention) * x


class CBAM(nn.Module):
    """Convolutional Block Attention Module (1D).

    Combines Channel Attention (CAM1D) and Spatial Attention (SAM1D)
    with a residual skip connection gated by alpha.
    """

    def __init__(self, channels, reduction_ratio=16, alpha=0.1):
        super().__init__()
        self.alpha = alpha
        self.cam = CAM1D(channels=channels, reduction_ratio=reduction_ratio)
        self.sam = SAM1D(bias=False)

    def forward(self, x):
        x_att = self.sam(self.cam(x))
        return x + self.alpha * (x_att - x)


class AE(nn.Module):
    """Attention-Guided Denoising Convolutional AutoEncoder (AG-DCAE).

    Architecture: Encoder (3x Conv1d) -> CBAM attention -> Bottleneck (FC) -> Decoder (3x ConvTranspose1d)

    Args:
        in_dim: Input sequence length (number of genomic bins).
        z_dim: Latent space dimensionality.
        k_size: Kernel size for all convolutional layers.
        alpha: CBAM residual gate strength.
        device: Device string ('cpu' or 'cuda').
    """

    def __init__(self, in_dim, z_dim, k_size, alpha=0.1, device="cuda"):
        super().__init__()
        self.device = device

        # Compute flattened dimension after 3 conv layers (stride=1)
        d = in_dim
        for _ in range(3):
            d = int(np.floor((d - k_size) / 1 + 1))

        self.encoder = nn.Sequential(
            nn.Conv1d(1, 128, k_size, stride=1),
            nn.LeakyReLU(),
            nn.Conv1d(128, 64, k_size, stride=1),
            nn.LeakyReLU(),
            nn.Conv1d(64, 32, k_size, stride=1),
            nn.LeakyReLU(),
        )

        self.cbam = CBAM(channels=32, reduction_ratio=4, alpha=alpha)

        self.fc1 = nn.Linear(d * 32, z_dim)
        self.fc2 = nn.Linear(z_dim, d * 32)

        self.decoder = nn.Sequential(
            nn.LeakyReLU(),
            nn.ConvTranspose1d(32, 64, k_size, stride=1),
            nn.LeakyReLU(),
            nn.ConvTranspose1d(64, 128, k_size, stride=1),
            nn.LeakyReLU(),
            nn.ConvTranspose1d(128, 1, k_size, stride=1),
        )

        initialize_weights(self)

    def bottleneck(self, h):
        return self.fc1(h)

    def encode(self, x):
        h = self.encoder(x)
        h = self.cbam(h)
        h = h.view(h.size(0), -1)
        return self.bottleneck(h)

    def decode(self, z):
        z = self.fc2(z)
        z = z.view(z.size(0), 32, int(z.size(1) / 32))
        z = self.decoder(z)
        return z.squeeze(1)

    def forward(self, x):
        z = self.encode(x)
        y = self.decode(z)
        return z, y
