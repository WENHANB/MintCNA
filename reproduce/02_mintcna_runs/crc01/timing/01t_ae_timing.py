# CRC01 re-run 2026-10-14: production AG-DCAE (Code/MultiCNA/CAE_SCAtt_DNAmeth.py) with latent_dim = 250,
# on the 197 QC-passing cells. Outputs carry the suffix _d250_197; nothing published is overwritten.
#!/usr/bin/env python
# coding: utf-8

# # Convolutional Autoencoder(CAE) with Spatial and Channel attention layers

# In[1]:


import numpy as np
import math
import torch
import pandas as pd
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
import time 
import torch.nn.functional as F
import os
from copy import deepcopy
#os.chdir('<MINTCNA_ROOT>/Data/CRC/DNAmeth/processed_l2r')
os.chdir('<MINTCNA_ROOT>/Data/CRC/DNAmeth/processed_l2r')
print(os.getcwd())


# In[2]:


# Function to initialize weights
def initialize_weights(net):
    for m in net.modules():
        if isinstance(m, nn.Conv1d) or isinstance(m, nn.ConvTranspose1d):
            m.weight.data.normal_(0, 0.02)
            if m.bias is not None:
                m.bias.data.zero_()
        elif isinstance(m, nn.Linear):
            m.weight.data.normal_(0, 0.02)
            if m.bias is not None:
                m.bias.data.zero_()
                
criterion = nn.MSELoss()

def edge_preserving_loss(y, x, lam=0.2):
    # y, x: (B, L)
    rec = criterion(y, x)
    dy = y[:, 1:] - y[:, :-1]
    dx = x[:, 1:] - x[:, :-1]
    edge = criterion(dy, dx)
    return rec + lam * edge

# Generator function to yield batches of data
def xs_gen(X, batch_size, shuffle=True):
    n = X.shape[0]
    idx = np.random.permutation(n) if shuffle else np.arange(n)
    steps = math.ceil(n / batch_size)
    for i in range(steps):
        sl = idx[i*batch_size:(i+1)*batch_size]
        yield i, X[sl]

class CBAM(nn.Module):
    def __init__(self, channels, reduction_ratio=16, alpha=0.1):
        super(CBAM, self).__init__()
        self.channels = channels
        self.reduction_ratio = reduction_ratio
        self.alpha = alpha
        self.cam = CAM1D(channels=self.channels, reduction_ratio=self.reduction_ratio)  # Channel Attention
        self.sam = SAM1D(bias=False)  # Spatial Attention

    def forward(self, x):
        x_att = self.sam(self.cam(x))
        return x + self.alpha * (x_att - x)

# Channel Attention Module
class CAM1D(nn.Module):
    def __init__(self, channels, reduction_ratio=16):
        super(CAM1D, self).__init__()
        self.channels = channels
        self.reduction_ratio = reduction_ratio
        
        # Fully connected layers for channel attention
        self.linear = nn.Sequential(
            nn.Linear(in_features=self.channels, out_features=self.channels // self.reduction_ratio, bias=True),
            nn.ReLU(inplace=True),
            nn.Linear(in_features=self.channels // self.reduction_ratio, out_features=self.channels, bias=True)
        )

    def forward(self, x):
        # Input shape: (B, C, L)
        b, c, l = x.size()

        # Global Max Pooling and Avg Pooling along the spatial dimension (L)
        max_pool = F.adaptive_max_pool1d(x, output_size=1).view(b, c)  # Shape: (B, C)
        avg_pool = F.adaptive_avg_pool1d(x, output_size=1).view(b, c)  # Shape: (B, C)

        # Fully connected layers for attention
        linear_max = self.linear(max_pool)  # Shape: (B, C)
        linear_avg = self.linear(avg_pool)  # Shape: (B, C)

        # Combine and apply Sigmoid
        attention = torch.sigmoid(linear_max + linear_avg).unsqueeze(-1)  # Shape: (B, C, 1)

        # Broadcast attention weights along the sequence dimension (L)
        output = x * attention  # Element-wise multiplication
        return output
        
# Spatial Attention Module
class SAM1D(nn.Module):
    def __init__(self, bias=False):
        super(SAM1D, self).__init__()
        self.bias = bias
        # 1D Convolution: Replace Conv2d with Conv1d
        self.conv = nn.Conv1d(in_channels=2, out_channels=1, kernel_size=7, stride=1, padding=3, bias=self.bias)

    def forward(self, x):
        # Apply max pooling along the channel dimension
        max_pool = torch.max(x, dim=1, keepdim=True)[0]  # Shape: (B, 1, L)
        # Apply average pooling along the channel dimension
        avg_pool = torch.mean(x, dim=1, keepdim=True)  # Shape: (B, 1, L)
        # Concatenate max and avg pooling results along the channel dimension
        concat = torch.cat((max_pool, avg_pool), dim=1)  # Shape: (B, 2, L)
        # Apply 1D convolution to learn attention weights
        attention = self.conv(concat)  # Shape: (B, 1, L)
        # Sigmoid activation and element-wise multiplication with input
        output = torch.sigmoid(attention) * x
        return output

# Autoencoder Model
class AE(torch.nn.Module):
    def __init__(self, in_dim, z_dim, k_size, alpha=0.1,device='cuda'):
        super(AE, self).__init__()
        self.device = device

        d = in_dim
        for i in range(3):
            d = np.floor((d - k_size) / 1 + 1)
        d = np.int32(d)

        self.encoder = nn.Sequential(
            nn.Conv1d(1, 128, k_size, stride=1),
            nn.LeakyReLU(),
            nn.Conv1d(128, 64, k_size, stride=1),
            nn.LeakyReLU(),
            nn.Conv1d(64, 32, k_size, stride=1),
            nn.LeakyReLU()
        )

        # Attention Mechanism
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
        z = self.fc1(h)
        return z

    def encode(self, x):
        h = self.encoder(x)
        h = self.cbam(h)
        h = h.view(h.size(0), -1)  # Flatten the output
        z = self.bottleneck(h)
        return z

    def decode(self, z):
        z = self.fc2(z)
        z = z.view(z.size(0), 32, int(z.size(1) / 32))  # Reshape the output
        z = self.decoder(z)
        z = z.squeeze(1)     # Squeeze to ensure correct shape (samples, markers)
        return z

    def forward(self, x):
        z = self.encode(x)
        y = self.decode(z)
        return z, y


# In[3]:


data_df = pd.read_csv('./l2r_smoothed_all_cells_100kb.csv',index_col=0)
data_df = data_df.drop(columns=['CRC01_MP1_209'])          # 197 QC-passing cells only (2026-10-14)
assert data_df.shape[1] == 197, data_df.shape
torch.manual_seed(30); np.random.seed(30)                   # reproducible re-run
data    = data_df.values.T  
print(data.shape)

class Args:
    latent_dim = 250         # Latent space size
    kernel_size = 5        # Kernel size for convolution (as in CAE)
    lr = 0.0001            # Learning rate
    epochs = 100           # Number of epochs to train
    batch_size = 16       # Batch size for training
    alpha = 0.1   # CBAM gate strength
    lam   = 0.2   # edge loss weight

args = Args()
sigma = 0.08

# Define and create AE architecture
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
ae_model = AE(in_dim=data.shape[1], z_dim=args.latent_dim, k_size=args.kernel_size, alpha=args.alpha).to(device)
optimizer = torch.optim.Adam(ae_model.parameters(), lr=args.lr)
# Define learning rate scheduler (drops LR every 30 epochs by a factor of 0.1)
#scheduler = optim.lr_scheduler.StepLR(optimizer, step_size=40, gamma=0.1)
print("has alpha?", hasattr(ae_model.cbam, "alpha"))
print("cbam:", ae_model.cbam)

# Split the data into training and test sets
train_data, test_data = train_test_split(data, test_size=0.2, random_state=30)


# In[4]:


# ---- training with simple early stop + save best model ----
epochs = []
train_loss = []
test_loss = []

patience = 15          # how many epochs to wait after no improvement
min_delta = 0.001      # minimum improvement to count
best_loss = float('inf')
bad_epochs = 0
best_state = None      # store the best model weights

for epoch in range(args.epochs):
    start_time = time.time()
    print("epoch", epoch + 1)
    ae_model.train()
    total_train_loss = 0

    for step, x in xs_gen(train_data, args.batch_size, shuffle=True):
        x = torch.from_numpy(x).float().to(device).unsqueeze(1)
        x_noisy = x + sigma * torch.randn_like(x)
        _, y = ae_model(x_noisy)
        loss = edge_preserving_loss(y, x.squeeze(1), lam=args.lam)
        total_train_loss += loss.item()
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

    ae_model.eval()
    total_test_loss = 0
    with torch.no_grad():
        for step, x in xs_gen(test_data, args.batch_size, shuffle=False):
            x = torch.from_numpy(x).float().to(device).unsqueeze(1)
            _, y = ae_model(x)
            loss = edge_preserving_loss(y, x.squeeze(1), lam=args.lam)
            total_test_loss += loss.item()

    train_loss.append(total_train_loss)
    test_loss.append(total_test_loss)
    epoch_duration = time.time() - start_time
    print(f"epoch: {epoch+1} | train loss: {total_train_loss:.4f} | test loss: {total_test_loss:.4f} | time: {epoch_duration:.2f}s")

    # ----- early stopping & save best -----
    if total_test_loss < best_loss - min_delta:
        best_loss = total_test_loss
        bad_epochs = 0
        best_state = deepcopy(ae_model.state_dict())  # save best weights
    else:
        bad_epochs += 1
        if bad_epochs >= patience:
            print(f"Early stopping at epoch {epoch+1}: no improvement in {patience} epochs (best test loss={best_loss:.4f})")
            break

# restore best model at the end
if best_state is not None:
    ae_model.load_state_dict(best_state)
    print("Best model restored (lowest test loss)")

# ---------------------------
#  SAVE DENOISED DNA OUTPUT
# ---------------------------

ae_model.eval()

zs = []
x_denoised = []   # store outputs batch by batch

# ---- run model on each batch ----
for step, x in xs_gen(data, args.batch_size, shuffle=False):   # data = numpy (cells × bins)
    x = torch.from_numpy(x).float().to(device).unsqueeze(1)  # [B,1,L]

    with torch.no_grad():
        z, y = ae_model(x)       # y shape: [B, L]
        zs.append(z.detach().cpu().numpy())
        x_denoised.append(y.cpu().numpy())

# ---- combine all batches ----
Z = np.concatenate(zs, axis=0)             # (cells, latent_dim)
x_denoised = np.concatenate(x_denoised, axis=0)   # (cells, bins)

# ---- convert to DataFrame: bins × cells, KEEP row + col names ----
# data_df is the original DataFrame you loaded (bins × cells)
x_denoised_df = pd.DataFrame(
    x_denoised.T,                # (bins × cells)
    index=data_df.index,         # bin / region names
    columns=data_df.columns      # cell IDs
)
z_df = pd.DataFrame(Z, index=data_df.columns)   # rows=cells

print("Denoised DF shape:", x_denoised_df.shape)
print("Denoised DF shape:", z_df.shape)

# ---- save WITH index (row names) ----
x_denoised_df.to_csv("<REPRODUCE_DIR>/02_mintcna_runs/crc01/timing/ae_dna_d250_197_turin.csv", index=True)
z_df.to_csv("<REPRODUCE_DIR>/02_mintcna_runs/crc01/timing/latent_z_dna_d250_197_turin.csv", index=True)


# In[5]:


# Loss plot
epochs = list(range(1, len(train_loss) + 1))
plt.figure(figsize=(8, 6))
plt.plot(epochs, train_loss, label="Training Loss", marker='o', linestyle='--', color='b', markersize=4)
plt.plot(epochs, test_loss, label="Test Loss", marker='s', linestyle='-', color='r', markersize=4)
plt.xlabel('Epochs')
plt.ylabel('Loss')
plt.title('Training and Test Loss Over Epochs')
plt.grid(True)
plt.legend()

plt.tight_layout()
os.makedirs("Result", exist_ok=True)
plt.savefig("<REPRODUCE_DIR>/02_mintcna_runs/crc01/timing/loss_curve_turin.pdf", format="pdf", bbox_inches="tight")  # saves to current working directory
plt.show()


# In[6]:


# sanity check one named cell matches
cell = data_df.columns[0]
x_orig = data_df[cell].values
x_den  = x_denoised_df[cell].values
print("corr saved:", np.corrcoef(x_orig, x_den)[0,1])
