# -*- coding: utf-8 -*-
"""
Created on Fri Aug 15 14:05:03 2025

@author: tracy
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


MIN_TIME = 80.0  # seconds
MAX_TIME = 110.0  # seconds


# Load CSV, skip header row
print("Loading pure data...")
df = pd.read_csv("pure_data.csv", skiprows=1, header=None)

# Remove every other line if it's NaN (blank row)
print("Cleaning data...")
df = df.dropna(how='all')

# Drop unwanted columns (node label = col 1, z = col 4)
df = df.drop(columns=[1, 4])

# Reset index
df = df.reset_index(drop=True)

# Rename columns just for code clarity
df.columns = ["time", "x", "y", "temperature"]


# ===== Optional: Plot temperature vs time for a single node =====
print("Plotting temperature vs. time for node with highest average temperature...")

# Compute average temperature per node
avg_temp_per_node = df.groupby(["x", "y"])["temperature"].mean()

# Select node with highest average temperature
node_max_avg = avg_temp_per_node.idxmax()
print(
    f"Selected node (x, y) = {node_max_avg} with avg temp = {avg_temp_per_node[node_max_avg]:.2f}")

# Extract time series for this node
node_df = df[(df["x"] == node_max_avg[0]) & (
    df["y"] == node_max_avg[1])].sort_values("time")

# Plot
plt.figure(figsize=(6, 4))
plt.plot(node_df["time"], node_df["temperature"],
         marker='o', color='tab:red', markersize=1.0)
plt.xlabel("Time (s)")
plt.ylabel("Temperature")
plt.title(f"Temperature vs Time at Node {node_max_avg}")
plt.grid(True)
plt.tight_layout()
plt.show()


# ===== STEP 2: Filter time window =====
df = df[(df["time"] >= MIN_TIME) & (
    df["time"] <= MAX_TIME)].reset_index(drop=True)


# ===== STEP 3: Group by node and FFT =====
print("Performing FFT...")
grouped = df.groupby(["x", "y"])

fft_results = {}  # (x, y) -> {freqs, phase}

for (x, y), group in grouped:
    group = group.sort_values("time")
    t = group["time"].values
    T = group["temperature"].values

    fft_vals = np.fft.fft(T)
    freqs = np.fft.fftfreq(len(t), d=(t[1] - t[0]))

    mask = freqs > 0
    freqs = freqs[mask]
    fft_vals = fft_vals[mask]

    phase = np.angle(fft_vals)  # radians

    fft_results[(x, y)] = {
        "freqs": freqs,
        "phase": phase
    }


# ===== STEP 4: Build frequency-grouped DataFrame =====
any_node = next(iter(fft_results))
freqs = fft_results[any_node]["freqs"]

phase_list = []

for f_idx, f in enumerate(freqs):
    for (x, y), data in fft_results.items():
        p = data["phase"][f_idx]
        phase_list.append([f, x, y, p])

fft_df = pd.DataFrame(phase_list, columns=["frequency", "x", "y", "phase"])
fft_df = fft_df.sort_values("frequency").reset_index(drop=True)


# ===== STEP 5: Phase maps with baseline subtraction + normalization =====
print("Plotting phase maps...")

# Define custom colormap: black -> orange -> yellow
colors = ["black", "red", "orange", "yellow", "white"]
custom_cmap = LinearSegmentedColormap.from_list(
    "black_orange_yellow", colors)

freqs = fft_df["frequency"].unique()

counter = 1

for f in freqs:
    freq_df = fft_df[fft_df["frequency"] == f].copy()

    # Subtract matrix node baseline
    # Here we take the max phase node as the "matrix"
    baseline = freq_df["phase"].max()
    freq_df["phase_corrected"] = - \
        (freq_df["phase"] - baseline)  # multiply by -1

    # Normalize 0 -> 1
    phase_min = freq_df["phase_corrected"].min()
    phase_max = freq_df["phase_corrected"].max()
    freq_df["phase_norm"] = 1 -\
        (freq_df["phase_corrected"] - phase_min) / (phase_max - phase_min)

    # Pivot to 2D grid
    phase_grid = freq_df.pivot(index="y", columns="x", values="phase_norm")

    # Plot
    plt.figure(figsize=(6, 5))
    plt.imshow(phase_grid, origin='lower', aspect='auto',
               cmap=custom_cmap, vmin=0, vmax=1)
    plt.colorbar(label='Normalized Phase')
    plt.title(f'Phase Map at {f:.2f} Hz')
    plt.xlabel('x')
    plt.ylabel('y')

    plt.tight_layout()
    plt.savefig(f'phase_map_{f:.2f}Hz.png', dpi=300)
    plt.show()

    print(f"Complete: Map {counter} of {len(freqs)}")

    counter += 1

# ===== STEP 6: Calculate SNR for each phase map using top/bottom percentiles =====
PERCENTILE = 20  # top/bottom X% to define defect/background

print("\nCalculating SNR for each phase map...")

snr_results = []

for f in freqs:
    freq_df = fft_df[fft_df["frequency"] == f].copy()

    # Subtract baseline and normalize as before
    baseline = freq_df["phase"].max()
    freq_df["phase_corrected"] = -(freq_df["phase"] - baseline)
    phase_min = freq_df["phase_corrected"].min()
    phase_max = freq_df["phase_corrected"].max()
    freq_df["phase_norm"] = 1 - \
        (freq_df["phase_corrected"] - phase_min) / (phase_max - phase_min)

    # Flatten normalized phase values
    phase_vals = freq_df["phase_norm"].values

    # Define defect and background regions
    top_threshold = np.percentile(phase_vals, 100 - PERCENTILE)  # defect
    bottom_threshold = np.percentile(phase_vals, PERCENTILE)     # background

    defect_vals = phase_vals[phase_vals >= top_threshold]
    background_vals = phase_vals[phase_vals <= bottom_threshold]

    # Compute SNR
    mean_defect = np.mean(defect_vals)
    mean_background = np.mean(background_vals)
    std_background = np.std(background_vals)

    snr = (mean_defect - mean_background) / std_background
    snr_results.append((f, snr))

    print(f"Frequency {f:.2f} Hz -> SNR: {snr:.3f}")

# Optional: plot SNR vs frequency
freqs_list, snr_list = zip(*snr_results)
plt.figure(figsize=(6, 4))
plt.plot(freqs_list, snr_list, marker='o', color='tab:blue')
plt.xlabel("Frequency (Hz)")
plt.ylabel(f"SNR (top/bottom {PERCENTILE}%)")
plt.title("SNR vs Frequency")
plt.grid(True)
plt.tight_layout()
plt.show()
