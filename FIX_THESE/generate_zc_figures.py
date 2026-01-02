#!/usr/bin/env python3
"""Generate figures for Zc-like rank-1 test."""
import json
import numpy as np
import matplotlib.pyplot as plt
from sim_generate import generate_dataset
from sim_fit_v3 import run_calibration_trial

# Load config
with open('tests_top3.json') as f:
    config = json.load(f)

zc_config = config['tests'][1]  # Zc-like
print(f"Test: {zc_config['name']}")
print(f"Channel A: {zc_config['channelA']['name']}")
print(f"Channel B: {zc_config['channelB']['name']}")

# Generate M0 dataset
print("\nGenerating M0 dataset...")
dataset = generate_dataset(zc_config, 'M0', scale_factor=1.0, seed=456)

# Run fit
print("Running rank-1 test (100 bootstrap, 60 starts)...")
result = run_calibration_trial(dataset, n_bootstrap=100, n_starts=60)

lambda_obs = result.get('lambda_obs', 0)
p_boot = result.get('p_boot', 0)
bootstrap_lambdas = result.get('bootstrap_lambdas', [])

print(f"\nLambda_obs: {lambda_obs:.4f}")
print(f"p_boot: {p_boot:.4f}")

# Figure 1: Bootstrap Lambda Distribution
fig, ax = plt.subplots(figsize=(8, 5))
if bootstrap_lambdas:
    ax.hist(bootstrap_lambdas, bins=30, density=True, alpha=0.7, color='steelblue', edgecolor='black')
    ax.axvline(lambda_obs, color='red', linestyle='--', linewidth=2, label=f'Observed Λ = {lambda_obs:.2f}')
    ax.axvline(5.99, color='orange', linestyle=':', linewidth=2, label='χ²(2) 95% = 5.99')
    ax.set_xlabel('Lambda (Λ)', fontsize=12)
    ax.set_ylabel('Density', fontsize=12)
    ax.set_title(f'Zc-like Bootstrap Distribution (p = {p_boot:.2f})', fontsize=14)
    ax.legend()
    ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('out/zc_bootstrap_hist.png', dpi=150)
print("Saved: out/zc_bootstrap_hist.png")

# Figure 2: Channel spectra
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Channel A
chA = dataset['channelA']
axes[0].errorbar(chA['x'], chA['y'], yerr=chA['sigma'], fmt='o', markersize=4, capsize=2, color='blue')
axes[0].set_xlabel('Mass (GeV)', fontsize=12)
axes[0].set_ylabel('Counts', fontsize=12)
axes[0].set_title(f"Channel A: {zc_config['channelA']['name']} (πJ/ψ)", fontsize=12)
axes[0].grid(True, alpha=0.3)

# Channel B
chB = dataset['channelB']
axes[1].errorbar(chB['x'], chB['y'], yerr=chB['sigma'], fmt='o', markersize=4, capsize=2, color='green')
axes[1].set_xlabel('Mass (GeV)', fontsize=12)
axes[1].set_ylabel('Counts', fontsize=12)
axes[1].set_title(f"Channel B: {zc_config['channelB']['name']} (DD*)", fontsize=12)
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('out/zc_channel_spectra.png', dpi=150)
print("Saved: out/zc_channel_spectra.png")

print("\nDone!")
