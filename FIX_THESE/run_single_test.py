#!/usr/bin/env python3
"""Single rank-1 test on generated data."""
import sys
import os
import json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from sim_generate import generate_dataset
from sim_fit_v3 import run_calibration_trial

# Load di-charmonium config
with open('tests_top3.json') as f:
    config = json.load(f)

dicharm_config = config['tests'][2]  # Di-charmonium
print(f"Test: {dicharm_config['name']}")
print(f"Channel A: {dicharm_config['channelA']['name']}")
print(f"Channel B: {dicharm_config['channelB']['name']}")
print()

# Generate M0 dataset (rank-1 true)
print("Generating M0 (rank-1 true) dataset...")
dataset = generate_dataset(dicharm_config, 'M0', scale_factor=1.0, seed=42)
print(f"  Channel A bins: {len(dataset['channelA']['y'])}")
print(f"  Channel B bins: {len(dataset['channelB']['y'])}")
print()

# Run single fit trial
print("Running rank-1 test (100 bootstrap, 60 starts)...")
result = run_calibration_trial(dataset, n_bootstrap=100, n_starts=60)

print()
print("=" * 50)
print("RESULT")
print("=" * 50)
print(f"Full result: {result}")
lam = result.get('Lambda_obs')
pval = result.get('p_boot')
print(f"Lambda_obs: {lam:.4f}" if lam is not None else "Lambda_obs: N/A")
print(f"p_boot: {pval:.4f}" if pval is not None else "p_boot: N/A")
print(f"Pass gates: {result.get('pass_trial', 'N/A')}")
print(f"Rejected (p<0.05): {result.get('rejected', 'N/A')}")
print()
if result.get('pass_trial'):
    if result.get('rejected'):
        print("VERDICT: Rank-1 REJECTED (false positive under M0)")
    else:
        print("VERDICT: Rank-1 NOT rejected (correct under M0)")
else:
    print("VERDICT: Fit failed gates")
