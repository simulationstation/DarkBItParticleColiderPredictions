#!/usr/bin/env python3
"""
BESIII Y(4220)/Y(4320) Rank-1 Factorization Test

Cross-channel test: pi+pi- J/psi vs pi+pi- hc

Tests whether the complex coupling ratio R = g(Y4320)/g(Y4220) is shared
between the two decay channels.
"""

import argparse
import csv
import json
import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2
from pathlib import Path
from multiprocessing import Pool, cpu_count
import warnings
warnings.filterwarnings('ignore')

# Physics constants (from arXiv:1611.01317 and HEPData)
Y4220_MASS = 4.222   # GeV
Y4220_WIDTH = 0.044  # GeV
Y4320_MASS = 4.320   # GeV
Y4320_WIDTH = 0.101  # GeV

# Fit window
E_MIN = 4.15
E_MAX = 4.45


def load_channel_data(filepath):
    """Load cross section data from CSV."""
    data = []
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split(',')
            if len(parts) >= 3:
                sqrt_s = float(parts[0])
                sigma = float(parts[1])
                stat_err = float(parts[2])
                sys_err = float(parts[3]) if len(parts) > 3 else 0
                # Total error (add in quadrature)
                total_err = np.sqrt(stat_err**2 + sys_err**2)
                if E_MIN <= sqrt_s <= E_MAX and sigma > 0:
                    data.append({
                        'sqrt_s': sqrt_s,
                        'sigma': sigma,
                        'err': max(total_err, 0.1)  # Floor to prevent zero errors
                    })
    return data


def breit_wigner(s, M, Gamma):
    """Relativistic Breit-Wigner amplitude (complex)."""
    return M * Gamma / (s - M**2 + 1j * M * Gamma)


def cross_section_model(sqrt_s, params, channel='A'):
    """
    Two-resonance cross section model.

    params for unconstrained fit:
        [norm_A, r_A, phi_A, norm_B, r_B, phi_B, M1, G1, M2, G2]

    params for constrained fit (shared R):
        [norm_A, norm_B, r_shared, phi_shared, M1, G1, M2, G2]
    """
    s = sqrt_s**2

    if len(params) == 10:
        # Unconstrained: different R per channel
        norm_A, r_A, phi_A, norm_B, r_B, phi_B, M1, G1, M2, G2 = params
        if channel == 'A':
            norm, r, phi = norm_A, r_A, phi_A
        else:
            norm, r, phi = norm_B, r_B, phi_B
    else:
        # Constrained: shared R
        norm_A, norm_B, r_shared, phi_shared, M1, G1, M2, G2 = params
        r, phi = r_shared, phi_shared
        norm = norm_A if channel == 'A' else norm_B

    # Complex amplitudes
    A1 = breit_wigner(s, M1, G1)
    A2 = breit_wigner(s, M2, G2)

    # Total amplitude with interference
    R = r * np.exp(1j * phi)
    A_total = A1 + R * A2

    # Cross section proportional to |A|^2
    sigma = norm * np.abs(A_total)**2

    return sigma


def neg_log_likelihood(params, data_A, data_B, constrained=True):
    """Negative log-likelihood for joint fit."""
    nll = 0

    # Penalty for unphysical parameters
    if constrained:
        norm_A, norm_B, r, phi, M1, G1, M2, G2 = params
        if norm_A < 0 or norm_B < 0 or r < 0 or G1 < 0.01 or G2 < 0.01:
            return 1e10
    else:
        norm_A, r_A, phi_A, norm_B, r_B, phi_B, M1, G1, M2, G2 = params
        if norm_A < 0 or norm_B < 0 or r_A < 0 or r_B < 0 or G1 < 0.01 or G2 < 0.01:
            return 1e10

    # Channel A (J/psi)
    for d in data_A:
        model = cross_section_model(d['sqrt_s'], params, 'A')
        if model <= 0:
            return 1e10
        residual = (d['sigma'] - model) / d['err']
        nll += 0.5 * residual**2

    # Channel B (hc)
    for d in data_B:
        model = cross_section_model(d['sqrt_s'], params, 'B')
        if model <= 0:
            return 1e10
        residual = (d['sigma'] - model) / d['err']
        nll += 0.5 * residual**2

    return nll


def fit_model(data_A, data_B, constrained=True, n_starts=50):
    """Fit two-resonance model with multi-start optimization."""
    best_nll = np.inf
    best_params = None

    for _ in range(n_starts):
        # Random initial parameters
        if constrained:
            # [norm_A, norm_B, r_shared, phi_shared, M1, G1, M2, G2]
            x0 = [
                np.random.uniform(100, 5000),   # norm_A
                np.random.uniform(100, 5000),   # norm_B
                np.random.uniform(0.1, 3.0),    # r_shared
                np.random.uniform(-np.pi, np.pi),  # phi_shared
                Y4220_MASS + np.random.uniform(-0.02, 0.02),  # M1
                Y4220_WIDTH + np.random.uniform(-0.01, 0.01), # G1
                Y4320_MASS + np.random.uniform(-0.02, 0.02),  # M2
                Y4320_WIDTH + np.random.uniform(-0.02, 0.02), # G2
            ]
            bounds = [
                (1, 50000), (1, 50000), (0.01, 10), (-np.pi, np.pi),
                (4.15, 4.28), (0.02, 0.15), (4.28, 4.40), (0.03, 0.20)
            ]
        else:
            # [norm_A, r_A, phi_A, norm_B, r_B, phi_B, M1, G1, M2, G2]
            x0 = [
                np.random.uniform(100, 5000),   # norm_A
                np.random.uniform(0.1, 3.0),    # r_A
                np.random.uniform(-np.pi, np.pi),  # phi_A
                np.random.uniform(100, 5000),   # norm_B
                np.random.uniform(0.1, 3.0),    # r_B
                np.random.uniform(-np.pi, np.pi),  # phi_B
                Y4220_MASS + np.random.uniform(-0.02, 0.02),
                Y4220_WIDTH + np.random.uniform(-0.01, 0.01),
                Y4320_MASS + np.random.uniform(-0.02, 0.02),
                Y4320_WIDTH + np.random.uniform(-0.02, 0.02),
            ]
            bounds = [
                (1, 50000), (0.01, 10), (-np.pi, np.pi),
                (1, 50000), (0.01, 10), (-np.pi, np.pi),
                (4.15, 4.28), (0.02, 0.15), (4.28, 4.40), (0.03, 0.20)
            ]

        try:
            result = minimize(
                neg_log_likelihood,
                x0,
                args=(data_A, data_B, constrained),
                method='L-BFGS-B',
                bounds=bounds,
                options={'maxiter': 2000}
            )

            if result.fun < best_nll:
                best_nll = result.fun
                best_params = result.x

        except Exception:
            continue

    return best_params, best_nll


def compute_chi2_dof(params, data, channel, constrained):
    """Compute chi2/dof for a single channel."""
    chi2_val = 0
    for d in data:
        model = cross_section_model(d['sqrt_s'], params, channel)
        residual = (d['sigma'] - model) / d['err']
        chi2_val += residual**2

    # dof = n_points - n_params_per_channel
    # constrained: 4 params shared + 1 norm per channel = 5
    # unconstrained: 3 params per channel + 4 shared = 7
    n_params = 5 if constrained else 7
    dof = len(data) - n_params
    return chi2_val / max(dof, 1), chi2_val, dof


def bootstrap_trial(args):
    """Single bootstrap trial."""
    data_A, data_B, seed, n_starts = args
    np.random.seed(seed)

    # Resample with replacement
    data_A_boot = [data_A[i] for i in np.random.randint(0, len(data_A), len(data_A))]
    data_B_boot = [data_B[i] for i in np.random.randint(0, len(data_B), len(data_B))]

    # Fit both models
    _, nll_con = fit_model(data_A_boot, data_B_boot, constrained=True, n_starts=n_starts)
    _, nll_unc = fit_model(data_A_boot, data_B_boot, constrained=False, n_starts=n_starts)

    Lambda = 2 * (nll_con - nll_unc)
    return max(Lambda, 0)


def run_rank1_test(data_A, data_B, n_bootstrap=100, n_starts=60):
    """Run complete rank-1 factorization test."""
    print("=" * 60)
    print("BESIII Y(4220)/Y(4320) Rank-1 Factorization Test")
    print("=" * 60)
    print(f"\nChannel A (pi+pi- J/psi): {len(data_A)} points")
    print(f"Channel B (pi+pi- hc): {len(data_B)} points")
    print(f"Energy window: {E_MIN:.2f} - {E_MAX:.2f} GeV")

    # Fit constrained model (shared R)
    print("\n[1/3] Fitting constrained model (shared R)...")
    params_con, nll_con = fit_model(data_A, data_B, constrained=True, n_starts=n_starts)

    # Fit unconstrained model (independent R per channel)
    print("[2/3] Fitting unconstrained model (independent R)...")
    params_unc, nll_unc = fit_model(data_A, data_B, constrained=False, n_starts=n_starts)

    # Likelihood ratio test statistic
    Lambda_obs = 2 * (nll_con - nll_unc)
    Lambda_obs = max(Lambda_obs, 0)

    print(f"\n--- Fit Results ---")
    print(f"NLL constrained:   {nll_con:.2f}")
    print(f"NLL unconstrained: {nll_unc:.2f}")
    print(f"Lambda_obs:        {Lambda_obs:.4f}")

    # Extract coupling ratios
    if params_con is not None:
        r_shared = params_con[2]
        phi_shared = params_con[3]
        print(f"\nShared |R|: {r_shared:.3f}")
        print(f"Shared arg(R): {phi_shared:.2f} rad ({np.degrees(phi_shared):.1f} deg)")

    if params_unc is not None:
        r_A, phi_A = params_unc[1], params_unc[2]
        r_B, phi_B = params_unc[4], params_unc[5]
        print(f"\nChannel A |R|: {r_A:.3f}, arg(R): {np.degrees(phi_A):.1f} deg")
        print(f"Channel B |R|: {r_B:.3f}, arg(R): {np.degrees(phi_B):.1f} deg")

    # Chi2/dof for fit health
    chi2_dof_A, _, _ = compute_chi2_dof(params_con, data_A, 'A', True)
    chi2_dof_B, _, _ = compute_chi2_dof(params_con, data_B, 'B', True)

    print(f"\n--- Fit Health ---")
    print(f"chi2/dof (A): {chi2_dof_A:.3f}")
    print(f"chi2/dof (B): {chi2_dof_B:.3f}")

    # Health gates
    gate_A = 0.3 < chi2_dof_A < 4.0
    gate_B = 0.3 < chi2_dof_B < 4.0
    gates_pass = gate_A and gate_B

    # Bootstrap p-value
    print(f"\n[3/3] Running bootstrap ({n_bootstrap} replicates)...")

    with Pool(max(1, cpu_count() - 1)) as pool:
        args = [(data_A, data_B, seed, max(20, n_starts // 3))
                for seed in range(n_bootstrap)]
        Lambda_boot = list(pool.map(bootstrap_trial, args))

    Lambda_boot = np.array(Lambda_boot)
    p_boot = np.mean(Lambda_boot >= Lambda_obs)
    p_chi2 = 1 - chi2.cdf(Lambda_obs, df=2)

    # Verdict
    if not gates_pass:
        verdict = "INCONCLUSIVE"
    elif p_boot < 0.05:
        verdict = "DISFAVORED"
    else:
        verdict = "NOT_REJECTED"

    print(f"\n" + "=" * 60)
    print("RESULT")
    print("=" * 60)
    print(f"Lambda_obs: {Lambda_obs:.4f}")
    print(f"p_boot:     {p_boot:.4f} ({int(p_boot * n_bootstrap)}/{n_bootstrap} exceedances)")
    print(f"p_chi2(2):  {p_chi2:.4f}")
    print(f"Gates:      {'PASS' if gates_pass else 'FAIL'}")
    print(f"Verdict:    {verdict}")

    return {
        'Lambda_obs': Lambda_obs,
        'p_boot': p_boot,
        'p_chi2': p_chi2,
        'chi2_dof_A': chi2_dof_A,
        'chi2_dof_B': chi2_dof_B,
        'gates_pass': gates_pass,
        'verdict': verdict,
        'r_shared': r_shared if params_con is not None else None,
        'phi_shared': phi_shared if params_con is not None else None,
        'Lambda_boot': Lambda_boot.tolist(),
        'params_con': params_con.tolist() if params_con is not None else None,
        'params_unc': params_unc.tolist() if params_unc is not None else None,
    }


def main():
    parser = argparse.ArgumentParser(description='BESIII Y Rank-1 Test')
    parser.add_argument('--channel-a', default='besiii_y_rank1/extracted/channelA_jpsi_xsec.csv')
    parser.add_argument('--channel-b', default='besiii_y_rank1/extracted/channelB_hc_xsec.csv')
    parser.add_argument('--bootstrap', type=int, default=200)
    parser.add_argument('--starts', type=int, default=80)
    parser.add_argument('--outdir', default='besiii_y_rank1/out')
    args = parser.parse_args()

    # Load data
    data_A = load_channel_data(args.channel_a)
    data_B = load_channel_data(args.channel_b)

    if len(data_A) < 5 or len(data_B) < 5:
        print(f"ERROR: Insufficient data points (A={len(data_A)}, B={len(data_B)})")
        return

    # Run test
    result = run_rank1_test(data_A, data_B, args.bootstrap, args.starts)

    # Save results
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    with open(outdir / 'result.json', 'w') as f:
        json.dump(result, f, indent=2)

    # Write markdown report
    report = f"""# BESIII Y(4220)/Y(4320) Rank-1 Test Result

## Verdict: **{result['verdict']}**

| Metric | Value |
|--------|-------|
| Lambda_obs | {result['Lambda_obs']:.4f} |
| p_boot | {result['p_boot']:.4f} |
| p_chi2(2) | {result['p_chi2']:.4f} |
| chi2/dof (A) | {result['chi2_dof_A']:.3f} |
| chi2/dof (B) | {result['chi2_dof_B']:.3f} |
| Gates | {'PASS' if result['gates_pass'] else 'FAIL'} |

## Shared Coupling Ratio

| Parameter | Value |
|-----------|-------|
| |R| | {result['r_shared']:.3f if result['r_shared'] else 'N/A'} |
| arg(R) | {np.degrees(result['phi_shared']):.1f} deg |

## Data Sources

- Channel A: BESIII e+e- -> pi+pi- J/psi (arXiv:1611.01317)
- Channel B: BESIII e+e- -> pi+pi- hc (HEPData ins2908630)

## Interpretation

{"The rank-1 factorization hypothesis is NOT REJECTED. The Y(4220) and Y(4320) states show a consistent coupling ratio across both decay channels, supporting a common production mechanism." if result['verdict'] == 'NOT_REJECTED' else "The rank-1 factorization hypothesis is DISFAVORED. The coupling ratios differ significantly between channels."}
"""

    with open(outdir / 'REPORT.md', 'w') as f:
        f.write(report)

    print(f"\nResults saved to {outdir}/")


if __name__ == '__main__':
    main()
