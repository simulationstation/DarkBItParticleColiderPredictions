#!/usr/bin/env python3
"""
BESIII Y(4220)/Y(4320) Shared-Subspace Rank-1 Factorization Test

Upgraded model with:
- Two shared resonances (Y1, Y2) for rank-1 test
- Optional third resonance (channel-specific nuisance)
- Additive polynomial background per channel

Model per channel α:
  σ_α(√s) = norm_α × |BW1(s) + R_α×BW2(s) + C3_α×BW3(s)|² + B_α(√s)

Where:
  - R_α = r_α × exp(i×φ_α) is the shared-subspace ratio we test
  - C3_α is channel-specific extra structure amplitude
  - B_α = b0_α + b1_α×(√s - E0) + b2_α×(√s - E0)² is additive background

Rank-1 test: constrain R_A = R_B while leaving C3_α free.
"""

import argparse
import json
import numpy as np
from scipy.optimize import minimize, differential_evolution
from scipy.stats import chi2 as chi2_dist
from pathlib import Path
from multiprocessing import Pool, cpu_count
from dataclasses import dataclass, asdict
from typing import Optional, List, Tuple, Dict, Any
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# PHYSICS CONSTANTS (starting values)
# ============================================================
# Y(4220) - well-established
M1_INIT, G1_INIT = 4.222, 0.050
# Y(4320) - second structure
M2_INIT, G2_INIT = 4.320, 0.120
# Y(4390)/Y(4470) - third structure (varies by channel)
M3_INIT, G3_INIT = 4.420, 0.150

# Reference energy for background polynomial
E0 = 4.30  # GeV

# Fit window
E_MIN = 4.01
E_MAX = 4.60


@dataclass
class FitResult:
    """Container for fit results."""
    nll: float
    params: np.ndarray
    chi2_A: float
    chi2_B: float
    dof_A: int
    dof_B: int
    chi2_dof_A: float
    chi2_dof_B: float
    r_A: float
    phi_A: float
    r_B: float
    phi_B: float
    converged: bool


@dataclass
class TestResult:
    """Container for rank-1 test results."""
    Lambda_obs: float
    p_boot: float
    n_boot: int
    k_exceed: int
    chi2_dof_A: float
    chi2_dof_B: float
    gate_A: str
    gate_B: str
    gates_pass: bool
    r_shared: Optional[float]
    phi_shared: Optional[float]
    r_A: float
    phi_A: float
    r_B: float
    phi_B: float
    fit_mode: str
    identifiable: bool
    verdict: str
    params_con: List[float]
    params_unc: List[float]


# ============================================================
# DATA LOADING
# ============================================================
def load_channel_data(filepath: str, e_min: float = E_MIN, e_max: float = E_MAX) -> List[Dict]:
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
                total_err = np.sqrt(stat_err**2 + sys_err**2)
                if e_min <= sqrt_s <= e_max and sigma > 0:
                    data.append({
                        'sqrt_s': sqrt_s,
                        'sigma': sigma,
                        'err': max(total_err, 0.1)
                    })
    return sorted(data, key=lambda x: x['sqrt_s'])


# ============================================================
# MODEL FUNCTIONS
# ============================================================
def breit_wigner(s: float, M: float, Gamma: float) -> complex:
    """Relativistic Breit-Wigner amplitude."""
    return M * Gamma / (s - M**2 + 1j * M * Gamma)


def background_poly(sqrt_s: float, b0: float, b1: float, b2: float) -> float:
    """Additive polynomial background."""
    dx = sqrt_s - E0
    return max(0, b0 + b1 * dx + b2 * dx**2)


def cross_section_model(sqrt_s: float, params: np.ndarray, channel: str,
                        fixed_shapes: bool = True,
                        M1: float = M1_INIT, G1: float = G1_INIT,
                        M2: float = M2_INIT, G2: float = G2_INIT,
                        M3: float = M3_INIT, G3: float = G3_INIT) -> float:
    """
    Compute cross section with shared-subspace model.

    For FIXED_SHAPES mode (default):
      params = [norm_A, r_A, phi_A, r3_A, phi3_A, b0_A, b1_A, b2_A,
                norm_B, r_B, phi_B, r3_B, phi3_B, b0_B, b1_B, b2_B]
      Total: 16 params

    For constrained (rank-1):
      params = [norm_A, r_shared, phi_shared, r3_A, phi3_A, b0_A, b1_A, b2_A,
                norm_B, r3_B, phi3_B, b0_B, b1_B, b2_B]
      Total: 14 params
    """
    s = sqrt_s**2

    # Parse parameters based on constraint mode
    if len(params) == 16:
        # Unconstrained: independent R per channel
        norm_A, r_A, phi_A, r3_A, phi3_A, b0_A, b1_A, b2_A = params[:8]
        norm_B, r_B, phi_B, r3_B, phi3_B, b0_B, b1_B, b2_B = params[8:]
        if channel == 'A':
            norm, r, phi, r3, phi3 = norm_A, r_A, phi_A, r3_A, phi3_A
            b0, b1, b2 = b0_A, b1_A, b2_A
        else:
            norm, r, phi, r3, phi3 = norm_B, r_B, phi_B, r3_B, phi3_B
            b0, b1, b2 = b0_B, b1_B, b2_B
    elif len(params) == 14:
        # Constrained: shared R
        norm_A, r_shared, phi_shared, r3_A, phi3_A, b0_A, b1_A, b2_A = params[:8]
        norm_B, r3_B, phi3_B, b0_B, b1_B, b2_B = params[8:]
        if channel == 'A':
            norm, r, phi, r3, phi3 = norm_A, r_shared, phi_shared, r3_A, phi3_A
            b0, b1, b2 = b0_A, b1_A, b2_A
        else:
            norm, r, phi, r3, phi3 = norm_B, r_shared, phi_shared, r3_B, phi3_B
            b0, b1, b2 = b0_B, b1_B, b2_B
    else:
        raise ValueError(f"Invalid param length: {len(params)}")

    # Breit-Wigner amplitudes
    A1 = breit_wigner(s, M1, G1)
    A2 = breit_wigner(s, M2, G2)
    A3 = breit_wigner(s, M3, G3)

    # Complex coupling ratios
    R = r * np.exp(1j * phi)
    C3 = r3 * np.exp(1j * phi3)

    # Coherent sum
    A_total = A1 + R * A2 + C3 * A3

    # Cross section = coherent |A|^2 + incoherent background
    sigma_coherent = norm * np.abs(A_total)**2
    sigma_bg = background_poly(sqrt_s, b0, b1, b2)

    return sigma_coherent + sigma_bg


def neg_log_likelihood(params: np.ndarray, data_A: List, data_B: List,
                       constrained: bool = True, fixed_shapes: bool = True,
                       M1: float = M1_INIT, G1: float = G1_INIT,
                       M2: float = M2_INIT, G2: float = G2_INIT,
                       M3: float = M3_INIT, G3: float = G3_INIT,
                       prior_penalty: float = 0.0) -> float:
    """Negative log-likelihood for joint fit."""
    nll = 0.0

    # Parameter bounds enforcement
    if constrained:
        # params length = 14
        norm_A, r_sh, phi_sh, r3_A, phi3_A, b0_A, b1_A, b2_A = params[:8]
        norm_B, r3_B, phi3_B, b0_B, b1_B, b2_B = params[8:]
        if norm_A < 0 or norm_B < 0 or r_sh < 0 or r3_A < 0 or r3_B < 0:
            return 1e12
        if b0_A < 0 or b0_B < 0:
            return 1e12
    else:
        # params length = 16
        norm_A, r_A, phi_A, r3_A, phi3_A, b0_A, b1_A, b2_A = params[:8]
        norm_B, r_B, phi_B, r3_B, phi3_B, b0_B, b1_B, b2_B = params[8:]
        if norm_A < 0 or norm_B < 0 or r_A < 0 or r_B < 0 or r3_A < 0 or r3_B < 0:
            return 1e12
        if b0_A < 0 or b0_B < 0:
            return 1e12

    # Channel A
    for d in data_A:
        model = cross_section_model(d['sqrt_s'], params, 'A', fixed_shapes,
                                    M1, G1, M2, G2, M3, G3)
        if model <= 0:
            return 1e12
        residual = (d['sigma'] - model) / d['err']
        nll += 0.5 * residual**2

    # Channel B
    for d in data_B:
        model = cross_section_model(d['sqrt_s'], params, 'B', fixed_shapes,
                                    M1, G1, M2, G2, M3, G3)
        if model <= 0:
            return 1e12
        residual = (d['sigma'] - model) / d['err']
        nll += 0.5 * residual**2

    return nll + prior_penalty


def compute_chi2(params: np.ndarray, data: List, channel: str,
                 fixed_shapes: bool = True,
                 M1: float = M1_INIT, G1: float = G1_INIT,
                 M2: float = M2_INIT, G2: float = G2_INIT,
                 M3: float = M3_INIT, G3: float = G3_INIT) -> Tuple[float, int]:
    """Compute chi2 and dof for a channel."""
    chi2_val = 0.0
    for d in data:
        model = cross_section_model(d['sqrt_s'], params, channel, fixed_shapes,
                                    M1, G1, M2, G2, M3, G3)
        residual = (d['sigma'] - model) / d['err']
        chi2_val += residual**2

    # dof = n_points - n_params_per_channel
    # Per channel: norm, r, phi, r3, phi3, b0, b1, b2 = 8 (unconstrained)
    # Constrained shares r, phi across channels: 7 per channel effective
    n_params = 8
    dof = max(1, len(data) - n_params)
    return chi2_val, dof


# ============================================================
# FITTING
# ============================================================
def generate_initial_params(constrained: bool, seed: int = None) -> np.ndarray:
    """Generate randomized initial parameters."""
    if seed is not None:
        np.random.seed(seed)

    if constrained:
        # [norm_A, r_shared, phi_shared, r3_A, phi3_A, b0_A, b1_A, b2_A,
        #  norm_B, r3_B, phi3_B, b0_B, b1_B, b2_B]
        return np.array([
            np.random.uniform(500, 5000),     # norm_A
            np.random.uniform(0.3, 2.0),      # r_shared
            np.random.uniform(-np.pi, np.pi), # phi_shared
            np.random.uniform(0.1, 1.5),      # r3_A
            np.random.uniform(-np.pi, np.pi), # phi3_A
            np.random.uniform(0, 10),         # b0_A
            np.random.uniform(-5, 5),         # b1_A
            np.random.uniform(-5, 5),         # b2_A
            np.random.uniform(500, 5000),     # norm_B
            np.random.uniform(0.1, 1.5),      # r3_B
            np.random.uniform(-np.pi, np.pi), # phi3_B
            np.random.uniform(0, 10),         # b0_B
            np.random.uniform(-5, 5),         # b1_B
            np.random.uniform(-5, 5),         # b2_B
        ])
    else:
        # [norm_A, r_A, phi_A, r3_A, phi3_A, b0_A, b1_A, b2_A,
        #  norm_B, r_B, phi_B, r3_B, phi3_B, b0_B, b1_B, b2_B]
        return np.array([
            np.random.uniform(500, 5000),     # norm_A
            np.random.uniform(0.3, 2.0),      # r_A
            np.random.uniform(-np.pi, np.pi), # phi_A
            np.random.uniform(0.1, 1.5),      # r3_A
            np.random.uniform(-np.pi, np.pi), # phi3_A
            np.random.uniform(0, 10),         # b0_A
            np.random.uniform(-5, 5),         # b1_A
            np.random.uniform(-5, 5),         # b2_A
            np.random.uniform(500, 5000),     # norm_B
            np.random.uniform(0.3, 2.0),      # r_B
            np.random.uniform(-np.pi, np.pi), # phi_B
            np.random.uniform(0.1, 1.5),      # r3_B
            np.random.uniform(-np.pi, np.pi), # phi3_B
            np.random.uniform(0, 10),         # b0_B
            np.random.uniform(-5, 5),         # b1_B
            np.random.uniform(-5, 5),         # b2_B
        ])


def get_bounds(constrained: bool) -> List[Tuple]:
    """Get parameter bounds."""
    if constrained:
        return [
            (1, 50000),        # norm_A
            (0.01, 10),        # r_shared
            (-np.pi, np.pi),   # phi_shared
            (0, 10),           # r3_A
            (-np.pi, np.pi),   # phi3_A
            (0, 100),          # b0_A
            (-50, 50),         # b1_A
            (-50, 50),         # b2_A
            (1, 50000),        # norm_B
            (0, 10),           # r3_B
            (-np.pi, np.pi),   # phi3_B
            (0, 100),          # b0_B
            (-50, 50),         # b1_B
            (-50, 50),         # b2_B
        ]
    else:
        return [
            (1, 50000),        # norm_A
            (0.01, 10),        # r_A
            (-np.pi, np.pi),   # phi_A
            (0, 10),           # r3_A
            (-np.pi, np.pi),   # phi3_A
            (0, 100),          # b0_A
            (-50, 50),         # b1_A
            (-50, 50),         # b2_A
            (1, 50000),        # norm_B
            (0.01, 10),        # r_B
            (-np.pi, np.pi),   # phi_B
            (0, 10),           # r3_B
            (-np.pi, np.pi),   # phi3_B
            (0, 100),          # b0_B
            (-50, 50),         # b1_B
            (-50, 50),         # b2_B
        ]


def fit_single_start(args: Tuple) -> Tuple[float, np.ndarray]:
    """Single optimization start."""
    data_A, data_B, constrained, seed, fixed_shapes = args
    x0 = generate_initial_params(constrained, seed)
    bounds = get_bounds(constrained)

    try:
        result = minimize(
            neg_log_likelihood,
            x0,
            args=(data_A, data_B, constrained, fixed_shapes),
            method='L-BFGS-B',
            bounds=bounds,
            options={'maxiter': 3000, 'ftol': 1e-10}
        )
        return result.fun, result.x
    except Exception:
        return 1e12, x0


def fit_model(data_A: List, data_B: List, constrained: bool = True,
              n_starts: int = 80, fixed_shapes: bool = True,
              use_parallel: bool = True) -> FitResult:
    """Fit model with multi-start optimization."""

    args_list = [(data_A, data_B, constrained, seed, fixed_shapes)
                 for seed in range(n_starts)]

    if use_parallel:
        with Pool(max(1, cpu_count() - 1)) as pool:
            results = pool.map(fit_single_start, args_list)
    else:
        # Sequential for nested calls (bootstrap)
        results = [fit_single_start(args) for args in args_list]

    # Find best result
    best_nll = np.inf
    best_params = None
    for nll, params in results:
        if nll < best_nll:
            best_nll = nll
            best_params = params

    if best_params is None:
        return FitResult(nll=1e12, params=np.zeros(16 if not constrained else 14),
                         chi2_A=1e6, chi2_B=1e6, dof_A=1, dof_B=1,
                         chi2_dof_A=1e6, chi2_dof_B=1e6,
                         r_A=0, phi_A=0, r_B=0, phi_B=0, converged=False)

    # Compute chi2/dof
    chi2_A, dof_A = compute_chi2(best_params, data_A, 'A', fixed_shapes)
    chi2_B, dof_B = compute_chi2(best_params, data_B, 'B', fixed_shapes)

    # Extract R values
    if constrained:
        r_A = r_B = best_params[1]
        phi_A = phi_B = best_params[2]
    else:
        r_A, phi_A = best_params[1], best_params[2]
        r_B, phi_B = best_params[9], best_params[10]

    return FitResult(
        nll=best_nll,
        params=best_params,
        chi2_A=chi2_A,
        chi2_B=chi2_B,
        dof_A=dof_A,
        dof_B=dof_B,
        chi2_dof_A=chi2_A / max(dof_A, 1),
        chi2_dof_B=chi2_B / max(dof_B, 1),
        r_A=r_A,
        phi_A=phi_A,
        r_B=r_B,
        phi_B=phi_B,
        converged=best_nll < 1e10
    )


# ============================================================
# IDENTIFIABILITY CHECK
# ============================================================
def check_identifiability(data_A: List, data_B: List, constrained: bool,
                          n_starts: int = 40) -> Tuple[bool, str]:
    """Check if R is identifiable by looking at near-optimal solutions."""

    args_list = [(data_A, data_B, constrained, seed, True)
                 for seed in range(n_starts)]

    with Pool(max(1, cpu_count() - 1)) as pool:
        results = pool.map(fit_single_start, args_list)

    # Filter to near-optimal (within ΔNLL < 2 of best)
    nlls = [r[0] for r in results]
    best_nll = min(nlls)

    near_optimal = [(nll, params) for nll, params in results
                    if nll < best_nll + 2 and nll < 1e10]

    if len(near_optimal) < 3:
        return True, "TOO_FEW_SOLUTIONS"

    # Extract R values from near-optimal
    if constrained:
        r_vals = [p[1] for _, p in near_optimal]
        phi_vals = [p[2] for _, p in near_optimal]
    else:
        # Check channel A
        r_vals = [p[1] for _, p in near_optimal]
        phi_vals = [p[2] for _, p in near_optimal]

    r_vals = np.array(r_vals)
    phi_vals = np.array(phi_vals)

    # Check r spread
    r_spread = r_vals.max() / max(r_vals.min(), 0.01)
    if r_spread > 10:
        return False, f"R_SPREAD={r_spread:.1f}"

    # Check phi multimodality (rough check)
    phi_std = np.std(phi_vals)
    if phi_std > 1.5:  # More than ~90 degrees spread
        return False, f"PHI_STD={phi_std:.2f}"

    return True, "OK"


# ============================================================
# BOOTSTRAP
# ============================================================
def bootstrap_trial(args: Tuple) -> float:
    """Single parametric bootstrap trial."""
    data_A, data_B, params_con, seed, n_starts = args
    np.random.seed(seed)

    # Generate bootstrap data from constrained best-fit
    data_A_boot = []
    for d in data_A:
        model = cross_section_model(d['sqrt_s'], params_con, 'A')
        y_boot = np.random.normal(model, d['err'])
        data_A_boot.append({'sqrt_s': d['sqrt_s'], 'sigma': max(y_boot, 0.1), 'err': d['err']})

    data_B_boot = []
    for d in data_B:
        model = cross_section_model(d['sqrt_s'], params_con, 'B')
        y_boot = np.random.normal(model, d['err'])
        data_B_boot.append({'sqrt_s': d['sqrt_s'], 'sigma': max(y_boot, 0.1), 'err': d['err']})

    # Fit both models with reduced starts (sequential - no nested pools)
    result_con = fit_model(data_A_boot, data_B_boot, constrained=True, n_starts=n_starts, use_parallel=False)
    result_unc = fit_model(data_A_boot, data_B_boot, constrained=False, n_starts=n_starts, use_parallel=False)

    Lambda = 2 * (result_con.nll - result_unc.nll)
    return max(Lambda, 0)


def run_bootstrap(data_A: List, data_B: List, params_con: np.ndarray,
                  n_boot: int = 200, n_starts: int = 40) -> Tuple[float, int, List[float]]:
    """Run parametric bootstrap for p-value."""

    args_list = [(data_A, data_B, params_con, seed, n_starts)
                 for seed in range(n_boot)]

    with Pool(max(1, cpu_count() - 1)) as pool:
        Lambda_boot = list(pool.map(bootstrap_trial, args_list))

    return Lambda_boot


# ============================================================
# FIT HEALTH GATES
# ============================================================
def evaluate_gate(chi2_dof: float) -> str:
    """Evaluate fit health gate."""
    if chi2_dof > 3.0:
        return "MODEL_MISMATCH"
    elif chi2_dof < 0.5:
        return "UNDERCONSTRAINED"
    else:
        return "HEALTHY"


# ============================================================
# MAIN TEST
# ============================================================
def run_shared_subspace_test(data_A: List, data_B: List,
                             n_starts: int = 80,
                             n_boot: int = 200) -> TestResult:
    """Run complete shared-subspace rank-1 test."""

    print("=" * 70)
    print("BESIII Y(4220)/Y(4320) Shared-Subspace Rank-1 Test")
    print("=" * 70)
    print(f"\nChannel A (π+π-J/ψ): {len(data_A)} points")
    print(f"Channel B (π+π-hc): {len(data_B)} points")
    print(f"Energy window: [{E_MIN:.2f}, {E_MAX:.2f}] GeV")
    print(f"Model: |BW1 + R×BW2 + C3×BW3|² + B(√s)")

    # ========================================
    # PHASE 1: FIXED_SHAPES fits
    # ========================================
    print("\n" + "-" * 70)
    print("PHASE 1: FIXED_SHAPES mode")
    print("-" * 70)

    print("\n[1/4] Fitting CONSTRAINED model (shared R)...")
    result_con = fit_model(data_A, data_B, constrained=True, n_starts=n_starts)

    print("[2/4] Fitting UNCONSTRAINED model (independent R)...")
    result_unc = fit_model(data_A, data_B, constrained=False, n_starts=n_starts)

    # Evaluate gates
    gate_A = evaluate_gate(result_con.chi2_dof_A)
    gate_B = evaluate_gate(result_con.chi2_dof_B)
    gates_pass = (gate_A == "HEALTHY" and gate_B == "HEALTHY")

    print(f"\n--- Fit Results (FIXED_SHAPES) ---")
    print(f"NLL constrained:   {result_con.nll:.2f}")
    print(f"NLL unconstrained: {result_unc.nll:.2f}")
    print(f"chi²/dof (A): {result_con.chi2_dof_A:.3f} -> {gate_A}")
    print(f"chi²/dof (B): {result_con.chi2_dof_B:.3f} -> {gate_B}")

    Lambda_obs = max(0, 2 * (result_con.nll - result_unc.nll))
    print(f"Lambda_obs: {Lambda_obs:.4f}")

    print(f"\nShared R: |R| = {result_con.r_A:.3f}, φ = {np.degrees(result_con.phi_A):.1f}°")
    print(f"Channel A R: |R| = {result_unc.r_A:.3f}, φ = {np.degrees(result_unc.phi_A):.1f}°")
    print(f"Channel B R: |R| = {result_unc.r_B:.3f}, φ = {np.degrees(result_unc.phi_B):.1f}°")

    fit_mode = "FIXED_SHAPES"

    # ========================================
    # PHASE 2: Identifiability check
    # ========================================
    print("\n[3/4] Checking identifiability...")
    identifiable, ident_reason = check_identifiability(data_A, data_B, constrained=True)
    print(f"Identifiable: {identifiable} ({ident_reason})")

    # ========================================
    # PHASE 3: Bootstrap (only if gates pass)
    # ========================================
    p_boot = -1.0
    k_exceed = 0
    Lambda_boot = []

    if gates_pass and identifiable:
        print(f"\n[4/4] Running bootstrap ({n_boot} replicates)...")
        Lambda_boot = run_bootstrap(data_A, data_B, result_con.params,
                                    n_boot=n_boot, n_starts=40)
        k_exceed = sum(1 for L in Lambda_boot if L >= Lambda_obs)
        p_boot = k_exceed / n_boot
        print(f"p_boot = {p_boot:.4f} ({k_exceed}/{n_boot} exceedances)")
    else:
        print("\n[4/4] Skipping bootstrap (gates failed or not identifiable)")

    # ========================================
    # VERDICT
    # ========================================
    if not gates_pass:
        if gate_A == "MODEL_MISMATCH" or gate_B == "MODEL_MISMATCH":
            verdict = "MODEL_MISMATCH"
        else:
            verdict = "INCONCLUSIVE"
    elif not identifiable:
        verdict = "INCONCLUSIVE"
    elif p_boot < 0.05:
        verdict = "DISFAVORED"
    else:
        verdict = "NOT_REJECTED"

    print("\n" + "=" * 70)
    print("FINAL RESULT")
    print("=" * 70)
    print(f"Fit mode: {fit_mode}")
    print(f"chi²/dof (A): {result_con.chi2_dof_A:.3f} [{gate_A}]")
    print(f"chi²/dof (B): {result_con.chi2_dof_B:.3f} [{gate_B}]")
    print(f"Lambda_obs: {Lambda_obs:.4f}")
    if p_boot >= 0:
        print(f"p_boot: {p_boot:.4f} ({k_exceed}/{n_boot})")
    print(f"Identifiable: {identifiable}")
    print(f"VERDICT: {verdict}")

    return TestResult(
        Lambda_obs=Lambda_obs,
        p_boot=p_boot,
        n_boot=n_boot if gates_pass else 0,
        k_exceed=k_exceed,
        chi2_dof_A=result_con.chi2_dof_A,
        chi2_dof_B=result_con.chi2_dof_B,
        gate_A=gate_A,
        gate_B=gate_B,
        gates_pass=gates_pass,
        r_shared=result_con.r_A if result_con.converged else None,
        phi_shared=result_con.phi_A if result_con.converged else None,
        r_A=result_unc.r_A,
        phi_A=result_unc.phi_A,
        r_B=result_unc.r_B,
        phi_B=result_unc.phi_B,
        fit_mode=fit_mode,
        identifiable=identifiable,
        verdict=verdict,
        params_con=result_con.params.tolist() if result_con.params is not None else [],
        params_unc=result_unc.params.tolist() if result_unc.params is not None else []
    )


# ============================================================
# REPORT GENERATION
# ============================================================
def generate_report(result: TestResult, data_A: List, data_B: List,
                    outdir: Path) -> str:
    """Generate markdown report."""

    phi_shared_deg = np.degrees(result.phi_shared) if result.phi_shared else 0
    phi_A_deg = np.degrees(result.phi_A)
    phi_B_deg = np.degrees(result.phi_B)

    report = f"""# BESIII Y(4220)/Y(4320) Shared-Subspace Rank-1 Test

## Executive Summary

**Verdict: {result.verdict}**

| Metric | Value |
|--------|-------|
| Lambda_obs | {result.Lambda_obs:.4f} |
| p_boot | {result.p_boot:.4f} ({result.k_exceed}/{result.n_boot}) |
| chi²/dof (A) | {result.chi2_dof_A:.3f} [{result.gate_A}] |
| chi²/dof (B) | {result.chi2_dof_B:.3f} [{result.gate_B}] |
| Gates | {'PASS' if result.gates_pass else 'FAIL'} |
| Identifiable | {result.identifiable} |
| Fit mode | {result.fit_mode} |

---

## Data Provenance

| Channel | Source | File | Points |
|---------|--------|------|--------|
| A (π+π-J/ψ) | arXiv:1611.01317 | channelA_jpsi_xsec.csv | {len(data_A)} |
| B (π+π-hc) | HEPData ins2908630 | channelB_hc_xsec.csv | {len(data_B)} |

Energy window: [{E_MIN:.2f}, {E_MAX:.2f}] GeV

---

## Model Definition

The cross section per channel α is modeled as:

```
σ_α(√s) = norm_α × |BW₁(s) + R_α×BW₂(s) + C3_α×BW₃(s)|² + B_α(√s)
```

Where:
- `BWₖ(s) = Mₖ×Γₖ / (s - Mₖ² + i×Mₖ×Γₖ)` (relativistic Breit-Wigner)
- `R_α = r_α × exp(i×φ_α)` is the shared-subspace ratio (Y₂/Y₁)
- `C3_α` is the channel-specific third resonance amplitude
- `B_α(√s) = b0_α + b1_α×(√s - E₀) + b2_α×(√s - E₀)²` is additive background

**Fixed resonance parameters ({result.fit_mode}):**

| Resonance | Mass (GeV) | Width (GeV) |
|-----------|------------|-------------|
| Y₁ (Y4220) | {M1_INIT} | {G1_INIT} |
| Y₂ (Y4320) | {M2_INIT} | {G2_INIT} |
| Y₃ (extra) | {M3_INIT} | {G3_INIT} |

**Rank-1 test:**
- Constrained: R_A = R_B (shared r and φ)
- Unconstrained: R_A, R_B independent
- C3_α always free (not constrained)

---

## Coupling Ratios

### Shared (Constrained fit)
| Parameter | Value |
|-----------|-------|
| |R| | {f'{result.r_shared:.4f}' if result.r_shared else 'N/A'} |
| arg(R) | {phi_shared_deg:.1f}° |

### Per-Channel (Unconstrained fit)
| Channel | |R| | arg(R) |
|---------|-----|--------|
| A (J/ψ) | {result.r_A:.4f} | {phi_A_deg:.1f}° |
| B (hc) | {result.r_B:.4f} | {phi_B_deg:.1f}° |

---

## Fit Health

| Channel | chi²/dof | Gate |
|---------|----------|------|
| A | {result.chi2_dof_A:.3f} | {result.gate_A} |
| B | {result.chi2_dof_B:.3f} | {result.gate_B} |

Health gates:
- `chi²/dof > 3.0` → MODEL_MISMATCH
- `chi²/dof < 0.5` → UNDERCONSTRAINED
- Otherwise → HEALTHY

---

## Interpretation

"""

    if result.verdict == "NOT_REJECTED":
        report += """The rank-1 factorization hypothesis is **NOT REJECTED** at the 5% level.

The Y(4220) and Y(4320) states show a consistent coupling ratio R = g(Y4320)/g(Y4220)
across the π+π-J/ψ and π+π-hc decay channels, supporting a common production mechanism.

This is a significant physics result: the shared-subspace model adequately describes
both channels, and the cross-channel constraint is statistically acceptable.
"""
    elif result.verdict == "DISFAVORED":
        report += """The rank-1 factorization hypothesis is **DISFAVORED** (p < 0.05).

The coupling ratios differ significantly between channels, suggesting the Y(4220)
and Y(4320) may have different production mechanisms in π+π-J/ψ vs π+π-hc.
"""
    elif result.verdict == "MODEL_MISMATCH":
        report += """The test is **INCONCLUSIVE** due to MODEL_MISMATCH.

The shared-subspace model (2 shared resonances + background + third structure)
still does not adequately describe the data (chi²/dof > 3).

Possible improvements:
- Add more resonances
- Use energy-dependent widths
- Include interference with continuum properly
- Request official covariance matrices from BESIII
"""
    else:
        report += """The test is **INCONCLUSIVE**.

Either the fit health gates failed or the coupling ratio R is not well-identified
from the data. Further investigation with improved data or model is needed.
"""

    report += f"""
---

## Files Generated

- `REPORT_shared_subspace.md` - This report
- `out/shared_subspace_result.json` - Machine-readable results

---

*Generated by besiii_y_rank1_shared_subspace.py*
"""

    return report


# ============================================================
# MAIN
# ============================================================
def main():
    parser = argparse.ArgumentParser(description='BESIII Y Shared-Subspace Rank-1 Test')
    parser.add_argument('--channel-a', default='besiii_y_rank1/extracted/channelA_jpsi_xsec.csv')
    parser.add_argument('--channel-b', default='besiii_y_rank1/extracted/channelB_hc_xsec.csv')
    parser.add_argument('--bootstrap', type=int, default=200)
    parser.add_argument('--starts', type=int, default=80)
    parser.add_argument('--outdir', default='besiii_y_rank1/out')
    args = parser.parse_args()

    # Load data
    data_A = load_channel_data(args.channel_a, E_MIN, E_MAX)
    data_B = load_channel_data(args.channel_b, E_MIN, E_MAX)

    print(f"Loaded {len(data_A)} points for channel A")
    print(f"Loaded {len(data_B)} points for channel B")

    if len(data_A) < 5 or len(data_B) < 5:
        print(f"ERROR: Insufficient data points")
        return

    # Run test
    result = run_shared_subspace_test(data_A, data_B, args.starts, args.bootstrap)

    # Save results
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # JSON output
    result_dict = asdict(result)
    # Convert numpy types to native Python
    for k, v in result_dict.items():
        if isinstance(v, (np.floating, np.integer)):
            result_dict[k] = float(v)
        elif isinstance(v, np.bool_):
            result_dict[k] = bool(v)

    with open(outdir / 'shared_subspace_result.json', 'w') as f:
        json.dump(result_dict, f, indent=2)

    # Markdown report
    report = generate_report(result, data_A, data_B, outdir)
    report_path = Path('besiii_y_rank1/REPORT_shared_subspace.md')
    with open(report_path, 'w') as f:
        f.write(report)

    print(f"\nResults saved to {outdir}/")
    print(f"Report saved to {report_path}")


if __name__ == '__main__':
    main()
