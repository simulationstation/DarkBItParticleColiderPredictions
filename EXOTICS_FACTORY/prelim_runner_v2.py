#!/usr/bin/env python3
"""
prelim_runner_v2.py - Run preliminary rank-1 fits for all EXOTICS_FACTORY families

Features:
- Uses validated sim_generate + sim_fit_v3 pipeline
- Per-family BW parameter configurations
- Incremental results (writes after each family completes)
- Per-family timeout (8 min max)
- Resumable (skips already-completed families)

Exclusions: Zc-like, CMS X(6900)/X(7100), di-charmonium (already done)
"""

import sys
import os
import json
import csv
import time
import signal
from datetime import datetime
from pathlib import Path
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Dict, Optional, List

# Add FIX_THESE to path for sim modules
sys.path.insert(0, str(Path(__file__).parent.parent / 'FIX_THESE'))

import numpy as np
from sim_generate import generate_dataset, BWParams
from sim_fit_v3 import (
    run_calibration_trial, fit_joint_unconstrained, fit_joint_constrained,
    compute_lambda, check_fit_health, bootstrap_pvalue_calibrated,
    CHI2_DOF_MIN, CHI2_DOF_MAX, DEVIANCE_DOF_MAX, DEFAULT_BOOTSTRAP
)

# ============================================================
# FAMILY CONFIGURATIONS
# ============================================================

FAMILY_CONFIGS = {
    "belle_zb": {
        "name": "Belle Zb(10610)/Zb(10650)",
        "category": "exotic",
        "states": ["Zb(10610)", "Zb(10650)"],
        "channels": ["π Υ(1S)", "π Υ(2S)"],
        "bw_params": {"m1": 0.0, "m2": 1.0, "gamma1": 0.41, "gamma2": 0.26},
        "R_true": {"r": 0.7, "phi_deg": -60},
        "deltaR_M1": {"dr": 0.3, "dphi_deg": 45},
        "deltaR_M4": {"dr": 0.15, "dphi_deg": 20},
        "channelA": {"type": "poisson", "nbins": 40, "x_range": [-0.5, 1.5], "mean_count_scale": 80, "bg_level": 0.25},
        "channelB": {"type": "poisson", "nbins": 40, "x_range": [-0.5, 1.5], "mean_count_scale": 60, "bg_level": 0.30},
        "expected_verdict": "SUPPORTED",
    },

    "lhcb_pc_doublet": {
        "name": "LHCb Pc(4440)/Pc(4457)",
        "category": "exotic",
        "states": ["Pc(4440)", "Pc(4457)"],
        "channels": ["J/ψ p (full)", "J/ψ p (tight)"],
        "bw_params": {"m1": 0.0, "m2": 1.0, "gamma1": 1.21, "gamma2": 0.38},
        "R_true": {"r": 0.5, "phi_deg": -45},
        "deltaR_M1": {"dr": 0.25, "dphi_deg": 60},
        "deltaR_M4": {"dr": 0.12, "dphi_deg": 25},
        "channelA": {"type": "poisson", "nbins": 35, "x_range": [-0.5, 1.5], "mean_count_scale": 100, "bg_level": 0.40},
        "channelB": {"type": "poisson", "nbins": 35, "x_range": [-0.5, 1.5], "mean_count_scale": 70, "bg_level": 0.35},
        "expected_verdict": "SUPPORTED",
    },

    "strange_pcs": {
        "name": "Strange Pcs(4459)/Pcs(4338)",
        "category": "exotic",
        "states": ["Pcs(4459)", "Pcs(4338)"],
        "channels": ["J/ψ Λ (primary)", "J/ψ Λ (alt)"],
        "bw_params": {"m1": 0.0, "m2": 1.0, "gamma1": 0.06, "gamma2": 0.14},
        "R_true": {"r": 0.6, "phi_deg": -30},
        "deltaR_M1": {"dr": 0.2, "dphi_deg": 50},
        "deltaR_M4": {"dr": 0.1, "dphi_deg": 20},
        "channelA": {"type": "poisson", "nbins": 30, "x_range": [-0.5, 1.5], "mean_count_scale": 50, "bg_level": 0.35},
        "channelB": {"type": "poisson", "nbins": 30, "x_range": [-0.5, 1.5], "mean_count_scale": 40, "bg_level": 0.40},
        "expected_verdict": "SUPPORTED",
    },

    "besiii_y_pipijpsi_hc": {
        "name": "BESIII Y(4220)/Y(4320)",
        "category": "exotic",
        "states": ["Y(4220)", "Y(4320)"],
        "channels": ["π+π- J/ψ", "π+π- hc"],
        "bw_params": {"m1": 0.0, "m2": 1.0, "gamma1": 0.45, "gamma2": 1.03},
        "R_true": {"r": 0.8, "phi_deg": -75},
        "deltaR_M1": {"dr": 0.35, "dphi_deg": 55},
        "deltaR_M4": {"dr": 0.18, "dphi_deg": 25},
        "channelA": {"type": "poisson", "nbins": 45, "x_range": [-0.5, 1.5], "mean_count_scale": 90, "bg_level": 0.20},
        "channelB": {"type": "poisson", "nbins": 40, "x_range": [-0.5, 1.5], "mean_count_scale": 55, "bg_level": 0.30},
        "expected_verdict": "SUPPORTED",
    },

    "besiii_belle_isr_y": {
        "name": "BESIII/Belle ISR Y(4260)/Y(4360)",
        "category": "exotic",
        "states": ["Y(4260)", "Y(4360)"],
        "channels": ["ISR π+π- J/ψ", "ISR π+π- ψ(2S)"],
        "bw_params": {"m1": 0.0, "m2": 1.0, "gamma1": 0.40, "gamma2": 0.70},
        "R_true": {"r": 0.65, "phi_deg": -55},
        "deltaR_M1": {"dr": 0.30, "dphi_deg": 40},
        "deltaR_M4": {"dr": 0.15, "dphi_deg": 18},
        "channelA": {"type": "poisson", "nbins": 40, "x_range": [-0.5, 1.5], "mean_count_scale": 70, "bg_level": 0.25},
        "channelB": {"type": "poisson", "nbins": 35, "x_range": [-0.5, 1.5], "mean_count_scale": 45, "bg_level": 0.35},
        "expected_verdict": "SUPPORTED",
    },

    "control_babar_omega": {
        "name": "BaBar ω(1420)/ω(1650) [CONTROL]",
        "category": "control",
        "states": ["ω(1420)", "ω(1650)"],
        "channels": ["ω π+π-", "ω f0"],
        "bw_params": {"m1": 0.0, "m2": 1.0, "gamma1": 1.16, "gamma2": 1.26},
        "R_true": {"r": 1.2, "phi_deg": -150},
        "deltaR_M1": {"dr": 0.8, "dphi_deg": 90},
        "deltaR_M4": {"dr": 0.4, "dphi_deg": 45},
        "channelA": {"type": "poisson", "nbins": 40, "x_range": [-0.5, 1.5], "mean_count_scale": 85, "bg_level": 0.20},
        "channelB": {"type": "poisson", "nbins": 40, "x_range": [-0.5, 1.5], "mean_count_scale": 65, "bg_level": 0.25},
        "expected_verdict": "DISFAVORED",
        "use_M1": True,
    },

    "control_babar_phi": {
        "name": "BaBar φ(1680)/φ(2170) [CONTROL]",
        "category": "control",
        "states": ["φ(1680)", "φ(2170)"],
        "channels": ["φ f0", "φ π+π-"],
        "bw_params": {"m1": 0.0, "m2": 1.0, "gamma1": 0.31, "gamma2": 0.17},
        "R_true": {"r": 0.9, "phi_deg": -120},
        "deltaR_M1": {"dr": 0.7, "dphi_deg": 85},
        "deltaR_M4": {"dr": 0.35, "dphi_deg": 40},
        "channelA": {"type": "poisson", "nbins": 40, "x_range": [-0.5, 1.5], "mean_count_scale": 75, "bg_level": 0.22},
        "channelB": {"type": "poisson", "nbins": 40, "x_range": [-0.5, 1.5], "mean_count_scale": 60, "bg_level": 0.28},
        "expected_verdict": "DISFAVORED",
        "use_M1": True,
    },

    "x3872_like": {
        "name": "X(3872)-like single pole",
        "category": "exotic",
        "states": ["X(3872)"],
        "channels": ["J/ψ π+π-", "J/ψ π+π-π0"],
        "bw_params": {"m1": 0.0, "m2": 1.0, "gamma1": 0.001, "gamma2": 0.05},
        "R_true": {"r": 0.1, "phi_deg": 0},
        "deltaR_M1": {"dr": 0.05, "dphi_deg": 30},
        "deltaR_M4": {"dr": 0.02, "dphi_deg": 15},
        "channelA": {"type": "poisson", "nbins": 50, "x_range": [-0.5, 1.5], "mean_count_scale": 120, "bg_level": 0.15},
        "channelB": {"type": "poisson", "nbins": 45, "x_range": [-0.5, 1.5], "mean_count_scale": 80, "bg_level": 0.20},
        "expected_verdict": "SUPPORTED",
    },
}

# ============================================================
# TIMEOUT HANDLER
# ============================================================

class TimeoutError(Exception):
    pass

@contextmanager
def timeout(seconds):
    def signal_handler(signum, frame):
        raise TimeoutError(f"Timed out after {seconds} seconds")
    old_handler = signal.signal(signal.SIGALRM, signal_handler)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old_handler)

# ============================================================
# PRELIM FIT RUNNER
# ============================================================

@dataclass
class PrelimResult:
    family_id: str
    category: str
    status: str
    gates: str
    chi2_dof_A: float
    chi2_dof_B: float
    r_shared: float
    phi_shared: float
    Lambda_obs: float
    p_boot: float
    n_boot: int
    verdict: str
    expected: str
    runtime_sec: float
    error_msg: str = ""


def run_prelim_fit(family_id: str, config: Dict, n_bootstrap: int = 100,
                   n_starts: int = 120, seed: int = 42) -> PrelimResult:
    """Run preliminary fit for a single family using run_calibration_trial."""

    start_time = time.time()

    try:
        # Determine mechanism: M0 for exotics, M1 for controls
        mechanism = "M1" if config.get("use_M1", False) else "M0"

        # Generate synthetic data
        dataset = generate_dataset(config, mechanism, scale_factor=1.0, seed=seed)

        # Run calibration trial (this does everything)
        result = run_calibration_trial(dataset, n_bootstrap=n_bootstrap, n_starts=n_starts)

        runtime = time.time() - start_time

        # Extract results
        if not result.get('converged', False):
            return PrelimResult(
                family_id=family_id,
                category=config["category"],
                status="FIT_FAILED",
                gates="FIT_FAILED",
                chi2_dof_A=0.0,
                chi2_dof_B=0.0,
                r_shared=0.0,
                phi_shared=0.0,
                Lambda_obs=0.0,
                p_boot=0.0,
                n_boot=0,
                verdict="FIT_FAILED",
                expected=config.get("expected_verdict", "UNKNOWN"),
                runtime_sec=runtime
            )

        gates = result.get('gates', 'UNKNOWN')
        chi2_dof_a = result.get('chi2_dof_a', 0)
        chi2_dof_b = result.get('chi2_dof_b', 0)
        lambda_obs = result.get('lambda_obs', 0)
        p_boot = result.get('p_boot', 0)
        r_shared = result.get('r_con', 0)
        phi_shared = result.get('phi_con', 0)

        # Determine verdict
        if gates == 'MISMATCH':
            verdict = "MODEL_MISMATCH"
        elif gates == 'UNDERCONSTRAINED':
            verdict = "INCONCLUSIVE"
        elif p_boot >= 0.05:
            verdict = "NOT_REJECTED"
        else:
            verdict = "DISFAVORED"

        return PrelimResult(
            family_id=family_id,
            category=config["category"],
            status="COMPLETED",
            gates=gates,
            chi2_dof_A=chi2_dof_a,
            chi2_dof_B=chi2_dof_b,
            r_shared=r_shared,
            phi_shared=phi_shared,
            Lambda_obs=lambda_obs,
            p_boot=p_boot if not np.isnan(p_boot) else 0.0,
            n_boot=n_bootstrap if gates == 'PASS' else 0,
            verdict=verdict,
            expected=config.get("expected_verdict", "UNKNOWN"),
            runtime_sec=runtime
        )

    except Exception as e:
        runtime = time.time() - start_time
        return PrelimResult(
            family_id=family_id,
            category=config["category"],
            status="ERROR",
            gates="N/A",
            chi2_dof_A=0.0,
            chi2_dof_B=0.0,
            r_shared=0.0,
            phi_shared=0.0,
            Lambda_obs=0.0,
            p_boot=0.0,
            n_boot=0,
            verdict="ERROR",
            expected=config.get("expected_verdict", "UNKNOWN"),
            runtime_sec=runtime,
            error_msg=str(e)
        )


def write_result_to_csv(result: PrelimResult, csv_path: Path, write_header: bool = False):
    """Append a single result to CSV (incremental)."""
    fieldnames = [
        "family_id", "category", "status", "gates",
        "chi2_dof_A", "chi2_dof_B",
        "r_shared", "phi_shared", "Lambda_obs", "p_boot", "n_boot",
        "verdict", "expected", "runtime_sec", "error_msg"
    ]

    mode = 'w' if write_header else 'a'
    with open(csv_path, mode, newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if write_header:
            writer.writeheader()
        writer.writerow({
            "family_id": result.family_id,
            "category": result.category,
            "status": result.status,
            "gates": result.gates,
            "chi2_dof_A": f"{result.chi2_dof_A:.3f}",
            "chi2_dof_B": f"{result.chi2_dof_B:.3f}",
            "r_shared": f"{result.r_shared:.3f}",
            "phi_shared": f"{result.phi_shared:.1f}",
            "Lambda_obs": f"{result.Lambda_obs:.3f}",
            "p_boot": f"{result.p_boot:.4f}",
            "n_boot": result.n_boot,
            "verdict": result.verdict,
            "expected": result.expected,
            "runtime_sec": f"{result.runtime_sec:.1f}",
            "error_msg": result.error_msg
        })


def load_completed_families(csv_path: Path) -> set:
    """Load already-completed families for resumability."""
    if not csv_path.exists():
        return set()
    completed = set()
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            completed.add(row['family_id'])
    return completed


def generate_summary_md(results: List[PrelimResult], output_dir: Path):
    """Generate markdown summary."""
    md_path = output_dir / "PRELIM_SUMMARY_v2.md"

    with open(md_path, 'w') as f:
        f.write("# EXOTICS_FACTORY Preliminary Fit Results (v2)\n\n")
        f.write(f"Generated: {datetime.now().isoformat()}\n\n")

        f.write("## Results Summary\n\n")
        f.write("| Family | Category | Verdict | p_boot | Λ | |R|_shared | Expected | Match |\n")
        f.write("|--------|----------|---------|--------|---|----------|----------|-------|\n")

        for r in results:
            match = "✓" if r.verdict == r.expected or (r.verdict == "NOT_REJECTED" and r.expected == "SUPPORTED") else "✗"
            f.write(f"| {r.family_id} | {r.category} | {r.verdict} | {r.p_boot:.3f} | {r.Lambda_obs:.2f} | {r.r_shared:.3f} | {r.expected} | {match} |\n")

        f.write("\n## Families Supporting Rank-1 Factorization\n\n")
        supported = [r for r in results if r.verdict == "NOT_REJECTED"]
        if supported:
            for r in supported:
                f.write(f"- **{r.family_id}**: p={r.p_boot:.3f}, |R|={r.r_shared:.3f}, φ={r.phi_shared:.1f}°\n")
        else:
            f.write("*None*\n")

        f.write("\n## Families Disfavoring Rank-1 (Controls)\n\n")
        disfavored = [r for r in results if r.verdict == "DISFAVORED"]
        if disfavored:
            for r in disfavored:
                f.write(f"- **{r.family_id}**: p={r.p_boot:.3f}\n")
        else:
            f.write("*None*\n")

        f.write("\n## Inconclusive / Errors\n\n")
        other = [r for r in results if r.verdict in ("INCONCLUSIVE", "MODEL_MISMATCH", "TIMEOUT", "ERROR", "FIT_FAILED")]
        if other:
            for r in other:
                f.write(f"- **{r.family_id}**: {r.verdict}")
                if r.error_msg:
                    f.write(f" ({r.error_msg})")
                f.write("\n")
        else:
            f.write("*None*\n")

        f.write("\n---\n*Generated by prelim_runner_v2.py*\n")


def main():
    """Main entry point."""
    print("=" * 60)
    print("EXOTICS_FACTORY Preliminary Fits v2")
    print("=" * 60)

    output_dir = Path("EXOTICS_FACTORY/out")
    output_dir.mkdir(parents=True, exist_ok=True)

    csv_path = output_dir / "PRELIM_MASTER_TABLE_v2.csv"

    completed = load_completed_families(csv_path)
    if completed:
        print(f"\nResuming: {len(completed)} families already completed")

    families_to_run = [fid for fid in FAMILY_CONFIGS.keys() if fid not in completed]

    if not families_to_run:
        print("\nAll families already completed!")
        return

    print(f"\nRunning {len(families_to_run)} families: {families_to_run}")
    print(f"Settings: n_bootstrap=100, n_starts=120, timeout=480s\n")

    write_header = len(completed) == 0
    results = []

    for i, family_id in enumerate(families_to_run):
        config = FAMILY_CONFIGS[family_id]
        print(f"\n[{i+1}/{len(families_to_run)}] Running {family_id} ({config['category']})...")
        sys.stdout.flush()

        try:
            with timeout(480):
                result = run_prelim_fit(
                    family_id=family_id,
                    config=config,
                    n_bootstrap=100,
                    n_starts=120,
                    seed=42 + i
                )
        except TimeoutError:
            result = PrelimResult(
                family_id=family_id,
                category=config["category"],
                status="TIMEOUT",
                gates="N/A",
                chi2_dof_A=0, chi2_dof_B=0,
                r_shared=0, phi_shared=0,
                Lambda_obs=0, p_boot=0, n_boot=0,
                verdict="TIMEOUT",
                expected=config.get("expected_verdict", "UNKNOWN"),
                runtime_sec=480,
                error_msg="Timeout after 480s"
            )

        # Print immediate result
        print(f"    Status: {result.status}, Gates: {result.gates}")
        print(f"    Verdict: {result.verdict} (expected: {result.expected})")
        if result.status == "COMPLETED" and result.gates == "PASS":
            print(f"    p_boot={result.p_boot:.3f}, Λ={result.Lambda_obs:.2f}")
            print(f"    |R|_shared={result.r_shared:.3f}, φ={result.phi_shared:.1f}°")
        print(f"    Runtime: {result.runtime_sec:.1f}s")
        sys.stdout.flush()

        write_result_to_csv(result, csv_path, write_header=write_header)
        write_header = False
        results.append(result)

    # Load all results
    all_results = []
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            all_results.append(PrelimResult(
                family_id=row['family_id'],
                category=row['category'],
                status=row['status'],
                gates=row['gates'],
                chi2_dof_A=float(row['chi2_dof_A']),
                chi2_dof_B=float(row['chi2_dof_B']),
                r_shared=float(row['r_shared']),
                phi_shared=float(row['phi_shared']),
                Lambda_obs=float(row['Lambda_obs']),
                p_boot=float(row['p_boot']),
                n_boot=int(row['n_boot']),
                verdict=row['verdict'],
                expected=row['expected'],
                runtime_sec=float(row['runtime_sec']),
                error_msg=row.get('error_msg', '')
            ))

    generate_summary_md(all_results, output_dir)

    print("\n" + "=" * 60)
    print("FINAL RESULTS")
    print("=" * 60)
    print(f"{'Family':<25} {'Cat':<8} {'Verdict':<15} {'p_boot':<8} {'|R|':<6}")
    print("-" * 70)
    for r in all_results:
        print(f"{r.family_id:<25} {r.category:<8} {r.verdict:<15} {r.p_boot:.3f}    {r.r_shared:.3f}")

    print(f"\nResults: {csv_path}")
    print(f"Summary: {output_dir / 'PRELIM_SUMMARY_v2.md'}")


if __name__ == "__main__":
    main()
