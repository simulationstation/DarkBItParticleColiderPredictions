# BESIII Y(4220)/Y(4320) Rank-1 Factorization Test - Comprehensive Report

## Executive Summary

**Verdict: INCONCLUSIVE (Fit Health Gates FAILED)**

The rank-1 factorization test could not reach a definitive conclusion because the two-resonance model does not adequately describe the data in the energy window 4.15-4.45 GeV. The chi²/dof values are extremely poor (24.9 for J/ψ, 9.8 for hc), indicating model mismatch.

---

## Test Results

| Metric | Value | Notes |
|--------|-------|-------|
| Lambda_obs | 331.6274 | Extremely high - would strongly reject if fit were healthy |
| p_boot | 0.1350 (27/200) | Surprisingly high given Lambda |
| p_chi2(2) | 0.0000 | Expected given huge Lambda |
| chi²/dof (A) | 24.930 | **FAR OUTSIDE healthy range [0.3, 4.0]** |
| chi²/dof (B) | 9.774 | **FAR OUTSIDE healthy range [0.3, 4.0]** |
| Gates | **FAIL** | Model mismatch detected |
| Verdict | **INCONCLUSIVE** | Cannot draw physics conclusions |

---

## Coupling Ratio Estimates (from poor fit)

### Constrained (Shared R):
- |R| = 0.666
- arg(R) = -146.3°

### Unconstrained (Independent R per channel):
| Channel | |R| | arg(R) |
|---------|-----|--------|
| A (π+π-J/ψ) | 0.435 | 174.2° |
| B (π+π-hc) | 1.911 | 83.0° |

The coupling ratios differ significantly, but this result is unreliable due to model mismatch.

---

## Data Summary

### Channel A: e+e- → π+π- J/ψ (arXiv:1611.01317)
- Source: BESIII 2017 paper, Table III
- Points in window [4.15, 4.45] GeV: **10 points**

```
sqrt_s (GeV), sigma (pb), stat_err, sys_err
4.1886,15.5,3.8,0.9
4.2077,53.4,5.4,3.1
4.2171,60.3,5.7,3.5
4.2263,85.1,1.5,4.9
4.2417,84.4,6.3,4.9
4.2580,59.5,1.4,3.4
4.3079,52.0,5.7,3.0
4.3583,25.4,1.2,1.5
4.3874,20.0,3.2,1.2
4.4156,12.1,0.6,0.7
```

### Channel B: e+e- → π+π- hc (HEPData ins2908630)
- Source: BESIII 2025 HEPData record
- Points in window [4.15, 4.45] GeV: **30 points**

```
sqrt_s (GeV), sigma (pb), stat_err, sys_err
4.1570,5.7,1.9,0.7
4.1780,13.8,0.8,0.8
4.1890,17.7,2.1,1.6
4.1990,21.3,2.3,2.8
4.2080,39.0,8.6,2.7
4.2090,29.1,2.7,1.6
4.2170,46.6,9.2,2.6
4.2190,42.4,3.0,2.4
4.2260,46.3,2.1,2.6
4.2360,43.9,2.9,2.4
4.2420,37.3,8.7,2.1
4.2440,40.7,2.8,2.3
4.2580,38.9,2.2,2.2
4.2670,39.8,2.8,2.2
4.2780,37.0,4.8,2.1
4.2870,36.7,2.9,2.0
4.3080,47.2,9.9,3.6
4.3110,39.7,3.1,2.2
4.3370,43.7,3.1,2.4
4.3580,48.2,2.9,3.0
4.3770,46.7,3.0,2.9
4.3870,44.8,8.9,3.3
4.3950,46.6,3.1,2.9
4.4160,42.6,2.1,2.7
4.4180,34.6,21.2,2.2
4.4230,93.5,31.3,5.8
4.4280,67.5,27.8,4.2
4.4360,42.9,2.9,2.7
4.4380,40.2,22.5,2.5
4.4480,3.7,17.3,0.2
```

---

## Model Description

### Two-Resonance Breit-Wigner Model

The cross section is modeled as:

```
σ(√s) = norm × |A₁ + R × A₂|²

where:
  A₁ = M₁Γ₁ / (s - M₁² + iM₁Γ₁)  [Y(4220) amplitude]
  A₂ = M₂Γ₂ / (s - M₂² + iM₂Γ₂)  [Y(4320) amplitude]
  R = |R| × exp(iφ)               [complex coupling ratio]
```

### Fixed/Fitted Parameters

| Parameter | Initial Value | Fit Range | Source |
|-----------|---------------|-----------|--------|
| M₁ (Y4220) | 4.222 GeV | [4.15, 4.28] GeV | arXiv:1611.01317 |
| Γ₁ (Y4220) | 0.044 GeV | [0.02, 0.15] GeV | arXiv:1611.01317 |
| M₂ (Y4320) | 4.320 GeV | [4.28, 4.40] GeV | arXiv:1611.01317 |
| Γ₂ (Y4320) | 0.101 GeV | [0.03, 0.20] GeV | arXiv:1611.01317 |

---

## NLL Values

| Model | NLL | n_params |
|-------|-----|----------|
| Constrained (shared R) | 184.49 | 8 |
| Unconstrained (independent R) | 18.68 | 10 |

ΔnLL = 165.81, which gives Lambda = 2×ΔnLL = 331.62

---

## Potential Issues Identified

### 1. Model Mismatch
The simple two-BW interference model doesn't capture the lineshape. Possible reasons:
- **Additional resonances**: Y(4260) fine structure, Y(4360), Y(4390)
- **Non-resonant background**: Continuum contribution not modeled
- **Energy-dependent widths**: Widths may vary with √s
- **Interference with ψ(3770)**: Not included in current model

### 2. Different Physics in Two Channels
The π+π-J/ψ and π+π-hc channels may have:
- Different intermediate states (Zc(3900) vs Zc(4020))
- Different phase space factors
- Different detector systematics

### 3. Energy Point Mismatch
- Channel A: High-luminosity XYZ data at specific energies
- Channel B: Denser scan data covering more energies
- Not all energies overlap perfectly

### 4. Line Shape Differences
From the BESIII papers:
- **π+π-J/ψ (arXiv:1611.01317)**: Two structures at 4222 MeV (Γ=44) and 4320 MeV (Γ=101)
- **π+π-hc (HEPData 2025)**: Three structures at 4224 MeV, 4327 MeV, and 4467 MeV

The presence of a THIRD resonance in the hc channel near 4467 MeV could explain the model failure.

---

## Comparison with Preliminary Simulation

From EXOTICS_FACTORY prelim results with simulated data:
- **Verdict**: NOT_REJECTED (p_boot = 0.07, Λ = 4.60)
- **chi²/dof**: 0.72 and 0.87 (healthy)

The simulation used:
- Y(4220): M=4.218 GeV, Γ=0.066 GeV
- Y(4320): M=4.320 GeV, Γ=0.101 GeV
- True R_ratio: |R|=0.7, φ=-35°

The real data behaves very differently from the simulated benchmark.

---

## Raw Console Output

```
============================================================
BESIII Y(4220)/Y(4320) Rank-1 Factorization Test
============================================================

Channel A (pi+pi- J/psi): 10 points
Channel B (pi+pi- hc): 30 points
Energy window: 4.15 - 4.45 GeV

[1/3] Fitting constrained model (shared R)...
[2/3] Fitting unconstrained model (independent R)...

--- Fit Results ---
NLL constrained:   184.49
NLL unconstrained: 18.68
Lambda_obs:        331.6274

Shared |R|: 0.666
Shared arg(R): -2.55 rad (-146.3 deg)

Channel A |R|: 0.435, arg(R): 174.2 deg
Channel B |R|: 1.911, arg(R): 83.0 deg

--- Fit Health ---
chi2/dof (A): 24.930
chi2/dof (B): 9.774

[3/3] Running bootstrap (200 replicates)...

============================================================
RESULT
============================================================
Lambda_obs: 331.6274
p_boot:     0.1350 (27/200 exceedances)
p_chi2(2):  0.0000
Gates:      FAIL
Verdict:    INCONCLUSIVE
```

---

## Data Sources

| Channel | Source | Record/Paper |
|---------|--------|--------------|
| π+π-J/ψ | arXiv | [1611.01317](https://arxiv.org/abs/1611.01317) |
| π+π-hc | HEPData | [ins2908630](https://www.hepdata.net/record/ins2908630) |

---

## Files Generated

```
besiii_y_rank1/
├── data/
│   ├── hepdata/
│   │   ├── table3_part1_hc.csv   # HEPData download
│   │   ├── table3_part2_hc.csv   # HEPData download
│   │   └── table4_hc.csv         # HEPData download
│   └── pdf/
│       └── besiii_jpsi_1611.01317.pdf  # Original paper
├── extracted/
│   ├── channelA_jpsi_xsec.csv    # Processed J/psi data
│   └── channelB_hc_xsec.csv      # Processed hc data
├── src/
│   ├── convert_hepdata_hc.py     # HEPData conversion script
│   └── besiii_y_rank1_test.py    # Main analysis script
├── out/
│   └── run.log                   # Full console output
└── REPORT.md                     # This report
```

---

## Recommendations for GPT Diagnosis

1. **Check if two-BW model is appropriate** for cross-section line shapes
2. **Consider adding ψ(3770) interference** as done in original paper
3. **Try three-resonance model** for hc channel (includes Y(4467))
4. **Match energy points** between channels more carefully
5. **Add continuum/background term** to improve fit quality
6. **Verify normalization** of cross sections between channels

---

## Conclusion

The rank-1 factorization test is **INCONCLUSIVE** because the simple two-Breit-Wigner model fails to describe either channel's lineshape (chi²/dof = 24.9 and 9.8). Before any physics conclusions can be drawn about whether Y(4220) and Y(4320) share common production mechanisms, the model must be improved to achieve reasonable fit quality.

The coupling ratios extracted from the poor fits show dramatic differences between channels (|R|=0.44 vs 1.91), but this cannot be trusted until the model describes the data adequately.
