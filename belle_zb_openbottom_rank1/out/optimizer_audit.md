# Optimizer Audit

## Analysis Method

For this analysis, we used **published Table I parameters** from Belle arXiv:1512.07419
rather than re-fitting the binned spectra. This approach was chosen because:

1. Direct spectral fitting failed fit-health gates (χ²/dof > 3)
2. Published parameters represent Belle's official analysis
3. Avoids normalization and phase-space modeling uncertainties

## Table-Based Analysis

### Input Parameters

From Belle Table I (BB*π, Model-2):

**Solution 1:**
| Parameter | Value | Error |
|-----------|-------|-------|
| f(Zb10610) | 1.01 | 0.13 |
| f(Zb10650) | 0.05 | 0.04 |
| φ(Zb10650) | -0.26 rad | 0.68 rad |

**Solution 2:**
| Parameter | Value | Error |
|-----------|-------|-------|
| f(Zb10610) | 1.18 | 0.15 |
| f(Zb10650) | 0.24 | 0.11 |
| φ(Zb10650) | -1.63 rad | 0.14 rad |

### Derived Quantities

Amplitude ratio approximated as |R| ≈ sqrt(f_2/f_1):

| Solution | |R| | σ(|R|) | arg(R) | σ(arg(R)) |
|----------|-----|--------|--------|---------|
| 1 | 0.222 | 0.090 | -15° | 39° |
| 2 | 0.451 | 0.107 | -93° | 8° |

## Cross-Family Test Results

### Solution 1 vs Υ-channel average

| Metric | Value |
|--------|-------|
| χ² (magnitude) | 24.96 |
| dof | 1 |
| p-value | 0.0000 |
| Verdict | DISFAVORED |

### Solution 2 vs Υ-channel average

| Metric | Value |
|--------|-------|
| χ² (magnitude) | 8.54 |
| dof | 1 |
| p-value | 0.0035 |
| Verdict | DISFAVORED |

## Stability Checks

### Λ ≥ 0 Constraint

Not applicable for table-based analysis (no likelihood ratio test performed).

### Multi-Start Consistency

Not applicable for table-based analysis (parameters taken directly from publication).

## Direct Spectral Fit Results (Failed)

For completeness, the direct spectral fits gave:

| Channel | Model | χ²/dof | Health |
|---------|-------|--------|--------|
| BB*π | Coherent two-BW | 18.72 | MODEL_MISMATCH |
| BB*π | Incoherent | 14414 | MODEL_MISMATCH |
| B*B*π | Single-BW | 27936 | MODEL_MISMATCH |

These failures indicate that the simple Breit-Wigner models without proper:
- Phase space factors
- Efficiency weighting
- Non-resonant contributions

cannot adequately describe the data. The published Belle analysis uses more
sophisticated models and fitting procedures.

## Conclusion

The table-based analysis provides the most reliable extraction of R from the
open-bottom data. The cross-family comparison shows tension with hidden-bottom
R values, which may indicate:

1. Physical differences between hidden and open-bottom decays
2. Threshold enhancement effects for Zb(10610) in BB* channel
3. Systematic effects in the approximation |R| ≈ sqrt(f_2/f_1)

---
*Generated for Belle open-bottom rank-1 analysis*
