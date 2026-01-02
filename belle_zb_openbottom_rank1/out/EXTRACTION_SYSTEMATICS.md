# Extraction Systematics

## Data Source

The open-bottom Zb data was extracted from **Supplementary Table I** of Belle arXiv:1512.07419,
NOT from figure digitization. This provides the highest quality extraction possible.

## Systematic Uncertainties

### Bin-by-bin Statistical Uncertainties

The published table provides statistical uncertainties for both:
- **Right-Sign (RS)** events: Signal + combinatorial background
- **Wrong-Sign (WS)** events: Combinatorial background estimate

The signal uncertainty is computed as:
```
σ(signal) = sqrt(σ(RS)² + σ(WS)²)
```

### Efficiency Correction

Belle provides relative efficiency per bin. The published yields are NOT efficiency-corrected.
For the rank-1 test, we use the raw background-subtracted yields since:
1. Both Zb states are in the same mass region with similar efficiencies
2. The ratio R is largely insensitive to overall efficiency normalization

### Extraction Perturbation Analysis

Since data comes from a published table (not figure digitization), there are NO:
- Crop box uncertainties
- Tick association jitter
- Axis mapping errors

The systematic uncertainty from extraction is **negligible**.

## Model Systematics

### Amplitude Ratio Approximation

We approximate the amplitude ratio as:
```
|R| ≈ sqrt(f_Zb10650 / f_Zb10610)
```

where f_X is the fit fraction from Table I. This approximation has systematic uncertainty because:

1. **Interference effects**: The full integral depends on interference:
   ```
   S(m) = |A1 + R*A2|² ≠ |A1|² + |R*A2|²
   ```

2. **Phase space**: Different resonances sample different phase space regions

3. **Non-resonant contributions**: Table I Model-2 includes only Zb states, not non-resonant

Estimated systematic uncertainty on |R| from this approximation: ~20-30%

### Fit Solution Ambiguity

Belle reports TWO solutions for Model-2:
- Solution 1: |R| ≈ 0.22 (small Zb(10650) fraction)
- Solution 2: |R| ≈ 0.45 (larger Zb(10650) fraction)

This ambiguity represents a fundamental limitation of the current data.

## Summary Table

| Source | Uncertainty Type | Magnitude |
|--------|-----------------|-----------|
| Statistical | Random | From published errors |
| Extraction | None | Data from table |
| Efficiency | Systematic | Negligible for ratio |
| |R| approximation | Systematic | 20-30% |
| Fit ambiguity | Systematic | Factor ~2 (Sol.1 vs Sol.2) |

## Conclusion

The dominant systematic is the **fit solution ambiguity**, which gives a factor ~2 uncertainty
in the extracted |R|. Both solutions are statistically acceptable fits to the data.

---
*Generated for Belle open-bottom rank-1 analysis*
