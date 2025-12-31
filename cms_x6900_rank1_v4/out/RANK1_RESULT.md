# CMS Rank-1 Bottleneck Test Results

**Version**: Publication Grade v2.0
**Date**: 2025-12-31 09:58:50

## Summary

| Metric | Value |
|--------|-------|
| **Verdict** | NOT_REJECTED |
| Reason | p=0.3987 >= 0.05 |
| NLL (constrained) | 16.87 |
| NLL (unconstrained) | 16.62 |
| Lambda | 0.50 |
| dof_diff | 2 |
| chi2(2) 95% threshold | 5.99 |

## P-values

| Method | p-value | Notes |
|--------|---------|-------|
| **Bootstrap (primary)** | 0.3987 | 319/800 exceedances |
| Minimum resolvable | 0.0013 | 1/N_bootstrap |
| Wilks (reference) | 0.7791 | chi2(2) approximation |

## Fit Health

| Channel | chi2/dof | deviance/dof | Status |
|---------|----------|--------------|--------|
| A | 1.21 | 1.19 | HEALTHY |
| B | 1.91 | 2.42 | HEALTHY |

*Thresholds: UNDERCONSTRAINED if chi2/dof < 0.5, MODEL_MISMATCH if chi2/dof > 3.0 or dev/dof > 3.0*

## Coupling Ratios

```
R_shared = 7.5094 * exp(i * 1.5131)
R_A      = 8.0087 * exp(i * 1.5079)
R_B      = 7.5924 * exp(i * 1.5543)
```

## Input Files

- Channel A: `/home/primary/DarkBItParticleColiderPredictions/cms_x6900_rank1_v4/extracted/channelA_trimmed.csv`
- Channel B: `/home/primary/DarkBItParticleColiderPredictions/cms_x6900_rank1_v4/extracted/channelB_jpsi_psi2S_bins.csv`
- Bootstrap replicates: 800
- Optimizer starts: 300

## Interpretation Guide

| Verdict | Meaning |
|---------|--------|
| NOT_REJECTED | Data consistent with shared R (rank-1) |
| DISFAVORED | Evidence against shared R |
| INCONCLUSIVE | Cannot draw conclusion (fit issues) |
| MODEL_MISMATCH | Model does not describe data |
| OPTIMIZER_FAILURE | Numerical issues in optimization |
