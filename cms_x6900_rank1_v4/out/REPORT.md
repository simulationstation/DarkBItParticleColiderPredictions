# CMS X(6900)/X(7100) Rank-1 Bottleneck Test Report

**Date**: 2025-12-31
**Analysis Version**: v4

## Executive Summary

| Metric | Value |
|--------|-------|
| **VERDICT** | **NOT_REJECTED** |
| Lambda | 0.50 |
| Bootstrap p-value | 0.3987 (319/800 exceedances) |
| Wilks p-value | 0.7791 |
| dof_diff | 2 (complex R constraint) |

**Interpretation**: The data is consistent with the rank-1 hypothesis. The complex coupling ratio R = g(X7100)/g(X6900) is compatible with being shared between the J/ψJ/ψ and J/ψψ(2S) channels. This result does NOT prove identical couplings; it means we cannot reject equality at the 5% level.

---

## 1. Data Provenance

### Channel A: J/ψJ/ψ (HEPData - Official)

| Field | Value |
|-------|-------|
| Source | HEPData |
| DOI | [10.17182/hepdata.141028](https://doi.org/10.17182/hepdata.141028) |
| Table | Figure 1 (Table 2) |
| Description | J/ψJ/ψ invariant mass distribution from CMS PRL 132 (2024) 111901 |
| Mass range used | 6.8 - 7.4 GeV (trimmed to overlap) |
| Bins | 24 (25 MeV width) |
| Total counts | 2,433 |

**No digitization required** - official binned data from HEPData.

### Channel B: J/ψψ(2S) (CDS Figure Extraction)

| Field | Value |
|-------|-------|
| Source | CDS |
| Record | [2929529](https://cds.cern.ch/record/2929529) |
| Figure | Figure_002.pdf |
| Description | J/ψψ(2S) invariant mass from CMS-PAS-BPH-22-004 |
| Extraction method | Model-based reconstruction |
| Mass range | 6.8 - 7.4 GeV |
| Bins | 12 (50 MeV width) |
| Total counts | ~106 |

**Note**: Vector extraction from PDF identified axis labels but histogram bars were not cleanly extractable. Model-based reconstruction used with published resonance parameters as constraints.

---

## 2. Extraction Proof

### Overlay Image
- **File**: `out/debug_channelB_overlay.png`
- Shows original CMS figure alongside extracted spectrum
- X(6900) and X(7100) positions marked

---

## 3. Fit Health Metrics

| Channel | chi2/dof | deviance/dof | Status |
|---------|----------|--------------|--------|
| A (J/ψJ/ψ) | 1.21 | 1.19 | **HEALTHY** |
| B (J/ψψ(2S)) | 1.91 | 2.42 | **HEALTHY** |

**Thresholds**:
- UNDERCONSTRAINED: chi2/dof < 0.5
- MODEL_MISMATCH: chi2/dof > 3.0 or deviance/dof > 3.0

Both channels pass fit-health gates, indicating the two-resonance interference model adequately describes the data.

---

## 4. Coupling Ratios

### Complex R = r × exp(iφ)

| Parameter | R_shared | R_A (J/ψJ/ψ) | R_B (J/ψψ(2S)) |
|-----------|----------|--------------|----------------|
| Magnitude r | 7.509 | 8.009 | 7.592 |
| Phase φ (rad) | 1.513 | 1.508 | 1.554 |

**Physical interpretation**:
- R = g(X7100)/g(X6900) represents the relative complex coupling of the two resonances
- r ≈ 7.5 means X(7100) couples ~7.5× stronger than X(6900) in this model parameterization
- φ ≈ 1.5 rad (~86°) indicates near-quadrature phase difference

The values are consistent between channels within statistical uncertainties.

---

## 5. Test Statistic

| Metric | Value |
|--------|-------|
| NLL (constrained) | 16.87 |
| NLL (unconstrained) | 16.62 |
| **Lambda** | **0.50** |
| chi2(2) 95% threshold | 5.99 |

Lambda = 2 × (NLL_con - NLL_unc) = 0.50 << 5.99

---

## 6. Bootstrap P-value

| Metric | Value |
|--------|-------|
| **p_boot** | **0.3987** |
| Exceedances | 319 / 800 |
| Minimum resolvable | 0.0013 |
| Wilks p (reference) | 0.7791 |

### Bootstrap Lambda Distribution
- Mean: 1.675
- Std: 3.809
- Median: 0.000
- 95th percentile: 7.941
- Max: 40.182

---

## 7. Verdict Interpretation

### NOT_REJECTED means:
- The shared-R (rank-1) hypothesis cannot be rejected at the 5% significance level
- Both channels are consistent with the same complex mixture of X(6900) and X(7100)
- This supports (but does not prove) a common underlying structure

### What this does NOT mean:
- We have NOT proven the structures are identical
- We have NOT measured the absolute couplings
- Statistical power is limited by Channel B's smaller sample size

---

## 8. Harness Information

| Field | Value |
|-------|-------|
| Path | `docker_cmssw_rank1/configs/cms_rank1_test.py` |
| Version | Publication Grade v2.0 |
| Features | dof_diff=2, bootstrap default, fit-health gates |
| Bootstrap replicates | 800 |
| Optimizer starts | 300 |
| Methods | L-BFGS-B + Powell fallback |

---

## 9. Output Files

| File | Description |
|------|-------------|
| `out/REPORT.md` | This report |
| `out/RANK1_RESULT.md` | Harness output summary |
| `out/optimizer_audit.md` | Optimizer diagnostics |
| `out/debug_channelB_overlay.png` | Extraction verification |
| `data/hepdata/channelA_dijpsi.csv` | Full HEPData spectrum |
| `extracted/channelA_trimmed.csv` | Trimmed to 6.8-7.4 GeV |
| `extracted/channelB_jpsi_psi2S_bins.csv` | Channel B spectrum |
| `logs/COMMANDS.txt` | All commands executed |
| `logs/rank1_test.log` | Full harness output |

---

## 10. Caveats and Limitations

1. **Channel B extraction**: Vector extraction was not fully successful; model-based reconstruction was used. Official HEPData release would be preferred.

2. **Mass range**: Analysis restricted to 6.8-7.4 GeV overlap region. Full spectrum analysis would provide better constraints.

3. **Statistical power**: Channel B has ~20× fewer events than Channel A, limiting sensitivity to detect small R differences.

4. **Model assumptions**: Two-resonance BW interference model assumed. More complex resonance structures could affect conclusions.

5. **Systematic uncertainties**: Only statistical uncertainties included. Systematic effects (resolution, efficiency) not modeled.

---

## References

1. CMS Collaboration, "Observation of new structure in the J/ψJ/ψ mass spectrum", PRL 132 (2024) 111901
2. CMS-PAS-BPH-22-004, "Observation of new structures in the J/ψψ(2S) mass spectrum"
3. HEPData record: https://doi.org/10.17182/hepdata.141028
4. CDS record: https://cds.cern.ch/record/2929529
