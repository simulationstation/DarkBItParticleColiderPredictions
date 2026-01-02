# CMSSW Synthetic Feasibility Audit

Classification legend:

- **A) Direct synthetic resonance decay (good)**: two-state family represented as
  two resonances decaying to the channel final states.
- **B) Proxy synthetic (ok but labeled)**: line-shape / production mismatch with
  original experiment; still valid as a toy rank-1 identifiability test.
- **C) Not feasible (skip)**: requires detector simulation or external
  generator beyond CMSSW GEN-only.

## Family classifications

| Family | Category | CMSSW Synth Class | Notes |
|--------|----------|-------------------|-------|
| belle_zb | exotic | A | Synthetic Zb(10610/10650) → Υ(nS)π and h_b(mP)π. h_b PDG IDs assumed; if unavailable in Pythia, restrict to Υ channels only. |
| lhcb_pc_doublet | exotic | A | Synthetic Pc narrow doublet decays to J/ψ p. Production mode is toy-only. |
| strange_pcs | exotic | A | Synthetic Pcs states decaying to J/ψ Λ. |
| besiii_y_pipijpsi_hc | exotic | B | Proxy for e+e- line shapes; toy resonances decaying to ππJ/ψ and ππh_c. |
| besiii_belle_isr_y | exotic | B | Proxy for ISR line shapes; toy resonances decaying to ππJ/ψ and ππψ(2S). |
| control_babar_phi | control | B | Proxy for e+e- cross sections; toy resonances decaying to ϕ f0 and ϕ ππ. |
| control_babar_omega | control | B | Proxy for e+e- cross sections; toy resonances decaying to ω f0 and ω ππ. |
| x3872_like | exotic | A | Synthetic X(3872)-like states decaying to J/ψ ππ and J/ψ πππ0. |

## Exclusions

- `cms_atlas_dicharmonium_other`: removed per instructions (di-charmonium locked down).
- `cms_x6900_x7100`, `zc_like`: already handled elsewhere.
