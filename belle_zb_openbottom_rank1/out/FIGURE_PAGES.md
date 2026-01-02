# Belle Open-Bottom Zb Figure Pages

## Source Paper

- **arXiv**: 1512.07419
- **Title**: "Study of e+e−→B(∗)B̄(∗)π± at √s = 10.866 GeV"
- **Collaboration**: Belle

## Relevant Figure Pages

### Page 4 (Main spectra - Figure 2)
- **Content**: Mmiss(Bπ) and Mmiss(π) distributions
- **Panels**:
  - (a) M*miss(Bπ) for BB*π candidates
  - (b) M*miss(Bπ) for B*B*π candidates
  - (c) Mmiss(π) for BB*π signal region
  - (d) Mmiss(π) for B*B*π signal region
- **Rendered**: data/figures/page04_600dpi.png

### Page 5 (Table I - Fit results)
- **Content**: Summary of fit results with different models
- **Models**: Model-0 through Model-3 with various Zb combinations
- **Key parameters**: fZb(10610), fZb(10650), φZb(10650), fnr
- **Rendered**: data/figures/page05_600dpi.png

### Page 8 (Supplementary - Binned data)
- **Content**: Bin-by-bin Mmiss(π) distributions
- **Data**: Right-sign events, wrong-sign events, relative efficiency
- **Channels**: BB*π (26 bins), B*B*π (17 bins)
- **Rendered**: data/figures/page08_600dpi.png

## Data Extraction Method

**Primary method**: Direct extraction from Supplementary Table I (Page 8)

The paper provides binned data in the supplementary material, eliminating the need for figure digitization. This is the gold standard for extraction:

- Bin width: 5 MeV/c²
- Mass range: 10590-10720 MeV/c²
- Columns: RS events, WS events, relative efficiency
- Signal: RS - WS (background-subtracted)

**Extracted files**:
- extracted/bb_star_pi.csv
- extracted/b_star_b_star_pi.csv
- extracted/binned_data.json

## Physical Notes

1. **BB* threshold**: ~10604 MeV/c² (M_B + M_B* ≈ 5279 + 5325)
2. **B*B* threshold**: ~10650 MeV/c² (2 × M_B* ≈ 2 × 5325)
3. **Zb(10610) mass**: 10607.2 MeV/c² (just above BB* threshold)
4. **Zb(10650) mass**: 10652.2 MeV/c² (just above B*B* threshold)

This means:
- BB*π can probe both Zb(10610) and Zb(10650)
- B*B*π can only probe Zb(10650) since Zb(10610) is below the B*B* threshold

---
*Generated for Belle open-bottom rank-1 analysis*
