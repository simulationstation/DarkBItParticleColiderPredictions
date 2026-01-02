#!/usr/bin/env python3
"""
Extract binned data from Belle arXiv:1512.07419 supplementary Table I.
Data represents Mmiss(π) distributions for BB*π and B*B*π final states.
"""

import os
import json

# Binned data from Supplementary Table I of arXiv:1512.07419
# Format: (bin_low, bin_high, RS_events, RS_err, WS_events, WS_err, efficiency)
# RS = Right-Sign (signal + background), WS = Wrong-Sign (background estimate)

BB_STAR_PI_DATA = [
    (10590, 10595, 1.0, 1.0, 2.3, 1.7, 0.937),
    (10595, 10600, 15.0, 3.9, 4.7, 2.4, 0.958),
    (10600, 10605, 38.0, 6.2, 14.3, 4.1, 0.979),
    (10605, 10610, 64.0, 8.0, 28.6, 5.8, 0.999),
    (10610, 10615, 81.0, 9.0, 46.4, 7.4, 1.017),
    (10615, 10620, 75.0, 8.7, 48.8, 7.6, 1.034),
    (10620, 10625, 81.0, 9.0, 45.2, 7.3, 1.050),
    (10625, 10630, 65.0, 8.1, 41.7, 7.0, 1.064),
    (10630, 10635, 47.0, 6.9, 44.1, 7.2, 1.075),
    (10635, 10640, 45.0, 6.7, 23.8, 5.3, 1.084),
    (10640, 10645, 39.0, 6.2, 35.7, 6.5, 1.089),
    (10645, 10650, 43.0, 6.6, 39.3, 6.8, 1.091),
    (10650, 10655, 40.0, 6.3, 34.5, 6.4, 1.089),
    (10655, 10660, 31.0, 5.6, 26.2, 5.6, 1.082),
    (10660, 10665, 35.0, 5.9, 30.9, 6.1, 1.070),
    (10665, 10670, 40.0, 6.3, 29.8, 5.9, 1.051),
    (10670, 10675, 27.0, 5.2, 19.0, 4.7, 1.025),
    (10675, 10680, 34.0, 5.8, 29.8, 5.9, 0.991),
    (10680, 10685, 26.0, 5.1, 27.4, 5.7, 0.946),
    (10685, 10690, 32.0, 5.7, 20.2, 4.9, 0.891),
    (10690, 10695, 17.0, 4.1, 9.5, 3.4, 0.821),
    (10695, 10700, 18.0, 4.2, 17.8, 4.6, 0.735),
    (10700, 10705, 8.0, 2.8, 9.5, 3.4, 0.629),
    (10705, 10710, 8.0, 2.8, 10.7, 3.6, 0.495),
    (10710, 10715, 6.0, 2.4, 5.9, 2.7, 0.322),
    (10715, 10720, 2.0, 1.4, 0.0, 1.1, 0.056),
]

# B*B*π data (starts at higher mass due to B*B* threshold)
B_STAR_B_STAR_PI_DATA = [
    (10635, 10640, 1.0, 1.0, 0.0, 1.1, 1.084),
    (10640, 10645, 11.0, 3.3, 0.0, 1.1, 1.089),
    (10645, 10650, 21.0, 4.6, 10.8, 3.6, 1.091),
    (10650, 10655, 55.0, 7.4, 9.6, 3.4, 1.089),
    (10655, 10660, 45.0, 6.7, 32.2, 6.2, 1.082),
    (10660, 10665, 46.0, 6.8, 33.4, 6.3, 1.070),
    (10665, 10670, 51.0, 7.1, 40.5, 6.9, 1.051),
    (10670, 10675, 43.0, 6.6, 21.4, 5.0, 1.025),
    (10675, 10680, 32.0, 5.7, 27.5, 5.7, 0.991),
    (10680, 10685, 21.0, 4.6, 16.7, 4.4, 0.946),
    (10685, 10690, 18.0, 4.2, 16.7, 4.4, 0.891),
    (10690, 10695, 20.0, 4.5, 14.3, 4.1, 0.821),
    (10695, 10700, 11.0, 3.3, 15.5, 4.3, 0.735),
    (10700, 10705, 9.0, 3.0, 7.1, 2.9, 0.629),
    (10705, 10710, 5.0, 2.2, 3.6, 2.1, 0.495),
    (10710, 10715, 5.0, 2.2, 2.4, 1.7, 0.322),
    (10715, 10720, 1.0, 1.0, 0.0, 1.1, 0.056),
]

def write_csv(data, filename, channel_name):
    """Write binned data to CSV file."""
    with open(filename, 'w') as f:
        f.write(f"# Belle arXiv:1512.07419 - {channel_name} channel\n")
        f.write("# Mmiss(π) distribution, bin width = 5 MeV/c²\n")
        f.write("# RS = Right-Sign (signal+background), WS = Wrong-Sign (background)\n")
        f.write("# Signal = RS - WS (background-subtracted)\n")
        f.write("m_low_MeV,m_high_MeV,m_center_MeV,RS,RS_err,WS,WS_err,efficiency,signal,signal_err\n")

        for row in data:
            m_low, m_high, rs, rs_err, ws, ws_err, eff = row
            m_center = (m_low + m_high) / 2.0

            # Background-subtracted signal
            signal = rs - ws
            signal_err = (rs_err**2 + ws_err**2)**0.5

            f.write(f"{m_low:.1f},{m_high:.1f},{m_center:.1f},{rs:.1f},{rs_err:.1f},"
                   f"{ws:.1f},{ws_err:.1f},{eff:.3f},{signal:.1f},{signal_err:.2f}\n")

    print(f"Written: {filename}")


def main():
    os.makedirs("../extracted", exist_ok=True)

    # Write BB*π data
    write_csv(BB_STAR_PI_DATA, "../extracted/bb_star_pi.csv", "BB*π")

    # Write B*B*π data
    write_csv(B_STAR_B_STAR_PI_DATA, "../extracted/b_star_b_star_pi.csv", "B*B*π")

    # Also save raw data as JSON for reference
    data_dict = {
        "source": "Belle Collaboration arXiv:1512.07419",
        "title": "Study of e+e−→B(∗)B̄(∗)π± at √s = 10.866 GeV",
        "supplementary_table": "Table I",
        "bin_width_MeV": 5,
        "columns": ["m_low_MeV", "m_high_MeV", "RS_events", "RS_err", "WS_events", "WS_err", "efficiency"],
        "BB_star_pi": BB_STAR_PI_DATA,
        "B_star_B_star_pi": B_STAR_B_STAR_PI_DATA,
        "physics_note": "B*B* threshold (~10650 MeV) is above Zb(10610) mass, so Zb(10610)->B*B* is suppressed"
    }

    with open("../extracted/binned_data.json", "w") as f:
        json.dump(data_dict, f, indent=2)

    print("Written: ../extracted/binned_data.json")

    # Summary statistics
    print("\n=== Data Summary ===")
    bb_star_signal = sum(row[2] - row[4] for row in BB_STAR_PI_DATA)
    bsbs_signal = sum(row[2] - row[4] for row in B_STAR_B_STAR_PI_DATA)

    print(f"BB*π: {len(BB_STAR_PI_DATA)} bins, total signal ≈ {bb_star_signal:.0f} events")
    print(f"B*B*π: {len(B_STAR_B_STAR_PI_DATA)} bins, total signal ≈ {bsbs_signal:.0f} events")

    # Identify Zb region signal
    zb10610_region = sum(row[2] - row[4] for row in BB_STAR_PI_DATA if 10590 <= row[0] <= 10630)
    zb10650_region_bb = sum(row[2] - row[4] for row in BB_STAR_PI_DATA if 10635 <= row[0] <= 10670)
    zb10650_region_bsbs = sum(row[2] - row[4] for row in B_STAR_B_STAR_PI_DATA if 10635 <= row[0] <= 10670)

    print(f"\nZb(10610) region (10590-10630) in BB*π: ~{zb10610_region:.0f} signal events")
    print(f"Zb(10650) region (10635-10670) in BB*π: ~{zb10650_region_bb:.0f} signal events")
    print(f"Zb(10650) region (10635-10670) in B*B*π: ~{zb10650_region_bsbs:.0f} signal events")


if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    main()
