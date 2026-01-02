#!/usr/bin/env python3
"""Convert BESIII pi+pi- hc HEPData to standard CSV format."""

import csv
from pathlib import Path

def parse_hepdata_csv(filepath):
    """Parse HEPData CSV format extracting sqrt(s) and cross section."""
    data = []

    with open(filepath, 'r') as f:
        content = f.read()

    # Find the cross section block (last block in file)
    blocks = content.split('\n\n')

    for block in blocks:
        lines = block.strip().split('\n')
        if not lines:
            continue

        # Look for the dressed cross section block
        if 'dressed cross section' in lines[0].lower():
            for line in lines[1:]:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.split(',')
                if len(parts) >= 4:
                    try:
                        sqrt_s = float(parts[0])
                        sigma = float(parts[1])
                        stat_plus = float(parts[2])
                        stat_minus = abs(float(parts[3]))
                        # Average asymmetric errors
                        stat_err = (stat_plus + stat_minus) / 2

                        # Systematic errors if available
                        sys_err = 0
                        if len(parts) >= 6:
                            sys_plus = float(parts[4]) if parts[4] else 0
                            sys_minus = abs(float(parts[5])) if parts[5] else 0
                            sys_err = (sys_plus + sys_minus) / 2

                        # Skip negative cross sections
                        if sigma > 0:
                            data.append({
                                'sqrt_s': sqrt_s,
                                'sigma': sigma,
                                'stat_err': stat_err,
                                'sys_err': sys_err
                            })
                    except (ValueError, IndexError):
                        continue

    return data


def main():
    base = Path(__file__).parent.parent
    hepdata_dir = base / 'data' / 'hepdata'
    output_dir = base / 'extracted'
    output_dir.mkdir(exist_ok=True)

    # Combine all hc tables
    all_data = []

    for table_file in ['table3_part1_hc.csv', 'table3_part2_hc.csv', 'table4_hc.csv']:
        filepath = hepdata_dir / table_file
        if filepath.exists():
            data = parse_hepdata_csv(filepath)
            all_data.extend(data)
            print(f"Parsed {len(data)} points from {table_file}")

    # Sort by sqrt_s
    all_data.sort(key=lambda x: x['sqrt_s'])

    # Remove duplicates (keep first occurrence)
    seen = set()
    unique_data = []
    for d in all_data:
        key = round(d['sqrt_s'], 4)
        if key not in seen:
            seen.add(key)
            unique_data.append(d)

    # Write output
    output_file = output_dir / 'channelB_hc_xsec.csv'
    with open(output_file, 'w') as f:
        f.write("# BESIII e+e- -> pi+pi- hc cross section (HEPData ins2908630)\n")
        f.write("# Channel B: pi+pi- hc\n")
        f.write("# sqrt_s (GeV), sigma (pb), stat_err, sys_err\n")
        for d in unique_data:
            f.write(f"{d['sqrt_s']:.4f},{d['sigma']:.1f},{d['stat_err']:.1f},{d['sys_err']:.1f}\n")

    print(f"\nWrote {len(unique_data)} points to {output_file}")

    # Print energy overlap with J/psi channel
    jpsi_energies = [3.773, 3.808, 3.896, 4.008, 4.086, 4.189, 4.208, 4.217,
                     4.226, 4.242, 4.258, 4.308, 4.358, 4.387, 4.416, 4.467,
                     4.527, 4.575, 4.600]

    hc_energies = [d['sqrt_s'] for d in unique_data]

    # Find overlapping energy points (within 10 MeV)
    overlaps = []
    for e_jpsi in jpsi_energies:
        for e_hc in hc_energies:
            if abs(e_jpsi - e_hc) < 0.015:
                overlaps.append((e_jpsi, e_hc))

    print(f"\nOverlapping energy points (within 15 MeV): {len(overlaps)}")
    for e_jpsi, e_hc in overlaps[:10]:
        print(f"  J/psi: {e_jpsi:.3f} GeV, hc: {e_hc:.4f} GeV")


if __name__ == '__main__':
    main()
