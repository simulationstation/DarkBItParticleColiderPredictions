#!/usr/bin/env python3
"""Convert HEPData CSV to harness format."""
import csv

input_file = '../data/hepdata/figure1_spectrum.csv'
output_file = '../data/hepdata/channelA_dijpsi.csv'

with open(input_file, 'r') as f:
    lines = f.readlines()

# Skip header lines starting with #
data_lines = [l for l in lines if not l.startswith('#')]

# Parse CSV
reader = csv.reader(data_lines)
header = next(reader)

# Output format: mass_GeV, counts, stat_err
with open(output_file, 'w') as out:
    out.write("mass_GeV,counts,stat_err\n")
    for row in reader:
        if len(row) >= 6:
            m_center = float(row[0])
            m_low = float(row[1])
            m_high = float(row[2])
            counts = float(row[3])
            stat_plus = abs(float(row[4])) if row[4] else 0
            stat_minus = abs(float(row[5])) if row[5] else 0
            stat_err = max(stat_plus, stat_minus, 1.0)  # Min error of 1 for Poisson
            out.write(f"{m_center:.4f},{int(counts)},{stat_err:.4f}\n")

print(f"Converted {output_file}")
