#!/usr/bin/env python3
"""
Generate debug overlay plots showing extracted data vs original figures.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

def load_csv(filepath):
    """Load CSV data."""
    data = []
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('m_low'):
                continue
            parts = line.strip().split(',')
            if len(parts) >= 10:
                m_center = float(parts[2])
                signal = float(parts[8])
                signal_err = float(parts[9])
                data.append((m_center, signal, signal_err))
    return np.array(data)


def generate_overlay_plots():
    """Generate overlay plots."""
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    # Load data
    bb_star_data = load_csv("../extracted/bb_star_pi.csv")
    bsbs_data = load_csv("../extracted/b_star_b_star_pi.csv")

    # Create figure with two panels
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # BB*π channel
    ax = axes[0]
    m = bb_star_data[:, 0]
    y = bb_star_data[:, 1]
    yerr = bb_star_data[:, 2]

    ax.errorbar(m, y, yerr=yerr, fmt='o', color='blue', capsize=3, label='Extracted signal (RS-WS)')
    ax.axvline(10607.2, color='red', linestyle='--', alpha=0.7, label='Zb(10610)')
    ax.axvline(10652.2, color='green', linestyle='--', alpha=0.7, label='Zb(10650)')
    ax.axhline(0, color='gray', linestyle='-', alpha=0.3)
    ax.set_xlabel('$M_{miss}(\\pi)$ [MeV/$c^2$]')
    ax.set_ylabel('Background-subtracted events / 5 MeV')
    ax.set_title('BB*π channel - Extracted Data')
    ax.legend()
    ax.set_xlim(10580, 10730)

    # Annotate peak regions
    ax.annotate('Zb(10610)', xy=(10607, 30), fontsize=9)
    ax.annotate('Zb(10650)\n(small)', xy=(10650, 5), fontsize=9)

    # B*B*π channel
    ax = axes[1]
    m = bsbs_data[:, 0]
    y = bsbs_data[:, 1]
    yerr = bsbs_data[:, 2]

    ax.errorbar(m, y, yerr=yerr, fmt='s', color='red', capsize=3, label='Extracted signal (RS-WS)')
    ax.axvline(10607.2, color='red', linestyle='--', alpha=0.3, label='Zb(10610) - below threshold')
    ax.axvline(10652.2, color='green', linestyle='--', alpha=0.7, label='Zb(10650)')
    ax.axhline(0, color='gray', linestyle='-', alpha=0.3)
    ax.set_xlabel('$M_{miss}(\\pi)$ [MeV/$c^2$]')
    ax.set_ylabel('Background-subtracted events / 5 MeV')
    ax.set_title('B*B*π channel - Extracted Data (Zb(10650) only)')
    ax.legend()
    ax.set_xlim(10630, 10730)

    ax.annotate('Zb(10650)', xy=(10652, 40), fontsize=9)

    plt.tight_layout()
    plt.savefig('../out/debug_extracted_overlay.png', dpi=150)
    plt.close()
    print("Saved: ../out/debug_extracted_overlay.png")

    # Also create individual channel overlays
    for name, data, title in [
        ('bb_star_pi', bb_star_data, 'BB*π'),
        ('b_star_b_star_pi', bsbs_data, 'B*B*π')
    ]:
        fig, ax = plt.subplots(figsize=(8, 5))
        m = data[:, 0]
        y = data[:, 1]
        yerr = data[:, 2]

        ax.bar(m, y, width=4.5, alpha=0.7, color='steelblue', edgecolor='black', label='Signal (RS-WS)')
        ax.errorbar(m, y, yerr=yerr, fmt='none', color='black', capsize=2)
        ax.axvline(10607.2, color='red', linestyle='--', linewidth=2, label='Zb(10610)')
        ax.axvline(10652.2, color='green', linestyle='--', linewidth=2, label='Zb(10650)')
        ax.axhline(0, color='gray', linestyle='-', alpha=0.5)
        ax.set_xlabel('$M_{miss}(\\pi)$ [MeV/$c^2$]')
        ax.set_ylabel('Background-subtracted events / 5 MeV')
        ax.set_title(f'{title} channel - Data from Belle arXiv:1512.07419')
        ax.legend()

        plt.tight_layout()
        plt.savefig(f'../out/debug_{name}_overlay.png', dpi=150)
        plt.close()
        print(f"Saved: ../out/debug_{name}_overlay.png")


if __name__ == "__main__":
    generate_overlay_plots()
