"""Family-specific GEN analyzer wrapper (SYNTHETIC, DRY_RUN by default)."""
from __future__ import annotations

import argparse
from pathlib import Path

from EXOTICS_FACTORY.cmssw_synth.analyzers import gen_mass_spectra


def main() -> None:
    parser = argparse.ArgumentParser(description="control_babar_omega GEN-only analyzer (SYNTHETIC)")
    parser.add_argument("--spec", default="EXOTICS_FACTORY/families/control_babar_omega/cmssw_synth/spec_cmssw.yaml")
    parser.add_argument("--input", default="control_babar_omega_GEN.root", help="Input GEN ROOT file")
    parser.add_argument(
        "--output-dir",
        default="EXOTICS_FACTORY/out/cmssw_synth/control_babar_omega",
        help="Output directory for synthetic CSVs",
    )
    parser.add_argument("--run", action="store_true", help="Execute analysis (default is dry-run)")
    args = parser.parse_args()

    if not args.run:
        print("DRY_RUN: No analysis executed.")
        print(f"Spec: {args.spec}")
        print(f"Input: {args.input}")
        print(f"Output dir: {args.output_dir}")
        print("Re-run with --run to execute.")
        return

    gen_mass_spectra.run_analysis(
        input_root=Path(args.input),
        output_dir=Path(args.output_dir),
        channels=gen_mass_spectra.parse_channels(gen_mass_spectra.load_spec(Path(args.spec))),
    )


if __name__ == "__main__":
    main()
