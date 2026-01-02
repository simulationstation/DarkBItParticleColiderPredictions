"""GEN-only mass spectrum analyzer (SYNTHETIC template).

This module reads a CMSSW GEN ROOT file, reconstructs invariant masses
from final-state GenParticles, and writes binned CSV spectra in the
rank-1 harness format.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import awkward as ak
import numpy as np
import uproot
import yaml


@dataclass(frozen=True)
class Binning:
    bins: int
    min: float
    max: float


@dataclass(frozen=True)
class PdgGroup:
    ids: tuple[int, ...]
    abs: bool


@dataclass(frozen=True)
class ChannelSpec:
    name: str
    label: str
    pdg_groups: tuple[PdgGroup, ...]
    combination: str
    binning: Binning
    output_csv: str


def load_spec(path: Path) -> dict[str, Any]:
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def resolve_branch_map(tree: uproot.behaviors.TBranch.TBranch, prefix: str | None = None) -> dict[str, str]:
    keys = set(tree.keys())
    prefixes = [prefix] if prefix else []
    prefixes.extend([
        "GenPart",
        "genParticles",
        "prunedGenParticles",
        "packedGenParticles",
    ])

    for candidate in prefixes:
        if not candidate:
            continue
        pdg_key = f"{candidate}_pdgId"
        if pdg_key in keys:
            return {
                "pdgId": pdg_key,
                "status": f"{candidate}_status",
                "px": f"{candidate}_px",
                "py": f"{candidate}_py",
                "pz": f"{candidate}_pz",
                "energy": f"{candidate}_energy",
                "pt": f"{candidate}_pt",
                "eta": f"{candidate}_eta",
                "phi": f"{candidate}_phi",
                "mass": f"{candidate}_mass",
            }

    for key in keys:
        if not key.endswith("pdgId"):
            continue
        base = key[: -len("pdgId")]
        return {
            "pdgId": f"{base}pdgId",
            "status": f"{base}status",
            "px": f"{base}px",
            "py": f"{base}py",
            "pz": f"{base}pz",
            "energy": f"{base}energy",
            "pt": f"{base}pt",
            "eta": f"{base}eta",
            "phi": f"{base}phi",
            "mass": f"{base}mass",
        }

    raise KeyError("Unable to resolve GenParticle branches; provide a prefix override.")


def build_vectors(tree: uproot.behaviors.TBranch.TBranch, prefix: str | None = None) -> dict[str, ak.Array]:
    branch_map = resolve_branch_map(tree, prefix)
    pdg = tree[branch_map["pdgId"]].array(library="ak")
    status = None
    if branch_map["status"] in tree.keys():
        status = tree[branch_map["status"]].array(library="ak")

    if all(branch_map[key] in tree.keys() for key in ("px", "py", "pz", "energy")):
        px = tree[branch_map["px"]].array(library="ak")
        py = tree[branch_map["py"]].array(library="ak")
        pz = tree[branch_map["pz"]].array(library="ak")
        energy = tree[branch_map["energy"]].array(library="ak")
    elif all(branch_map[key] in tree.keys() for key in ("pt", "eta", "phi", "mass")):
        pt = tree[branch_map["pt"]].array(library="ak")
        eta = tree[branch_map["eta"]].array(library="ak")
        phi = tree[branch_map["phi"]].array(library="ak")
        mass = tree[branch_map["mass"]].array(library="ak")
        px = pt * np.cos(phi)
        py = pt * np.sin(phi)
        pz = pt * np.sinh(eta)
        energy = np.sqrt((pt * np.cosh(eta)) ** 2 + mass**2)
    else:
        raise KeyError("No compatible momentum branches found for GenParticles.")

    return {
        "pdgId": pdg,
        "status": status,
        "px": px,
        "py": py,
        "pz": pz,
        "energy": energy,
    }


def select_group(arrays: dict[str, ak.Array], group: PdgGroup) -> dict[str, ak.Array]:
    pdg = arrays["pdgId"]
    if group.abs:
        pdg = abs(pdg)
    mask = ak.zeros_like(pdg, dtype=bool)
    for pid in group.ids:
        if group.abs:
            mask = mask | (pdg == abs(pid))
        else:
            mask = mask | (pdg == pid)
    return {
        "px": arrays["px"][mask],
        "py": arrays["py"][mask],
        "pz": arrays["pz"][mask],
        "energy": arrays["energy"][mask],
    }


def combinations_masses(arrays: dict[str, ak.Array], groups: tuple[PdgGroup, ...]) -> ak.Array:
    grouped = []
    for group in groups:
        selected = select_group(arrays, group)
        grouped.append(ak.zip(selected))

    cartesian_inputs = {f"g{idx}": group for idx, group in enumerate(grouped)}
    combos = ak.cartesian(cartesian_inputs, axis=1)

    sum_px = ak.zeros_like(combos["g0"].px)
    sum_py = ak.zeros_like(combos["g0"].py)
    sum_pz = ak.zeros_like(combos["g0"].pz)
    sum_energy = ak.zeros_like(combos["g0"].energy)

    for idx in range(len(grouped)):
        sum_px = sum_px + combos[f"g{idx}"].px
        sum_py = sum_py + combos[f"g{idx}"].py
        sum_pz = sum_pz + combos[f"g{idx}"].pz
        sum_energy = sum_energy + combos[f"g{idx}"].energy

    mass_sq = sum_energy**2 - (sum_px**2 + sum_py**2 + sum_pz**2)
    mass_sq = ak.where(mass_sq > 0, mass_sq, 0)
    return ak.sqrt(mass_sq)


def histogram_masses(masses: ak.Array, binning: Binning) -> tuple[np.ndarray, np.ndarray]:
    flat = ak.to_numpy(ak.flatten(masses, axis=None))
    counts, edges = np.histogram(flat, bins=binning.bins, range=(binning.min, binning.max))
    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, counts


def write_csv(path: Path, centers: np.ndarray, counts: np.ndarray) -> None:
    stat_err = np.sqrt(np.maximum(counts, 1))
    data = np.column_stack([centers, counts, stat_err])
    header = "mass_GeV,counts,stat_err"
    np.savetxt(path, data, fmt=["%.6f", "%d", "%.6f"], delimiter=",", header=header, comments="")


def parse_channels(spec: dict[str, Any]) -> list[ChannelSpec]:
    channels = []
    for channel in spec.get("channels", []):
        groups = []
        for group in channel.get("pdg_groups", []):
            ids = tuple(int(pid) for pid in group.get("ids", []))
            groups.append(PdgGroup(ids=ids, abs=bool(group.get("abs", False))))
        binning = channel.get("binning", {})
        channels.append(
            ChannelSpec(
                name=str(channel.get("name")),
                label=str(channel.get("label")),
                pdg_groups=tuple(groups),
                combination=str(channel.get("combination", "pair")),
                binning=Binning(
                    bins=int(binning.get("bins", 60)),
                    min=float(binning.get("min", 0.0)),
                    max=float(binning.get("max", 10.0)),
                ),
                output_csv=str(channel.get("output_csv")),
            )
        )
    return channels


def run_analysis(input_root: Path, output_dir: Path, channels: list[ChannelSpec], tree_name: str = "Events") -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    tree = uproot.open(str(input_root))[tree_name]
    arrays = build_vectors(tree)

    for channel in channels:
        masses = combinations_masses(arrays, channel.pdg_groups)
        centers, counts = histogram_masses(masses, channel.binning)
        output_path = output_dir / channel.output_csv
        write_csv(output_path, centers, counts)


def main() -> None:
    parser = argparse.ArgumentParser(description="GEN-only mass spectra analyzer (SYNTHETIC)")
    parser.add_argument("--spec", required=True, help="Path to spec_cmssw.yaml")
    parser.add_argument("--input", required=True, help="Input GEN ROOT file")
    parser.add_argument("--output-dir", required=True, help="Output directory for CSVs")
    parser.add_argument("--tree", default="Events", help="Tree name (default: Events)")
    parser.add_argument("--run", action="store_true", help="Execute analysis (default is dry-run)")
    args = parser.parse_args()

    spec = load_spec(Path(args.spec))
    channels = parse_channels(spec)

    if not args.run:
        print("DRY_RUN: No analysis executed.")
        print(f"Input GEN file: {args.input}")
        print(f"Output dir: {args.output_dir}")
        print(f"Tree: {args.tree}")
        print(f"Channels: {[channel.name for channel in channels]}")
        print("Re-run with --run to execute.")
        return

    run_analysis(Path(args.input), Path(args.output_dir), channels, tree_name=args.tree)


if __name__ == "__main__":
    main()
