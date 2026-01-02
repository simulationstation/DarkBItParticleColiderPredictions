"""Master launcher for EXOTICS_FACTORY pipelines (dry-run by default)."""
from __future__ import annotations

import argparse
import importlib
import json
from pathlib import Path
from typing import Any

import yaml

from EXOTICS_FACTORY.toolkit.common.logging import configure_logger


LOGGER = configure_logger(__name__)
SUPPORTED_BACKENDS = ("pdf", "hepdata", "workspace", "cmssw_synth")


def load_registry(path: str | Path) -> dict[str, Any]:
    return yaml.safe_load(Path(path).read_text(encoding="utf-8"))


def find_family(registry: dict[str, Any], family_id: str) -> dict[str, Any]:
    for entry in registry.get("families", []):
        if entry.get("id") == family_id:
            return entry
    raise KeyError(f"Family '{family_id}' not found in registry")


def load_spec(family_id: str) -> dict[str, Any]:
    spec_path = Path("EXOTICS_FACTORY/families") / family_id / "spec.yaml"
    return yaml.safe_load(spec_path.read_text(encoding="utf-8"))


def run_pipeline(family_id: str, mode: str) -> None:
    module_path = f"EXOTICS_FACTORY.families.{family_id}.pipeline"
    pipeline = importlib.import_module(module_path)
    spec = load_spec(family_id)
    pipeline.run_pipeline(spec, mode=mode)


def load_cmssw_spec(family_id: str) -> dict[str, Any]:
    spec_path = Path("EXOTICS_FACTORY/families") / family_id / "cmssw_synth" / "spec_cmssw.yaml"
    return yaml.safe_load(spec_path.read_text(encoding="utf-8"))


def verify_cmssw_assets(family_id: str, spec: dict[str, Any]) -> list[Path]:
    base = Path("EXOTICS_FACTORY/families") / family_id / "cmssw_synth"
    required = [
        base / "spec_cmssw.yaml",
        base / "gen_fragment.py",
        base / "analyze_gen.py",
        base / "RUNBOOK_CMSSW_SYNTH.md",
        base / "run_docker_genonly.sh",
    ]
    cmssw = spec.get("cmssw", {})
    for key in ("generator_fragment", "analyzer_module", "analyzer_script"):
        path = cmssw.get(key)
        if path:
            required.append(Path(path))
    missing = [path for path in required if not path.exists()]
    if missing:
        missing_str = ", ".join(str(path) for path in missing)
        raise FileNotFoundError(f"Missing CMSSW synth assets: {missing_str}")
    return required


def run_cmssw_synth(family_id: str, mode: str) -> None:
    spec = load_cmssw_spec(family_id)
    steps = spec.get("runtime", {}).get("steps", [])
    if mode != "RUN":
        LOGGER.info("DRY_RUN: CMSSW synth will not execute.")
        LOGGER.info("Planned steps: %s", json.dumps(steps, indent=2))
        verify_cmssw_assets(family_id, spec)
        LOGGER.info("Verified CMSSW synth assets for %s", family_id)
        return
    raise RuntimeError("Execution disabled in factory scaffolding. Use DRY_RUN only.")


def main() -> None:
    parser = argparse.ArgumentParser(description="EXOTICS_FACTORY pipeline launcher")
    parser.add_argument("--family", required=True, help="Family id from registry")
    parser.add_argument(
        "--backend",
        choices=SUPPORTED_BACKENDS,
        help="Backend to use (default: first preferred backend, or cmssw_synth)",
    )
    parser.add_argument("--list-backends", action="store_true", help="List backend options and exit")
    parser.add_argument("--run", action="store_true", help="Execute pipeline steps")
    args = parser.parse_args()

    registry = load_registry("EXOTICS_FACTORY/registry/families.yaml")
    entry = find_family(registry, args.family)
    LOGGER.info("Selected family: %s", json.dumps(entry, indent=2))

    preferred = entry.get("preferred_backends", [])
    available_backends = sorted(set(preferred) | {"cmssw_synth"})
    LOGGER.info("Available backends: %s", available_backends)

    if args.list_backends:
        print("\n".join(available_backends))
        return

    backend = args.backend or (preferred[0] if preferred else "cmssw_synth")
    mode = "RUN" if args.run else "DRY_RUN"

    if backend == "cmssw_synth":
        run_cmssw_synth(args.family, mode)
        return

    if backend not in preferred:
        LOGGER.warning("Backend '%s' not listed in preferred_backends for %s", backend, args.family)

    run_pipeline(args.family, mode)


if __name__ == "__main__":
    main()
