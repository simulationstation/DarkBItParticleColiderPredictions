# CMSSW Synthetic Backend (GEN-only)

This directory hosts **synthetic** CMSSW GEN-only templates that mimic the
rank-1 harness interface. These are **not** physics models of the experiments;
they exist to provide a consistent, reproducible GEN-only proxy for
identifiability studies.

## What lives here

- `analyzers/`: reusable GEN-only analyzer that writes CSV spectra.
- `generators/`: CMSSW generator fragment templates (SYNTHETIC).
- `docker_runner_templates/`: run scripts that create CMSSW areas, execute GEN,
  and run analyzers (DRY_RUN by default).
- `runbooks/`: general runbook patterns and notes.
- `family_bindings/`: index of family -> cmssw_synth spec bindings.
- `CMSSW_SYNTH_AUDIT.md`: feasibility classification for each family.

## Output format

All synthetic outputs must be clearly labeled `SYNTHETIC` and use the
rank-1 harness CSV schema:

```
mass_GeV,counts,stat_err
```

`stat_err` is computed as `sqrt(max(counts, 1))` (Poisson).

## Usage pattern (later)

1. Copy the family generator fragment into a CMSSW release area.
2. Run `cmsDriver.py` to generate a GEN-only ROOT file.
3. Run the analyzer script to produce per-channel CSVs.
4. Feed the CSVs into the rank-1 harness pipeline.

> **Important:** All scripts default to DRY_RUN and refuse to execute unless
> explicitly invoked with a `--run` flag or `RUN=1`.
