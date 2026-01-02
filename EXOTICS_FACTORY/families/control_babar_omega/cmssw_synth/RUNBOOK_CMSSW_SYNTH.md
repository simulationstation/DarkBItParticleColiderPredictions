# CMSSW SYNTH Runbook: control_babar_omega

> **SYNTHETIC ONLY**: This GEN-only pipeline is a toy proxy for rank-1
> identifiability tests. It is **not** a physics model of the experimental
> dataset.

## Preconditions

- Docker with `/cvmfs` mounted inside the container.
- CMSSW release available via CVMFS.
- Repo mounted at a known path (e.g. `/work`).

## Files

- Spec: `EXOTICS_FACTORY/families/control_babar_omega/cmssw_synth/spec_cmssw.yaml`
- Generator fragment: `EXOTICS_FACTORY/families/control_babar_omega/cmssw_synth/gen_fragment.py`
- Analyzer: `EXOTICS_FACTORY/families/control_babar_omega/cmssw_synth/analyze_gen.py`
- Docker runner: `EXOTICS_FACTORY/families/control_babar_omega/cmssw_synth/run_docker_genonly.sh`

## GEN-only workflow (template)

1. Create CMSSW area in the container.
2. Copy the generator fragment into `CMSSW_*/src/Configuration/Generator/python/`.
3. Run `cmsDriver.py` to build a GEN-only configuration.
4. Run `cmsRun` to produce `*_GEN.root`.
5. Run the analyzer to emit CSVs.

> All scripts default to **DRY_RUN**. Set `RUN=1` and `DRY_RUN=0` to execute.

## Outputs

Outputs are written to:

```
EXOTICS_FACTORY/out/cmssw_synth/control_babar_omega
```

Each CSV is labeled with `SYNTHETIC` and uses the rank-1 harness format:

```
mass_GeV,counts,stat_err
```

## Rank-1 harness

Once CSVs exist, feed them into the harness as usual:

```
python EXOTICS_FACTORY/toolkit/rank1_harness.py \
  --input-dir EXOTICS_FACTORY/out/cmssw_synth/control_babar_omega \
  --family control_babar_omega \
  --backend cmssw_synth \
  --run
```

## Disclaimer

This pipeline is a **synthetic identifiability test only** and should not be
used to draw physics conclusions about the experimental spectra.
