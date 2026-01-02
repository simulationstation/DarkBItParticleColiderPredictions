#!/usr/bin/env bash
set -euo pipefail

FAMILY_ID="belle_zb"
CMSSW_RELEASE="${CMSSW_RELEASE:-CMSSW_14_0_13}"
EVENTS="${EVENTS:-50000}"
DRY_RUN="${DRY_RUN:-1}"
RUN_RANK1="${RUN_RANK1:-0}"

if [[ "${DRY_RUN}" == "1" ]]; then
  echo "DRY_RUN=1: no commands will be executed."
  echo "Set DRY_RUN=0 and RUN=1 to execute."
  exit 0
fi

if [[ "${RUN:-0}" != "1" ]]; then
  echo "RUN=1 is required to execute." >&2
  exit 1
fi

if [[ ! -d /cvmfs ]]; then
  echo "/cvmfs is required and must be mounted in the container." >&2
  exit 1
fi

WORKDIR="${WORKDIR:-/work}"
FAMILY_DIR="${WORKDIR}/EXOTICS_FACTORY/families/${FAMILY_ID}/cmssw_synth"
SPEC_PATH="${FAMILY_DIR}/spec_cmssw.yaml"
GEN_FRAGMENT_SRC="${FAMILY_DIR}/gen_fragment.py"
ANALYZER_SRC="${FAMILY_DIR}/analyze_gen.py"
OUTPUT_DIR="${WORKDIR}/EXOTICS_FACTORY/out/cmssw_synth/${FAMILY_ID}"

source /cvmfs/cms.cern.ch/cmsset_default.sh

scram project CMSSW "${CMSSW_RELEASE}"
cd "${CMSSW_RELEASE}/src"
eval "$(scram runtime -sh)"

mkdir -p Configuration/Generator/python
cp "${GEN_FRAGMENT_SRC}" "Configuration/Generator/python/${FAMILY_ID}_gen.py"

cmsDriver.py Configuration/Generator/python/${FAMILY_ID}_gen.py \
  --python_filename "${FAMILY_ID}_GEN.py" \
  --eventcontent GEN \
  --datatier GEN \
  --conditions auto:phase1_2024_realistic \
  --step GEN \
  --geometry DB:Extended \
  --era Run3 \
  --no_exec \
  -n "${EVENTS}"

cmsRun "${FAMILY_ID}_GEN.py"

python "${ANALYZER_SRC}" \
  --spec "${SPEC_PATH}" \
  --input "${FAMILY_ID}_GEN.root" \
  --output-dir "${OUTPUT_DIR}" \
  --run

if [[ "${RUN_RANK1}" == "1" ]]; then
  python "${WORKDIR}/EXOTICS_FACTORY/toolkit/rank1_harness.py" \
    --input-dir "${OUTPUT_DIR}" \
    --family "${FAMILY_ID}" \
    --backend "cmssw_synth" \
    --run
fi
