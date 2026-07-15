#!/usr/bin/env bash
set -euo pipefail

ENV_FILE="envs/dnaseq.yml"

if command -v micromamba >/dev/null 2>&1; then
  micromamba create -f "${ENV_FILE}"
elif command -v conda >/dev/null 2>&1; then
  conda env create -f "${ENV_FILE}"
else
  cat >&2 <<'MSG'
Could not find micromamba or conda.
Please install Miniforge/Mambaforge, conda, or micromamba, then rerun this script.
MSG
  exit 127
fi
