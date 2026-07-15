#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_NAME="${DNASEQ_ENV_NAME:-koi-dnaseq}"

export XDG_CACHE_HOME="${PROJECT_ROOT}/.cache"

if [[ -x "${PROJECT_ROOT}/.local/bin/micromamba" ]]; then
  exec "${PROJECT_ROOT}/.local/bin/micromamba" run \
    -r "${PROJECT_ROOT}/.micromamba" \
    -n "${ENV_NAME}" \
    "$@"
elif command -v micromamba >/dev/null 2>&1; then
  exec micromamba run -n "${ENV_NAME}" "$@"
elif command -v conda >/dev/null 2>&1; then
  exec conda run -n "${ENV_NAME}" "$@"
else
  cat >&2 <<MSG
Could not find micromamba or conda.

Create the environment first with one of:
  conda env create -f envs/dnaseq.yml
  micromamba create -f envs/dnaseq.yml

Then rerun this command, or activate the environment manually and run the tool directly.
MSG
  exit 127
fi
