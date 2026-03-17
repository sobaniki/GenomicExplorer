#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${ENV_NAME:-GenomicExplorer}"
PROFILE="${1:-full}"   # core / full

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

# micromamba activate
#source "$(micromamba info --base)/etc/profile.d/micromamba.sh"
MM_BASE="$(micromamba info --base)"
MM_PROFILE="${MM_BASE}/etc/profile.d/micromamba.sh"
if [ -f "$MM_PROFILE" ]; then
  source "$MM_PROFILE"
else
  eval "$(micromamba shell hook --shell bash)"
fi
micromamba activate "$ENV_NAME"

# project-local R library (so GUI runners can see the same libs)
#export R_LIBS_USER="${ROOT_DIR}/_r_libs"
#mkdir -p "$R_LIBS_USER"

echo "[GE][install_optional_R] ENV=$ENV_NAME"
#echo "[GE][install_optional_R] R_LIBS_USER=$R_LIBS_USER"
echo "[GE][install_optional_R] PROFILE=$PROFILE"

Rscript "${SCRIPT_DIR}/install_optional_R.R" --profile "$PROFILE"

echo "[GE][install_optional_R] Done."
