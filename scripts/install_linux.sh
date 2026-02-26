#!/usr/bin/env bash
# scripts/install_linux.sh
# GenomicExplorer: create/update micromamba env from environment.yml (and optional micromamba-lock)
set -euo pipefail

ENV_NAME="${ENV_NAME:-GenomicExplorer}"
YML_PATH="${YML_PATH:-environment.yml}"
LOCK_LINUX="${LOCK_LINUX:-micromamba-linux-64.lock}"

log() { echo "[GE][install] $*"; }

fail() {
  echo "[GE][install][ERROR] $*" >&2
  exit 1
}

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || fail "Command not found: $1"
}

# ---- locate project root (script can be run from anywhere) ----
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

cd "$ROOT_DIR"

log "Project root: $ROOT_DIR"
log "Target env:    $ENV_NAME"

# ---- micromamba availability ----
need_cmd micromamba

# ---- make micromamba usable inside non-interactive script ----
# shellcheck disable=SC1091
#source "$(micromamba info --base)/etc/profile.d/micromamba.sh"

# ---- show basic info ----
log "micromamba base:    $(micromamba info --base)"
log "micromamba version: $(micromamba --version)"

# ---- validate input files ----
[ -f "$YML_PATH" ] || fail "Missing $YML_PATH (expected at project root)."

# ---- optional: prefer micromamba-lock if present ----
USE_LOCK="0"
if [ -f "$LOCK_LINUX" ] && command -v micromamba-lock >/dev/null 2>&1; then
  USE_LOCK="1"
  log "Found $LOCK_LINUX and micromamba-lock. Will install from lock (reproducible)."
else
  log "Installing from $YML_PATH (no lock or micromamba-lock not installed)."
fi

# ---- create/update env ----
if micromamba env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  log "micromamba env already exists: $ENV_NAME"
  if [ "$USE_LOCK" = "1" ]; then
    log "Updating env from lock: $LOCK_LINUX"
    micromamba-lock install -n "$ENV_NAME" "$LOCK_LINUX"
  else
    log "Updating env from yml:  $YML_PATH"
    micromamba env update -n "$ENV_NAME" -f "$YML_PATH" --prune
  fi
else
  log "Creating micromamba env: $ENV_NAME"
  if [ "$USE_LOCK" = "1" ]; then
    micromamba-lock install -n "$ENV_NAME" "$LOCK_LINUX"
  else
    micromamba env create -n "$ENV_NAME" -f "$YML_PATH"
  fi
fi

# ---- quick sanity checks ----
log "Activating env..."
MM_BASE="$(micromamba info --base)"
MM_PROFILE="${MM_BASE}/etc/profile.d/micromamba.sh"
if [ -f "$MM_PROFILE" ]; then
  source "$MM_PROFILE"
else
  eval "$(micromamba shell hook --shell bash)"
fi
micromamba activate "$ENV_NAME"

log "Python: $(python -V)"
log "R:      $(R --version | head -n 1 || true)"

# Plotly image export depends on kaleido (Python)
python - <<'PY'
import sys
ok = True
try:
    import plotly  # noqa
except Exception as e:
    ok = False
    print("[CHECK] plotly import failed:", e, file=sys.stderr)
try:
    import kaleido  # noqa
except Exception as e:
    ok = False
    print("[CHECK] kaleido import failed:", e, file=sys.stderr)
print("[CHECK] python deps:", "OK" if ok else "NG")
sys.exit(0 if ok else 2)
PY

log "Done."
log "Next: bash scripts/run_linux.sh"
