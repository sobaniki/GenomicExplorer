#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="GenomicExplorer"

# micromamba activate をスクリプト内で使う定番
#source "$(micromamba info --base)/etc/profile.d/micromamba.sh"
MM_BASE="$(micromamba info --base)"
MM_PROFILE="${MM_BASE}/etc/profile.d/micromamba.sh"
if [ -f "$MM_PROFILE" ]; then
  source "$MM_PROFILE"
else
  eval "$(micromamba shell hook --shell bash)"
fi
micromamba activate "$ENV_NAME"

export R_LIBS_USER="$(cd "$(dirname "$0")/.." && pwd)/_r_libs"
# 必要な人だけ：QtWebEngineが不安定な環境で試す用（READMEに“任意”として記載）
# export QTWEBENGINE_CHROMIUM_FLAGS="--disable-gpu --disable-gpu-compositing"

cd GE
python -m gui.app
