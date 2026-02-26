# scripts/run_windows.ps1（micromamba版）
$ErrorActionPreference = "Stop"

$ENV_NAME = $env:ENV_NAME; if ([string]::IsNullOrWhiteSpace($ENV_NAME)) { $ENV_NAME = "GenomicExplorer" }

$SCRIPT_DIR = Split-Path -Parent $MyInvocation.MyCommand.Path
$ROOT_DIR = Resolve-Path (Join-Path $SCRIPT_DIR "..")
Set-Location $ROOT_DIR

# project-local R libs (任意)
$env:R_LIBS_USER = (Join-Path $ROOT_DIR "_r_libs")
New-Item -ItemType Directory -Force -Path $env:R_LIBS_USER | Out-Null

# WebEngineは安定優先でOFF（必要なら1に）
if ([string]::IsNullOrWhiteSpace($env:GE_ENABLE_WEBENGINE)) { $env:GE_ENABLE_WEBENGINE = "0" }

# もしportable micromambaで root prefix が必要ならここで指定（必要な場合だけ）
# $env:MAMBA_ROOT_PREFIX = "C:\path\to\micromamba-root"

Set-Location (Join-Path $ROOT_DIR "GE_v0223b")
micromamba run -n $ENV_NAME -- python -m gui.app
