# scripts/install_windows.ps1
$ErrorActionPreference = "Stop"

$ENV_NAME  = $env:ENV_NAME;  if ([string]::IsNullOrWhiteSpace($ENV_NAME))  { $ENV_NAME  = "GenomicExplorer" }
$YML_PATH  = $env:YML_PATH;  if ([string]::IsNullOrWhiteSpace($YML_PATH))  { $YML_PATH  = "environment.windows.yml" }
$LOCK_FILE = $env:LOCK_FILE; if ([string]::IsNullOrWhiteSpace($LOCK_FILE)) { $LOCK_FILE = "micromamba-win-64.lock" }

function Log($msg) { Write-Host "[GE][install] $msg" }

# project root (script can be called anywhere)
$SCRIPT_DIR = Split-Path -Parent $MyInvocation.MyCommand.Path
$ROOT_DIR = Resolve-Path (Join-Path $SCRIPT_DIR "..")
Set-Location $ROOT_DIR

Log "Project root: $ROOT_DIR"
Log "Target env:   $ENV_NAME"

# micromamba hook for PowerShell
$micromamba = Get-Command micromamba -ErrorAction Stop
# & micromamba "shell.powershell" "hook" | Out-String | Invoke-Expression

if (!(Test-Path $YML_PATH) -and !(Test-Path $LOCK_FILE)) {
  throw "Missing $YML_PATH and $LOCK_FILE in project root."
}

# Prefer micromamba-lock if available + lockfile exists
$useLock = $false
if ((Test-Path $LOCK_FILE) -and (Get-Command micromamba-lock -ErrorAction SilentlyContinue)) {
  $useLock = $true
  Log "Using lockfile: $LOCK_FILE"
} else {
  Log "Using env yml:  $YML_PATH"
}

# create/update env
$envExists = micromamba env list | Select-String -Pattern "^\s*$ENV_NAME\s"
if ($envExists) {
  Log "Env exists: $ENV_NAME"
  if ($useLock) {
    micromamba-lock install -n $ENV_NAME $LOCK_FILE
  } else {
    micromamba env update -n $ENV_NAME -f $YML_PATH --prune
  }
} else {
  Log "Creating env: $ENV_NAME"
  if ($useLock) {
    micromamba-lock install -n $ENV_NAME $LOCK_FILE
  } else {
    micromamba env create -n $ENV_NAME -f $YML_PATH
  }
}

Log "Activating env..."
MM_BASE="$(micromamba info --base)"
MM_PROFILE="${MM_BASE}/etc/profile.d/micromamba.sh"
if [ -f "$MM_PROFILE" ]; then
  source "$MM_PROFILE"
else
  eval "$(micromamba shell hook --shell bash)"
fi
micromamba activate $ENV_NAME

python -V
Rscript --version

Log "Done. Next: scripts\run_windows.ps1"
