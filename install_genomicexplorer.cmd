@echo off
setlocal EnableExtensions

REM ============================================================
REM GenomicExplorer - Windows installer (micromamba)
REM Creates/updates conda-style env from environment.windows.yml
REM ============================================================

cd /d "%~dp0"

REM --- locate micromamba.exe (user-local / appdata / programdata) ---
set "MM=%USERPROFILE%\micromamba\micromamba.exe"
if not exist "%MM%" set "MM=%USERPROFILE%\AppData\Local\micromamba\micromamba.exe"
if not exist "%MM%" set "MM=%ProgramData%\micromamba\micromamba.exe"

if not exist "%MM%" (
  echo [GE][ERROR] micromamba.exe not found.
  echo Looked for:
  echo   %USERPROFILE%\micromamba\micromamba.exe
  echo   %USERPROFILE%\AppData\Local\micromamba\micromamba.exe
  echo   %ProgramData%\micromamba\micromamba.exe
  echo.
  echo Please install micromamba or adjust this script.
  pause
REM  exit /b 1
)

REM --- ensure "micromamba" is visible to child processes (QProcess etc.) ---
for %%I in ("%MM%") do set "MMDIR=%%~dpI"
set "PATH=%MMDIR%;%PATH%"

REM --- env name (change here if you want) ---
set "ENV_NAME=GenomicExplorer"

REM --- choose yml (windows-first) ---
set "YML=environment.windows.yml"
if not exist "%YML%" set "YML=environment.yml"

if not exist "%YML%" (
  echo [GE][ERROR] environment.windows.yml or environment.yml not found in project root.
  pause
REM  exit /b 1
)

echo [GE] micromamba: "%MM%"
echo [GE] env:       %ENV_NAME%
echo [GE] yml:       %YML%
echo.

echo [GE] Creating env (if exists, update)...
"%MM%" create -y -n "%ENV_NAME%" -f "%YML%"
pause
if errorlevel 1 (
  echo [GE] create failed (maybe already exists). Trying env update...
  "%MM%" env update -y -n "%ENV_NAME%" -f "%YML%"
  if errorlevel 1 (
    echo [GE][ERROR] Failed to create/update environment.
    pause
REM    exit /b 1
  )
)

echo.
echo [GE] Quick check: python & Rscript
"%MM%" run -n "%ENV_NAME%" python -V
if errorlevel 1 (
  echo [GE][ERROR] python check failed.
  pause
REM  exit /b 1
)

"%MM%" run -n "%ENV_NAME%" Rscript -e "cat('R OK\n'); sessionInfo()"
if errorlevel 1 (
  echo [GE][ERROR] Rscript check failed.
  pause
REM  exit /b 1
)

echo.
echo [GE] Install OK.
echo Next:
echo   1) (optional) install_optional_R.cmd  ^(if you ship it^)
echo   2) run_genomicexplorer.cmd
pause
endlocal
