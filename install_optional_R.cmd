@echo off
setlocal EnableExtensions
cd /d "%~dp0"

REM --- locate micromamba.exe ---
  set "MM=%USERPROFILE%\micromamba\micromamba.exe"
if not exist "%MM%" set "MM=%USERPROFILE%\AppData\Local\micromamba\micromamba.exe"
if not exist "%MM%" set "MM=%ProgramData%\micromamba\micromamba.exe"

if not exist "%MM%" (
  echo [GE][ERROR] micromamba.exe not found.
  pause
  exit /b 1
)

REM --- ensure micromamba is visible to child processes ---
  for %%I in ("%MM%") do set "MMDIR=%%~dpI"
set "PATH=%MMDIR%;%PATH%"

set "ENV_NAME=GenomicExplorer"

REM --- IMPORTANT: match GUI library location ---
set "R_LIBS_USER=%~dp0_r_libs"
if not exist "%R_LIBS_USER%" mkdir "%R_LIBS_USER%"

REM --- run your installer R script (adjust path if needed) ---
  "%MM%" run -n "%ENV_NAME%" Rscript "%~dp0scripts\install_optional_R_windows.R" --profile full

pause
endlocal
