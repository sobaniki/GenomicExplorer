@echo off
setlocal EnableExtensions

REM ============================================================
REM GenomicExplorer - Windows launcher (micromamba)
REM Starts GUI inside env and sets stable R_LIBS_USER
REM ============================================================

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

REM --- ensure "micromamba" is visible to child processes (plugin runners) ---
for %%I in ("%MM%") do set "MMDIR=%%~dpI"
set "PATH=%MMDIR%;%PATH%"

REM --- env name (must match installer) ---
set "ENV_NAME=GenomicExplorer"

REM --- IMPORTANT: project-local R library (so GUI and installers share the same lib) ---
set "R_LIBS_USER=%~dp0_r_libs"
if not exist "%R_LIBS_USER%" mkdir "%R_LIBS_USER%"

REM --- Stable default: embedded webengine OFF (set to 1 if you really need it) ---
if "%GE_ENABLE_WEBENGINE%"=="" set "GE_ENABLE_WEBENGINE=1"

REM --- start GUI inside env ---
pushd "GE"
"%MM%" run -n "%ENV_NAME%" python -m gui.app
set "RC=%ERRORLEVEL%"
popd

if not "%RC%"=="0" (
  echo [GE][ERROR] GUI exited with code %RC%
)
pause
endlocal
