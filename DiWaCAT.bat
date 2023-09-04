@echo off

:: Ensure this batch file stops if any of the following commands fail
setlocal enabledelayedexpansion
echo Checking if Python is installed...

:: Check for Python
where python >nul 2>nul
if errorlevel 1 (
    echo Python not found! Please ensure it's installed and added to PATH.
    exit /b 1
)

echo Python found! Checking dependencies...

:: Check and install dependencies (assuming you have a requirements.txt file)
python -m pip install --user -r dependencies.txt >NUL
set PIP_ERRORLEVEL=%ERRORLEVEL%

if %PIP_ERRORLEVEL% neq 0 (
    echo Failed to install some dependencies. Continuing anyway...
)

echo Attempting to Run DiWaCAT
python DiWaCAT_GUI.py

pause