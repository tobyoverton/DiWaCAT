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
	pause
)

echo Attempting to Build Cpp Library
cd Cpp_Source
mkdir build
cd build
setlocal enabledelayedexpansion
:: Find where pybind11 is installed and what architecture to build the python library to
for /f %%i in ('python -c "import pybind11; print(pybind11.get_cmake_dir())"') do set PYBIND11_PATH=%%i
for /f %%i in ('python -c "import struct; print(8 * struct.calcsize('P'))"') do set PYTHON_ARCH=%%i
if %PYTHON_ARCH%==64 (
    set CMAKE_ARCH=x64
) else (
    set CMAKE_ARCH=Win32
)

cmake -DCMAKE_BUILD_TYPE=Release -Dpybind11_DIR=!PYBIND11_PATH! -A !CMAKE_ARCH! ..
if errorlevel 1 (
    echo CMake configuration failed!
	pause
    exit /b 1
)
echo Building with CMake...
cmake --build . --config Release
if errorlevel 1 (
    echo CMake build failed!
	pause
    exit /b 1
)


cd ../..
echo Attempting to Run DiWaCAT
python DiWaCAT_GUI.py

pause