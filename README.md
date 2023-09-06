# DiWaCAT
Dielectric Wakefield Calculator and Tracker

[![DOI](https://zenodo.org/badge/686930732.svg)](https://zenodo.org/badge/latestdoi/686930732)

A Python and C++ based wakefield solver for relativistic electron bunches in dielectric lined waveguides (DLW).

For windows/linux, the FieldSolver_Executable folder must also be downloaded for all parts of DiWaCAT_GUI.exe to run. This folder contains an executable file and dependencies for the C++ based field solver. The full calculator and tracker is given by the DiWaCAT_GUI file. The .py file is included for running on iOS machines. This requires pyqt and pyqtgraph packages to be installed however the field calculator part of the code will not run.

A folder of C++ files are included alongside the field solver executable. These solutions can be used to independently calculate fields in planar and circular DLWs, using the functions outlined in FieldFunctions.h. Note using these files directly requires correctly linking to the HDF5 C++ library files.

For a full outline of capabilities/tips see the DiWaCAT Manual document. Full documentation to replace this is coming soon.
