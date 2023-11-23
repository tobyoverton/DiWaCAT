# DiWaCAT
Dielectric Wakefield Calculator and Tracker

[![DOI](https://zenodo.org/badge/686930732.svg)](https://zenodo.org/badge/latestdoi/686930732)

A Python and C++ based wakefield solver for relativistic electron bunches in dielectric lined waveguides (DLW). The files outlining each function of the DiWaCAT GUI is in the GUI_Source folder and the C++ implementation of the field calculations (planar using [doi.org/10.1103/PhysRevSTAB.16.051302] and circular using [doi.org/10.1103/PhysRevD.42.1819]) are given in Cpp_Source.

The python dependencies will be installed from dependecies.txt using pip. The following dependencies are also needed:
  * Python - minimum version 3.11 (untested for Python 3.8 - 3.11)
  * CMake - minimum version 3.12
  * C++ compiler - minimum compatability with C++11
DiWaCAT has been built on Windows and untested on Linux/macOS. As a minimum these platforms will require a change to Cpp_Source/CMakeLists.txt Line 26 to set the correct suffix for built Python library.

For a full outline of capabilities/tips see the DiWaCAT Manual document. Full documentation to replace this is coming soon.
