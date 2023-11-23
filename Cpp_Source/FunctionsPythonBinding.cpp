#if defined(_MSC_VER) || defined(__MINGW32__)
#define strdup _strdup
#endif

#define _USE_MATH_DEFINES
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "PlanarGreensFunctions.h"
#include "CircularGreenFunction.h"

#ifndef strdup
char* strdup(const char* s) {
    char* d = new char[strlen(s) + 1];
    if (d == NULL) return NULL;
    strcpy(d, s);
    return d;
}
#endif

namespace py = pybind11;

PYBIND11_MODULE(DiWaCAT_library, m) {
    m.def("EigenvaluesCalculator", &EigenvaluesCalculator, "Eigenvalues Calculator");
    m.def("CalcWakeElement", &CalcWakeElement, "Wake Potential Calculator");
    m.def("TotalForceMeshHDF5", &TotalForceMeshHDF5, "Force Calculator Horizontal");
    m.def("TotalForceMeshHDF5VerticalPlate", &TotalForceMeshHDF5VerticalPlate, "Force Calculator Vertical");
	m.def("FindModes", &FindModes, "Find Circular Modes");
	m.def("ModeConvergence", &ModeConvergence, "Circular Mode Convergence");
	m.def("CalcWakeCyl", &CalcWakeCyl, "Circular Wake Potential");
	m.def("TotalForceMeshCircular", &TotalForceMeshCircular, "Wakefield calculation for circular DLW");
	
}