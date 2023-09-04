#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example/Template Mesh Creation Script

Read in a beam from .h5 file and mesh using the parameters given

"""

import os
import sys
import scipy.constants as const

#Make sure the directory is set so the folder Python_Tools is visible
from Python_Tools.Modules import beam_tools as dbt, read_beam_file as rbf

# ---------------- Read In the File and Create Beam --------- #

# N.B. I don't think you can set the output file to overwrite the original
# Input beam file isn't closed in the script so wouldn't allow it
# You can overwrite other files
INPUT_FILE =  "./../Streaker_Work/Profile Reconstruction/Test_Streaks/Beams_Out/First_DLW/DoubleGausBeam.h5"
OUTPUT_FILE =  "./../Streaker_Work/Profile Reconstruction/Test_Streaks/Beams_Out/First_Meshed/DoubleGausBeam.h5"

# open a beam file and get it
myBeam = dbt.BeamFromFile()

#Reading a beam from DiWaCAT Field Solver:
myBeam.read_WakeCode_beam_file(INPUT_FILE)

#If you want to drift the beam do it here
#myBeam.driftBeam(1)


# -----------------  Set up the Mesh Parameters ------------ #

#This is using the RMS size of the beam to find mesh values
#If you know the beam is very not-Gaussian it's best to define the range yourself

#Setup parameters for the binning mesh
cps = 3 #cells per sigma - transverse
cpsLong = 3.5 #cells per sigma - longitudinal

#Total length of the mesh is L and length of each mesh cell is D
Lx = 4.5 * myBeam.Sx # define over +/- 3 sigma - this works nicely for our Gaussian beam but we might want to be smarter here
Dx = Lx/(4.5*cps)

Ly = 4.5 * myBeam.Sy
Dy = Ly/(4.5*cps)

Lz = 6 * myBeam.Sz
Dz = Lz/(6*cpsLong)

#e.g. if you knew Gaussian assumptions are wrong vertically and wanted a mesh across 1000um with each cell 50um
#This would make a mesh with 20 cells across
#Ly = 1000e-6;
#Dy = 50e-6


# -------------------  Perform the binning---------------------#

#Initalise instance of beamBinner with the myBeam object
binnedBeam = dbt.beamBinner(myBeam)

#and bin
binnedBeam.binTheBeam(Lx=Lx,Dx=Dx,Ly=Ly,Dy=Dy,Lz=Lz,Dz=Dz)
## Test hist full_kick_output
print("Shape of binned beam array is",binnedBeam.beamHist.shape)

#how much charge was captured on the mesh?
print("Fraction of total charge captured on this mesh is",binnedBeam.nMacroHist/myBeam.nMacros)

# ----------------- Save the Meshed Beam ------------------#

#Now need to make a h5 file for full_kick_output,including macros
print("Writing to full_kick_output file > ",OUTPUT_FILE,"at",os.getcwd())
outputFileWritten = binnedBeam.write_hd5f_mesh_file(OUTPUT_FILE,outputMacros=True, overwrite_existing_file=True,includeSmoothed=False,includeUnSmoothed=True)
if outputFileWritten:
    print("Output file written")
else:
    print("Output file not written")
