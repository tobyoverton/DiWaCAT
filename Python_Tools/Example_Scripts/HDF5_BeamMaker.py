#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example/Template Beam Creation Script

List the beam parameters you want the beam to have at the waist
    •Generate the beam and allow for beam manipulation
    •Mesh the beam for use with DiWaCAT
"""

import os
import os, sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const

#Make sure the directory is set so the folder Python_Tools is visible
from Python_Tools.Modules import beam_tools as dbt

# ---------------------- 1) create the beam --------------------------- #

OUTPUT_FILENAME = "./../DWA_FEBE/HV_BBU_Suppression/200fs_BeamIn.h5"

# How many macros?
N_MACROS = 50000

# Make empty dict to hold the parameters we want the constructor class to call
bp = {}

bp["pz_MeV"] = 250 # MeV/c
bp["eps_x_N"] = 1.5 * const.micro   # m rad
bp["eps_y_N"] = 1.5 * const.micro  # m rad
bp["sig_x_0"] = 50 * const.micro  # m
bp["sig_y_0"] = 50 * const.micro  # m
SIGMA_T = 200e-15 * 2.998e14 # s converted to micro
bp["sig_z_0"] = SIGMA_T * const.micro # m (I know its messy to convert to micro and then convert back to m - allows for either s or um sig_long)
bp["plat_rise"] =  0.0* SIGMA_T * const.speed_of_light #For plateau type beam (rise and fall time)
bp["sig_pz_0"] = 0.2  # MeV/c
bp["lCorr_Fac"] = 0.95  # An arbitrary number, but want to convert to chirp at some point
#If you want beam offset you can choose here (best to set these in DiWaCAT though)
#May be useful to change here if for example you want a random offset value
bp["x_0"] = 0
bp["y_0"] = 0
bp["xp_0"] = 0
bp["yp_0"] = 0
bp["z_0"] = 0
TOTAL_CHARGE = 250e-12 #Bunch charge in C
bp["charge_per_macro"] = TOTAL_CHARGE / N_MACROS #Charge per macro

#Profile type either 'Gaussian', 'SkewGaussian' 'Uniform', 'Plataeu', or 'DoubleGauss'
bp["LongitudinalProfile"] = 'Gaussian'
#bp["sig_z_2"] = 200e-15 * 2.998e8;
#bp["offset"] = 3.5*200e-15 * 2.998e8;
#bp["rel_amp"] = 1.25
#bp["skew"]=-6;

# initalise the beam instance
myBeam = dbt.GaussianBeam()

# apply the parameters we want
myBeam.setBeamParameters(bp)

# Calculate the optics functions at the waist
myBeam.calcWaistOpticFunctions()


#Can change a parameter like this:
#myBeam.beam["x_0"] = 200

#Can view or manipulate the entire array like:
#myBeam._beam['xp']

# Key parameters for the definition of a 6D Gaussian beam
#print(myBeam.meanArray)
#print(myBeam.covarianceMatrix)

# Generate the beam by specifying a number of macros
myBeam.generate6DMacros(N_MACROS)


#If you wanted to have the beam not at the waist (e.g. want the waist 1 m downstream)
#myBeam.driftBeam(-0.1)
#myBeam._beam['z'] += 0.1

#Example of how to print a variable if you wanted:
#print("Beam gamma is", myBeam.beam["gamma"])
#print("Beam sigma_x at waist is",myBeam.Sx/const.micro," um")
#print("beam eps_x_N is",myBeam.normalized_horizontal_emittance/const.micro," mm mrad")

#---------------------------------------------------------------------------
#---------------------- 2) Bin the data ------------------------------
#---------------------------------------------------------------------------


#Setup parameters for the binning mesh

cps = 3 # 3 cells per sigma for transverse
cps_long = 4 # 4 we'll use more cells longitudinally

Lx = 4.5 * myBeam.Sx # define over +/- 2.25 sigma - this works nicely for our Gaussian beam but we might want to be smarter here
Dx = myBeam.Sx/cps

Ly = 4.5 * myBeam.Sy
Dy = myBeam.Sy/(cps)

Lz = 6 * myBeam.Sz
Dz = myBeam.Sz/cps_long


# ------------------Perform the binning---------------------- #

#Initalise instance of beamBinner with the myBeam object
binnedBeam = dbt.beamBinner(myBeam)

#and bin
binnedBeam.binTheBeam(Lx=Lx,Dx=Dx,Ly=Ly,Dy=Dy,Lz=Lz,Dz=Dz)

# Test hist output
print("Shape of binned beam array is",binnedBeam.beamHist.shape)

#how much charge was captured on the mesh?
print("Fraction of total charge captured on this mesh is",binnedBeam.nMacroHist/myBeam.nMacros)


#---------------------------------------------------------------------------
#----------- 3) Smooth the data and compare smoothing methods -------------#
#---------------------------------------------------------------------------


binnedBeam.smoothGaussian(sigma=0.75)
print("The beam was smoothed with",binnedBeam.smoothMethod,"with kwargs",binnedBeam.smooth_kwargs)

binnedBeam.smoothMedian(size=3)
print("The beam was smoothed with",binnedBeam.smoothMethod,"with kwargs",binnedBeam.smooth_kwargs)

binnedBeam.smoothUniform(size=3)
print("The beam was smoothed with",binnedBeam.smoothMethod,"with kwargs",binnedBeam.smooth_kwargs)


#---------------------------------------------------------------------------
#---------------------- 4) Output the file ------------------------------
#---------------------------------------------------------------------------

#Now need to make a h5 file for output,including macros (optional)
print("Writing to output file: ")
outputFileWritten = binnedBeam.write_hd5f_mesh_file(OUTPUT_FILENAME,outputMacros=True, overwrite_existing_file=True,includeSmoothed=True,includeUnSmoothed=False)

if outputFileWritten:
    print("Output file written")

else:
    print("Output file not written")
