# -*- coding: utf-8 -*-
"""
Example script to read in a DiWaCAT beam/field and manipulate the class object
"""

import os, time, csv, sys, subprocess

import numpy as np
import math
import matplotlib.pyplot as plt
import scipy
import scipy.interpolate, scipy.optimize
import scipy.constants as const
import time
from Python_Tools.Modules import diwacat_tools as DWA
from Python_Tools.Modules import beam_tools as dbt
from fastkde import fastKDE



# ---------------- 1. Reading in a file and transporting through the DLW ----------------------- #

DWA_Beam = DWA.DiWaCATOutput()
#N.B. this example is for a 600fs long Gaussian bunch
#DWA_Beam.ReadFromFile(".\\DWA_Code\\DiWaCAT_Executable\\DefaultBeamName.h5")
DWA_Beam.ReadFromFile(".\\..\\..\\DWA_Experiments\\BA1_CircularSims\\Offset\\400fs\\Offset_400.h5")

M = np.hypot(DWA_Beam._FieldPoints['Fx'], DWA_Beam._FieldPoints['Fy'])
plt.quiver(DWA_Beam._FieldPoints['x']*1e6, DWA_Beam._FieldPoints['y']*1e6, DWA_Beam._FieldPoints['Fx'], DWA_Beam._FieldPoints['Fy'], width = 0.001, headwidth = 1)
plt.xlabel('x [um]')
plt.ylabel('y [um]')
plt.show()
'''
#Print of some characteristics
print("Dielectric Gap = ", DWA_Beam._DielectricParameter['a'] * 1e3, " mm")
print("Structure Length = ", DWA_Beam._DielectricParameter['L'], " m")

#Change the dielectric gap
DWA_Beam.SetDielectricLength(0.2);

#Accessing beam qualities
#.beam is exactly the same as a GeneralBeam in beam_tools so can use any of those functions
print("Normalised Horizontal Emittance = ", DWA_Beam.beam.normalized_horizontal_emittance * 1/const.micro, " mm mrad")

#Transport through the DLW
#If you give remove_losses = True there will also be a check for beam losses at each stage
DWA_Beam.TransportThroughDLW(nSteps = 3)
DWA_Beam.ClearLostParticles()

#See how the emittance has changed
print("Normalised Horizontal Emittance = ", DWA_Beam.beam.normalized_horizontal_emittance * 1/const.micro, " mm mrad")


# ------------------ 2. Using a different beam but fields from DiWaCat ------------------------- #

#If using a different sigmat to the original remember the beam field starts at 2.5*sigmat of the original
#z_0 would need to be used to make sure the beams start in the same place
#In this case 2.5sigmat = 2.5*400e-15*3e8 m

#The steps to make the beam are in HDF5_BeamMaker
N_MACROS = 50000
bp = {}
bp["pz_MeV"] = 250 # MeV/c
bp["eps_x_N"] = 2.5 * const.micro   # m rad
bp["eps_y_N"] = 2.5 * const.micro  # m rad
bp["sig_x_0"] = 20 * const.micro  # m
bp["sig_y_0"] = 20 * const.micro  # m
SIGMA_T = 400e-15 * 2.998e14 # s converted to micro
bp["sig_z_0"] = SIGMA_T * const.micro # m (I know its messy to convert to micro and then convert back to m - allows for either s or um sig_long)
bp["plat_rise"] =  0.0* SIGMA_T * const.speed_of_light #For plateau type beam (rise and fall time)
#We'll add a variation in longitudinal momentum so LPS plot differences are visible
bp["sig_pz_0"] = 1  # MeV/c
bp["lCorr_Fac"] = 0.95  # An arbitrary number, but want to convert to chirp at some point
bp["x_0"] = 0
bp["y_0"] = 0
bp["xp_0"] = 0
bp["yp_0"] = 0
bp["z_0"] = 0
TOTAL_CHARGE = 250e-12 #Bunch charge in C
bp["charge_per_macro"] = TOTAL_CHARGE / N_MACROS #Charge per macro

#Profile type either 'Gaussian', 'SkewGaussian' 'Uniform', 'Plataeu', or 'DoubleGauss'
bp["LongitudinalProfile"] = 'Gaussian'
myBeam = dbt.GaussianBeam()
myBeam.setBeamParameters(bp)
myBeam.calcWaistOpticFunctions()
myBeam.generate6DMacros(N_MACROS)

DWA_beam2 = DWA.DiWaCATOutput();
#We could either read the beam as before and then change the .beam object or:
DWA_beam2.beam = myBeam;

DWA_beam2.GetFieldValues(".\\Python_Tools\\Example_Scripts\\Example_DiWaCAT_Beam.h5")

#We then have an object just like before
#The DWA length can also be changed by the direct object
DWA_beam2._DielectricParameter['L'] = 0.2

DWA_beam2.TransportThroughDLW(nSteps = 3)
DWA_beam2.ClearLostParticles()


# ----------------- 3. Plotting DiWaCAT objects --------------------- #

#-------------------------- Fields -----------------------------#
#Accessing the fields is a little awkward - need to provide [x,y,z] positions
#Create an array of x=0, y=0 with z from -3 to +3
PlottingPositions = [];
Z = []
for i in range(100):
    Z.append((-3 + 6*i*1/100)*600)
    PlottingPositions.append([0, 0, (-3 + 6*i*1/100)*600e-15*3e8])
EzFunction = DWA_Beam.ReturnFieldFunctions()[2]
EzArray = EzFunction(PlottingPositions)

plt.plot(Z, EzArray * 1e-6)
plt.xlabel('t [fs]')
plt.ylabel('$E_z$ [MV/m]')
plt.show()

#------------------------- Beam Info -----------------------#
#Plotting beam information is done via BeamPlotter class in beam_tools or directly accessing .beam arrays (e.g. DWA_Beam.beam._beam['y'])

plotter1 = dbt.BeamPlotter();
plotter1.addBeamToPlotter('DiWaCAT Beam', DWA_Beam.beam)
#plotter1.addBeamToPlotter('Created Beam', DWA_beam2.beam)
plotter1.plotLPS()
plotter1.general_scatter_plot(x_val = 'x', x_scale = 'milli', y_val = 'y', y_scale = 'milli')

'''

