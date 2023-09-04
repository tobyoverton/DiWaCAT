#- * - coding: utf - 8 -*-
"""
Simple file that tests the slicing in beam tools

@author: Tom Pacey, August 2020
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const

from Python_Tools.Modules import beam_tools as dbt

import matplotlib.axes._axes as axes # this for PyCharm completion for axes objects
from matplotlib.lines import Line2D # this for custom legends


#from scipy.stats import gaussian_kde
#from scipy.ndimage import gaussian_filter,median_filter,uniform_filter
#import h5py

# ---------------------- 1) create the beam ---------------------------

# How many macros?!?
N_MACROS = 32000

# Make empty dict to hold the parameters we want the constructor class to call
bp = {}

bp["pz_MeV"] = 35  # MeV/c
bp["eps_x_N"] = 5 * const.micro   # m rad
bp["eps_y_N"] = 5 * const.micro  # m rad
bp["sig_x_0"] = 100 * const.micro  # m
bp["sig_y_0"] = 100 * const.micro  # m
SIGMA_T = 200e-15 # s
bp["sig_z_0"] = SIGMA_T * const.speed_of_light # m
bp["sig_pz_0"] = 0.05  # MeV/c
bp["lCorr_Fac"] = 0.5  # An arbitrary number, but want to convert to chirp at some point
bp["x_0"] = 0
bp["y_0"] = 0
bp["xp_0"] = 0
bp["yp_0"] = 0
bp["z_0"] = 0
bp["charge"] = 100e-12 #C

bp["LongitudinalProfile"] = "Gaussian"

# initalise the beam instance
myBeam = dbt.GaussianBeam()

# apply the parameters we want
myBeam.setBeamParameters(bp)

# Calculate the optics functions at the waist
myBeam.calcWaistOpticFunctions()

# Finally generate the beam by specifying a number of macros
myBeam.generate6DMacros(N_MACROS)

myBeam.driftBeam(1)

fig1, (ax1,ax2,ax3) = plt.subplots(1,3)  # type:axes.Axes

#TODO wrap this up into beam plotter!
plt.figure(1)
ax1.plot(myBeam.x / const.milli, myBeam.y / const.milli, '.', color='blue', alpha=0.3,markersize = 0.1)
ax2.plot(myBeam.x / const.milli, myBeam.xp / const.milli, '.', color='blue', alpha=0.3,markersize = 0.1)
ax3.plot(myBeam.y / const.milli, myBeam.yp / const.milli, '.', color='blue', alpha=0.3,markersize = 0.1)

# Make the plots prettier
ax1.set_xlabel('x [mm]')
ax1.set_ylabel('y [mm]')
ax1.set_xlim(-5,5)
ax1.set_ylim(-5,5)
ax1.set_aspect(1)

ax2.set_xlabel('x [mm]')
ax2.set_ylabel('x\' [mm]')
ax2.set_xlim(-5,5)
ax2.set_ylim(-5,5)
ax2.set_aspect(1)

ax3.set_xlabel('y [mm]')
ax3.set_ylabel('y\' [mm]')
ax3.set_xlim(-5,5)
ax3.set_ylim(-5,5)
ax3.set_aspect(1)

plt.show()

horizontal_slicer = dbt.BeamSlicer(myBeam, 'x')

horizontal_slicer.set_slice_parameters(min_p=-1e-3,max_p=1e-3,N_slices=10,slice_width=50e-6)
horizontal_slicer.perform_slicing()

print(len(horizontal_slicer.sliced_beams))

fig2,(ax1,ax2) = plt.subplots(1,2) # type:axes.Axes

for beam in horizontal_slicer.sliced_beams:
    ax1.scatter(beam.x / const.milli,
             beam.y / const.milli,color='blue', alpha=0.1,marker='.',s=1)

    beam.driftBeam(0.5)
    ax2.scatter(beam.x / const.milli,
             beam.y / const.milli,color='blue', alpha=0.1,marker='.',s=1)

ax1.set_xlabel('x [mm]')
ax1.set_ylabel('y [mm]')
ax1.set_xlim(-1,1)
ax1.set_ylim(-5,5)

plt.show()
