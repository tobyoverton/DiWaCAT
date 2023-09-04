#- * - coding: utf - 8 -*-
"""
Simple file that tests the DWA_beam_tools module
1) Creates a 6D gaussian beam
2) Transports the beam through drift spaces and quads
3) Output the data from the beam to a h5 file
4) read that file back in and plot it again

@author: Tom Pacey, August 2020
"""
import sys

import numpy as np
import os
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
bp["eps_x_N"] = 1 * const.micro   # m rad
bp["eps_y_N"] = 1 * const.micro  # m rad
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
bp["charge_per_macro"] = 100e-12 / N_MACROS #C per macro

#bp["LongitudinalProfile"] = "Gaussian"

bp["LongitudinalProfile"] = "Plateau"
bp["plat_rise"] =  0.2* SIGMA_T * const.speed_of_light

#bp["LongitudinalProfile"] = "DoubleGauss"
#bp["offset"] = 4 * SIGMA_T * const.speed_of_light
#bp["rel_amp"] = 1
#bp["sig_z_2"] = 4 * SIGMA_T * const.speed_of_light

#bp["LongitudinalProfile"] = "SkewGaussian"
#bp["skew"] = -5

# initalise the beam instance
myBeam = dbt.GaussianBeam()

# apply the parameters we want
myBeam.setBeamParameters(bp)



# Calculate the optics functions at the waist
myBeam.calcWaistOpticFunctions()

print("Beam gamma is", myBeam._beam["gamma"])

# can change a parameter like this
#myBeam.beam["x_0"] = 200


# Finally generate the beam by specifying a number of macros
myBeam.generate6DMacros(N_MACROS)


print(myBeam.z[100], ",", myBeam.t[100])
'''

myBeam.print_all_beam_properties_table()

myBeamPlotter = dbt.BeamPlotter()

print("---------Transverse properties---------")
print(f"Beam sigma_x at waist is {myBeam.Sx/const.micro : .1f} um")
print(f"beam RMS eps_x_N is {myBeam.normalized_horizontal_emittance/const.micro: .1f} mm mrad")
print(f"beam 90% Gaussian eps_x_N is {myBeam.normalized_horizontal_emittance_90/const.micro: .1f} mm mrad")
print(f"beam 100% eps_x_N is {myBeam.normalized_mve_horizontal_emittance/const.micro: .1f} mm mrad")


print(f"Beam sigma_x at waist is {myBeam.Sx/const.micro : .1f} um")

#Insert percentile sig x here

print(f"Beam sigma_x 95% at waist is {myBeam.get_2nd_moment_percentile('x',95)/const.micro : .1f} um")
print(f"Beam sigma_x 50% at waist is {myBeam.get_2nd_moment_percentile('x',50)/const.micro : .1f} um")

#Short hand version
print(f"Beam sigma_x 95% at waist is {myBeam.Sx_95/const.micro : .1f} um")


#sys.exit()

plt.hist(myBeam.z,bins=100)

myBeamPlotter.addBeamToPlotter('At Waist',myBeam)
inital_LPS_plot = myBeamPlotter.plotLPS()


#---------------------------------------------------------------------------
#---------------------- 2) Drift, focus drift ------------------------------
#---------------------------------------------------------------------------

#Using quads to create asymetric beam, just for the sake of it (and easy idenitfication of x,y,z planes

#Setup a simple set of plots, x-y, x-xp,y-yp

# create a static reference to the waist beam
init_beam = myBeam.returnBeam()

# Using simple to implement optics functions to move the beam around a bit
myBeam.driftBeam(0.5) # a drift in m
print(f"Beam sigma_x after 0.5m drift is {myBeam.Sx/const.micro : .1f} um")
print(f"Beam sigma_y after 0.5m drift is {myBeam.Sy/const.micro : .1f} um")

quad_beam = myBeam.returnBeam()

#myBeam.shiftBeam('y',500*const.micro)

#=======
#print("Beam sigma_x after 5m drift is",myBeam.Sx/const.micro," um")
#>>>>>>> Stashed changes

beam_thick_quad = myBeam.returnBeam()

myBeam.thinLensFQuadBeam(0.25) # f quad with focal length m - set to focus x back down to waist
myBeam.driftBeam(0.25) # drift again
print(f"Beam sigma_x after quad and 0.5m drift is {myBeam.Sx/const.micro : .1f} um")
print(f"Beam sigma_y after quad and 0.5m drift is {myBeam.Sy/const.micro : .1f} um")



myBeamPlotter.clearAllBeams()
myBeamPlotter.addBeamToPlotter('Inital beam',init_beam)
myBeamPlotter.addBeamToPlotter('Beam at quad',quad_beam)
myBeamPlotter.addBeamToPlotter('Final beam',myBeam)

myBeamPlotter.plot_transverse_properties()

plt.show()

#apply a thick quad and track
#TODO add here work

myBeamPlotter.clearAllBeams()
myBeamPlotter.addBeamToPlotter('Beam at quad',quad_beam)



#plt.close()

print("------------------------")
print("outputting to file")

test_outputfile ="export_files/test_gaussian_beam_output.h5"


check_output = myBeam.write_HDF5_beam_file(test_outputfile,overwrite_existing_file=False)

print("Was file written: ",check_output)

print("------------------------")
print("reading file back in")
print("File exists: ",os.path.isfile(test_outputfile))

readback_beam = dbt.BeamFromFile()
readback_beam.read_HDF5_beam_file(test_outputfile)

myBeamPlotter.clearAllBeams()
myBeamPlotter.addBeamToPlotter("full_kick_output beam",myBeam)
myBeamPlotter.addBeamToPlotter("read in beam", readback_beam)

myBeamPlotter.plot_transverse_properties()

plt.show()
'''