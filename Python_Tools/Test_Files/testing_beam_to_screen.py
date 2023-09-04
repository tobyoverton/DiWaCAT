"""
@author Tom Pacey
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
#from scipy.stats import gaussian_kde
#from scipy.ndimage import gaussian_filter

from Python_Tools.Modules import beam_tools as dbt
from Python_Tools.Modules import beam_to_screen as b2s

#region Building TestBeam
#--------------- TEST BEAM---------------------

# How many macros?!?
N_MACROS = 10000


# Make empty dict to hold the parameters we want the constructor class to call
bp = {}

bp["pz_MeV"] = 35  # MeV/c
bp["eps_x_N"] = 1 * const.micro   # m rad
bp["eps_y_N"] = 1 * const.micro  # m rad
bp["sig_x_0"] = 1000 * const.micro  # m
bp["sig_y_0"] = 1000 * const.micro  # m
SIGMA_T = 200e-15 # s
bp["sig_z_0"] = SIGMA_T * const.speed_of_light # m
bp["sig_pz_0"] = 0.05  # MeV/c
bp["lCorr_Fac"] = 0.5  # An arbitrary number, but want to convert to chirp at some point
bp["x_0"] = 0 * const.micro
bp["y_0"] = 0
bp["xp_0"] = 0
bp["yp_0"] = 0
bp["z_0"] = 0
bp["charge_per_macro"] = 100e-12 / N_MACROS #C per macro

#bp["LongitudinalProfile"] = "Gaussian"

bp["LongitudinalProfile"] = "Plateau"
bp["plat_rise"] =  0.2* SIGMA_T * const.speed_of_light

print("Initialising test beam")
# initalise the beam instance
myBeam = dbt.GaussianBeam()
# apply the parameters we want
myBeam.setBeamParameters(bp)
# Calculate the optics functions at the waist
myBeam.calcWaistOpticFunctions()
myBeam.generate6DMacros(N_MACROS)
print("beam complete")
#endregion

#region Screen dict

#Establish a screen dict
# Some example reasonable parameters

screenDict = {}
screenDict["n_pix_x"] = 2560
screenDict["n_pix_y"] = 2160
screenDict["pix_size_x"] = 20 * const.micro
screenDict["pix_size_y"] = 20 * const.micro
screenDict["central_pix_x"] = 1500
screenDict["central_pix_y"] = 1000
screenDict["base_constant"] = 0
screenDict["beam_pix_pad"] = 5

screenDict["physical_res"] = 50 * const.micro

#Optional
screenDict["KDE_method"] = "fastKDE"


#endregion

testScreen = b2s.BeamToScreen(screenDict,myBeam)

print("Screen range",testScreen.x_vals.ptp())
print("Beam range",testScreen.beam_x_vals.ptp())

print(len(testScreen.x_vals))
print("Num pix screen =",np.size(testScreen.xx_vals))
print("Num pix beam =", np.size(testScreen.beam_pix_vals_pdf))

#Do the digitisation
bit_depth = 12
peak_pix = 2**bit_depth
testScreen.set_digitise_pix_by_max_value(bit_depth=bit_depth, max_pix_val=peak_pix)

print("digi run 1")
print("---------------------")
print("---------------------")


print("Beam charge ",testScreen.beam.charge)

print(np.max(testScreen.pix_vals_q_dens))
print(np.max(testScreen.pix_vals_digitised))

print(peak_pix)

#This beam saturates - this is the density to saturate
saturation_charge_dens = np.max(testScreen.pix_vals_q_dens)
print(saturation_charge_dens)



#Make a second screen that takes charge density to calcualte peak pix

#make a second beam that has tighter focus

#region second beam with smaller size in vert plane
# Make empty dict to hold the parameters we want the constructor class to call
bp2 = bp.copy()

bp2["sig_y_0"] = 500 * const.micro

print("Init second beam")
# initalise the beam instance
myBeam2 = dbt.GaussianBeam()
# apply the parameters we want
myBeam2.setBeamParameters(bp2)
# Calculate the optics functions at the waist
myBeam2.calcWaistOpticFunctions()
myBeam2.generate6DMacros(N_MACROS)
print("beam complete")
#endregion

#region third beam with large size and reduced charge
# Make empty dict to hold the parameters we want the constructor class to call
bp3 = bp.copy()

bp3["sig_y_0"] = 2000 * const.micro
bp3["sig_x_0"] = 2000 * const.micro
bp3["sig_y_0"] = 2000 * const.micro
bp3["charge_per_macro"] = 50e-12 / N_MACROS #C per macro


print("Init second beam")
# initalise the beam instance
myBeam3 = dbt.GaussianBeam()
# apply the parameters we want
myBeam3.setBeamParameters(bp3)
# Calculate the optics functions at the waist
myBeam3.calcWaistOpticFunctions()
myBeam3.generate6DMacros(N_MACROS)
print("beam complete")
#endregion



secondScreen = b2s.BeamToScreen(screenDict,myBeam2)

thirdScreen = b2s.BeamToScreen(screenDict,myBeam3)



#Here we are digitising the pixels, but setting the charge density for saturation to what it was in the first screen
secondScreen.set_digitise_pix_by_charge_dens(bit_depth=bit_depth,saturation_charge_dens =saturation_charge_dens)
thirdScreen.set_digitise_pix_by_charge_dens(bit_depth=bit_depth,saturation_charge_dens =saturation_charge_dens)

print(np.max(testScreen.pix_vals_digitised))
print(np.max(secondScreen.pix_vals_digitised))
print(np.max(thirdScreen.pix_vals_digitised))


#region full_kick_output plots

#Plot of the beam values at each virtual pix in real units
fig1,ax1 = plt.subplots(1,3)
ax1[0].pcolormesh(testScreen.beam_x_vals, testScreen.beam_y_vals, testScreen.beam_processed_pix, shading='auto',cmap = 'jet')
ax1[1].pcolormesh(secondScreen.beam_x_vals, secondScreen.beam_y_vals, secondScreen.beam_processed_pix, shading='auto',cmap = 'jet')
ax1[2].pcolormesh(thirdScreen.beam_x_vals, thirdScreen.beam_y_vals, thirdScreen.beam_processed_pix, shading='auto',cmap = 'jet')

ax_limit = np.max(thirdScreen.beam_x_vals)
print("ax limit  ",ax_limit)
for ax in ax1:
    ax.set_xlim(-1*ax_limit,ax_limit)
    ax.set_ylim(-1*ax_limit,ax_limit)

#plot of beam values in the image, normalised units
fig2,ax2 = plt.subplots(2,3)
ax2[0,0].imshow(testScreen.pix_vals_normalised,aspect=testScreen.x_vals.ptp()/testScreen.y_vals.ptp(),cmap = 'jet')
ax2[0,1].imshow(secondScreen.pix_vals_normalised,aspect=secondScreen.x_vals.ptp()/secondScreen.y_vals.ptp(),cmap = 'jet')
ax2[0,2].imshow(thirdScreen.pix_vals_normalised,aspect=secondScreen.x_vals.ptp()/secondScreen.y_vals.ptp(),cmap = 'jet')

#plot of values where saturation level is set by image 1
# must use clim here to get correct visualisation, but the values are protected in the array full_kick_output
ax2[1,0].imshow(testScreen.pix_vals_digitised,aspect=testScreen.x_vals.ptp()/testScreen.y_vals.ptp(),cmap = 'jet',clim=(0,peak_pix))
ax2[1,1].imshow(secondScreen.pix_vals_digitised,aspect=secondScreen.x_vals.ptp()/secondScreen.y_vals.ptp(),cmap = 'jet',clim=(0,peak_pix))
ax2[1,2].imshow(thirdScreen.pix_vals_digitised,aspect=secondScreen.x_vals.ptp()/secondScreen.y_vals.ptp(),cmap = 'jet',clim=(0,peak_pix))

plt.show()

#endregion

