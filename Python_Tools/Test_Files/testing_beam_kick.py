#- * - coding: utf - 8 -*-
"""
Simple file that tests the DWA_beam_tools module sub class beam kick


@author: Tom Pacey, March 2021
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const

from Python_Tools.Modules import beam_tools as dbt

from scipy.stats import gaussian_kde



#Initalise a uniform-gauss beam

# How many macros?!?
N_MACROS = 1000

# Make empty dict to hold the parameters we want the constructor class to call
bp = {}

bp["pz_MeV"] = 35  # MeV/c
bp["eps_x_N"] = 5 * const.micro   # m rad
bp["eps_y_N"] = 1 * const.micro  # m rad
bp["sig_x_0"] = 100 * const.micro  # m
bp["sig_y_0"] = 100 * const.micro  # m
SIGMA_T = 500e-15 # s
bp["sig_z_0"] = SIGMA_T * const.speed_of_light # m
bp["sig_pz_0"] = 0.05  # MeV/c
bp["lCorr_Fac"] = 0.0  # An arbitrary number, but want to convert to chirp at some point
bp["x_0"] = 0
bp["y_0"] = 0
bp["xp_0"] = 0
bp["yp_0"] = 0
bp["z_0"] = 0
bp["charge"] = 100e-12 #C

bp["LongitudinalProfile"] = "Plateau"
bp["plat_rise"] =  0.2* SIGMA_T * const.speed_of_light


# initalise the beam instance
myBeam = dbt.GaussianBeam()

# apply the parameters we want
myBeam.setBeamParameters(bp)

# Calculate the optics functions at the waist
myBeam.calcWaistOpticFunctions()

# Finally generate the beam by specifying a number of macros
myBeam.generate6DMacros(N_MACROS)

#


kp = {}
kp['freq'] = 0.2 * const.speed_of_light/(SIGMA_T * const.speed_of_light) # Hz
kp['amplitude'] = 0.2 # MV/m
kp['length'] = 1 # m
kp['phase'] = const.pi / 2 #rad
#kp['transverse_profile'] = 'Uniform'

kp['transverse_profile'] = 'Gaussian'
kp['sig_r'] = 100 * const.micro #m

#just like this
kickedBeam = dbt.BeamKicker(myBeam)
kickedBeam.copyInputToOutput()
kickedBeam.setKickParameters(kp)
kickedBeam.kickOutput()
outputBeam = kickedBeam.returnOutputBeam()
del kickedBeam #one can do this now


fig2, ax1 = plt.subplots()  # A new graph
ax1.plot(outputBeam.zn /const.micro, outputBeam.cpz / const.mega, '.', color='blue', alpha=0.9,markersize=1)
ax1.plot(myBeam.zn /const.micro, myBeam.cpz / const.mega, '.', color='red', alpha=0.9,markersize=1)
ax1.set_xlabel("zeta [um]")
ax1.set_ylabel("pz [MeV/c]")
plt.show()

#TODO Shunt this into the beam objected method

test_output = outputBeam.cpz/const.mega
test_input = myBeam.cpz/const.mega

cpzVals = np.linspace(min(test_output),max(test_output),100)
inputKDE = gaussian_kde(test_input)
outputKDE = gaussian_kde(test_output)

inputSmoothed = inputKDE.evaluate(cpzVals)
outputSmoothed = outputKDE.evaluate(cpzVals)

fig3,ax1 = plt.subplots()
ax1.plot(cpzVals,outputSmoothed,color = 'blue')
ax1.plot(cpzVals,inputSmoothed,color = 'red')
plt.show()



