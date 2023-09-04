"""
@author T Pacey
"""


import sys

import numpy as np
import os
import matplotlib.pyplot as plt
import scipy.constants as const

from Python_Tools.Modules import beam_tools as dbt
from Python_Tools.Modules import tracking_tools
from Python_Tools.Modules import  general_tools as gent


#region create beam
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

# Finally generate the beam by specifying a number of macros
myBeam.generate6DMacros(N_MACROS)
#endregion

lattice=[
    {'drift':{'L':1,"step_size":0.2}},
    {'thick_quad':{'L':0.2,'K':5,"step_size":0.02}},
    {'drift':{'L':1.5,"step_size":0.2}},
    {'thin_quad':{'f':-0.5}},
    {'drift': {'L': 1.5, "step_size": 0.2}}

]

parameters = [
    'Sy','Sx','Mz',"beta_x","alpha_x","gamma_x","normalized_horizontal_emittance"
]
print(f"Inital sigma_x = {myBeam.Sx}")
tracker = tracking_tools.LatticeTracker(myBeam,lattice)
tracker.set_saved_parameters(parameters)

print(tracker.global_step_value)
tracker.global_step_value = 0.1
tracker.perform_track()

print(f"Final sigma_x = {myBeam.Sx}")

plt.figure(1)
plt.plot(tracker.saved_parameters['Mz'],np.asarray(tracker.saved_parameters['Sx'])/const.micro,marker='.',linestyle = '-')
plt.plot(tracker.saved_parameters['Mz'],np.asarray(tracker.saved_parameters['Sy'])/const.micro,marker='.',linestyle = '-')
#plt.plot(parameters['Mz'],np.asarray(parameters['Sy'])/const.milli,marker='.',linestyle = '')
plt.title("Beam Sig Vals")
plt.xlabel("s (m)")
plt.ylabel("$\sigma_i$ (um)")

plt.figure(2)
plt.plot(tracker.saved_parameters['Mz'],np.asarray(tracker.saved_parameters['alpha_x']),marker='.',linestyle = '-')
plt.plot(tracker.saved_parameters['Mz'],np.asarray(tracker.saved_parameters['beta_x']),marker='.',linestyle = '-')
plt.plot(tracker.saved_parameters['Mz'],np.asarray(tracker.saved_parameters['gamma_x']),marker='.',linestyle = '-')

plt.title("Twiss_vals")
plt.xlabel("s (m)")
plt.ylabel("m")

plt.figure(3)
plt.plot(tracker.saved_parameters['Mz'],np.asarray(tracker.saved_parameters['normalized_horizontal_emittance'])/const.micro,marker='.',linestyle = '-')


plt.show()