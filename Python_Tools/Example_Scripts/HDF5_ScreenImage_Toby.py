#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 14:08:11 2022

@author: mbcx4to2
"""

#Reads in beam file, kicks to a screen and outputs the screen image

import os
import os, sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const

import matplotlib.axes._axes as axes # this for PyCharm completion for axes objects
from matplotlib.lines import Line2D # this for custom legends

from Python_Tools.Modules import beam_tools as dbt, read_beam_file as rbf, beam_to_screen as b2s

INPUT_FILE = "/Users/user/Documents/PhD Documents/Streaker_Work/Streaker_Paper/Simulations/Beams_Out/Emittance1um_Width50um/Offset800/200fsGaus_First.h5"

# open a beam file and get it

myBeam = dbt.BeamFromFile()

myBeam.read_WakeCode_beam_file(filename=INPUT_FILE)

distanceToScreen = 5;


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

myBeam.driftBeam(distanceToScreen);

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
print(np.max(testScreen.pix_vals_digitised))

fig1,ax1 = plt.subplots()
ax1.pcolormesh(testScreen.beam_x_vals, testScreen.beam_y_vals, testScreen.beam_processed_pix, shading='auto',cmap = 'jet')

ax_limit = np.max(testScreen.beam_x_vals)
print("ax limit  ",ax_limit)
ax1.set_xlim(-0.025,0.025)
ax1.set_ylim(-0.005,0.045)

plt.show()

