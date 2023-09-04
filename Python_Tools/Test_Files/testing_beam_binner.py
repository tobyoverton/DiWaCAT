#- * - coding: utf - 8 -*-
"""
Simple file that tests the DWA_beam_tools module
1) Creates a 6D gaussian beam
3) Bin the data onto a mesh
4) Smooth the data using either Gaussian, Uniform, Median smoothing method and plot comparisons
5) Output the data from the mesh to a h5 file - lots of meta data in this file!

@author: Tom Pacey, May 2021
"""
import os.path

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const

<<<<<<< Updated upstream:Python_Tools/Test_files/testing_beam_binner.py
<<<<<<< Updated upstream:Python_Tools/Test_files/testing_beam_binner.py
from Python_Tools.Modules import beam_tools as dbt
=======
import DWA_beam_tools as dbt
>>>>>>> Stashed changes:DWA_Beam_Tools/testing_beam_tools.py
=======
import DWA_beam_tools as dbt
>>>>>>> Stashed changes:DWA_Beam_Tools/testing_beam_tools.py

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


#---------------------------------------------------


#---------------------------------------------------------------------------
#---------------------- 3) Bin the data ------------------------------
#---------------------------------------------------------------------------


#Setup parameters for the binning mesh
cps = 3 #cells per sigma for transverse
cps_long = 4 #we'll use more cells longitudinally

Lx = 4.5 * myBeam.Sx # define over +/- 2.25 sigma - this works nicely for our Gaussian beam but we might want to be smarter here
Dx = myBeam.Sx/cps

Ly = 4.5 * myBeam.Sy
Dy = myBeam.Sy/(cps)

Lz = 6 * myBeam.Sz
Dz = myBeam.Sz/cps_long


# ----Perform the binning-----------------------------------------------------------------------------------------------

#Initalise instance of beamBinner with the myBeam object
binnedBeam = dbt.beamBinner(myBeam)

#and bin
binnedBeam.binTheBeam(Lx=Lx,Dx=Dx,Ly=Ly,Dy=Dy,Lz=Lz,Dz=Dz)

# Test hist full_kick_output
print("Shape of binned beam array is",binnedBeam.beamHist.shape)

#how much charge was captured on the mesh?
print("Fraction of total charge captured on this mesh is",binnedBeam.nMacroHist/myBeam.nMacros)

# Plotting from the binned bunch
#most particles are at the middle, so define the middle bin
midX = int(binnedBeam.beamHist.shape[0] / 2)
midY = int(binnedBeam.beamHist.shape[1] / 2)
midZ = int(binnedBeam.beamHist.shape[2] / 2)



#2D Contour map of longitudinal slice
xx, yy = np.meshgrid(binnedBeam.cenValsX / const.milli,binnedBeam.cenValsY / const.milli) #I don't really know what this does but need it for plotting

xx = xx.transpose() #Also seem to need this if x and y not equal
yy = yy.transpose()

#fig5, ax51 = plt.subplots(1,1)
#fig6, (ax61,ax62,ax63,ax64) = plt.subplots(1,4)

fig5 = plt.figure(tight_layout=False)
gs = fig5.add_gridspec(2,4, hspace=0.6)

fig5_ax1 = fig5.add_subplot(gs[0,:]) # type:axes.Axes
fig5_ax21 = fig5.add_subplot(gs[1,0]) # type:axes.Axes
fig5_ax22 = fig5.add_subplot(gs[1,1]) # type:axes.Axes
fig5_ax23 = fig5.add_subplot(gs[1,2]) # type:axes.Axes
fig5_ax24 = fig5.add_subplot(gs[1,3]) # type:axes.Axes

fig5_ax1.set_xlabel("dz [um]")
fig5_ax1.set_ylabel("Counts")
fig5_ax1.set_title(("Lineout in z from centre (x,y)"))


fig5_ax21.set_ylabel("y [mm]")
fig5_ax21.set_xlabel("x [mm]")
fig5_ax22.set_xlabel("x [mm]")
fig5_ax23.set_xlabel("x [mm]")
fig5_ax24.set_xlabel("x [mm]")


plt.figtext(0.5,0.45,"x-y plane at central (dz)",ha='center',size=12)


fig5_ax1.plot(binnedBeam.cenValsZ / const.micro,binnedBeam.beamHist[midX,midY,:],color = "black")

fig5_ax21.contourf(xx,yy,binnedBeam.beamHist[:,:,midZ],cmap = "Greys",levels=7)


#---------------------------------------------------------------------------
#---------------------- 4) Smooth the data and compare smoothing methods ------------------------------
#---------------------------------------------------------------------------


binnedBeam.smoothGaussian(sigma=0.75)
print("The beam was smoothed with",binnedBeam.smoothMethod,"with kwargs",binnedBeam.smooth_kwargs)

fig5_ax1.plot(binnedBeam.cenValsZ / const.micro,binnedBeam.smoothedHist[midX,midY,:],color = "blue")
fig5_ax22.contourf(xx,yy,binnedBeam.smoothedHist[:,:,midZ],cmap = "Blues",levels=7)

binnedBeam.smoothMedian(size=3)
print("The beam was smoothed with",binnedBeam.smoothMethod,"with kwargs",binnedBeam.smooth_kwargs)

fig5_ax1.plot(binnedBeam.cenValsZ / const.micro,binnedBeam.smoothedHist[midX,midY,:],color = "orange")
fig5_ax23.contourf(xx,yy,binnedBeam.smoothedHist[:,:,midZ],cmap = "Oranges",levels=7)


binnedBeam.smoothUniform(size=3)
print("The beam was smoothed with",binnedBeam.smoothMethod,"with kwargs",binnedBeam.smooth_kwargs)

fig5_ax1.plot(binnedBeam.cenValsZ / const.micro,binnedBeam.smoothedHist[midX,midY,:],color = "green")
fig5_ax24.contourf(xx,yy,binnedBeam.smoothedHist[:,:,midZ],cmap = "Greens",levels=7)


for ax in fig5.axes[1:]:
    ax.label_outer()

fig5_ax1.legend(["Raw binning","Gaussian Smoothing","Median Smoothing","Uniform Smoothing"])

plt.show()

#---------------------------------------------------------------------------
#---------------------- 5) Output the file ------------------------------
#---------------------------------------------------------------------------

#Now need to make a h5 file for full_kick_output,including macros (optional)
test_output_location = "export_files/test_beam_binner.h5"

print("checking if file", os.path.isfile(test_output_location))
print("Writing to full_kick_output file > ",test_output_location)
outputFileWritten = binnedBeam.write_hd5f_mesh_file(test_output_location,outputMacros=True, overwrite_exisiting_file=True,includeSmoothed=True,includeUnSmoothed=False)
if outputFileWritten:
    print("Output file written")

else:
    print("Output file not written")


##