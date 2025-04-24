# -*- coding: utf-8 -*-
"""
Using the DiWaCAT tools to make a drive-main setup

Need to use a script since we want a door-step drive and gaussian main
Needs 3 beams and a custom mesh to make sure all info is captured
"""

import Python_Tools.Modules.diwacat_tools as DWA
import Python_Tools.Modules.beam_tools as dbt
import numpy as np
import scipy.constants as const


beam1 = dbt.GaussianBeam()
beam2 = dbt.GaussianBeam()
beam3 = dbt.GaussianBeam()

TOTAL_CHARGE=250e-12
TOTAL_MACROS=250e3
#Make empty dictionary for the beam parameters
'''
bp1 = {}
N_MACROS1 = round(TOTAL_MACROS*0.5)
bp1["pz_MeV"] = 250
bp1["eps_x_N"] = 1e-6
bp1["eps_y_N"] = 1e-6
bp1["sig_x_0"] = 50e-6
bp1["sig_y_0"] = 50e-6
SIGMA_LONG1 = 0.8e-12 * const.c
bp1["sig_z_0"] = SIGMA_LONG1
bp1["sig_pz_0"] = 0
bp1["lCorr_Fac"] = 0#1  # An arbitrary number, but want to convert to chirp at some point
bp1["x_0"] = 0
bp1["y_0"] = 0
bp1["xp_0"] = 0
bp1["yp_0"] = 0
bp1["z_0"] = (1e-12 + 0.15e-12)*const.c
TOTAL_CHARGE1 = TOTAL_CHARGE*0.5 #Bunch charge in C
bp1["charge_per_macro"] = TOTAL_CHARGE1 / N_MACROS1 #Charge per macro
bp1["LongitudinalProfile"] = 'SkewGaussian'
bp1["skew"] = -8.0

bp2 = {}
N_MACROS2 = round(TOTAL_MACROS*0.5)
bp2["pz_MeV"] = 250
bp2["eps_x_N"] = 1e-6
bp2["eps_y_N"] = 1e-6
bp2["sig_x_0"] = 50e-6
bp2["sig_y_0"] = 50e-6
SIGMA_LONG2 = 2e-12 * const.c
bp2["sig_z_0"] = SIGMA_LONG1
bp2["sig_pz_0"] = 0
bp2["lCorr_Fac"] = 0#1  # An arbitrary number, but want to convert to chirp at some point
bp2["x_0"] = 0
bp2["y_0"] = 0
bp2["xp_0"] = 0
bp2["yp_0"] = 0
bp2["z_0"] = 0
TOTAL_CHARGE2 = 1e-9#TOTAL_CHARGE*0.5 #Bunch charge in C
bp2["charge_per_macro"] = TOTAL_CHARGE2 / N_MACROS2 #Charge per macro
bp2["LongitudinalProfile"] = 'Plateau'
bp2["plat_rise"] = 100e-15 * const.c
'''
bp1 = {}
N_MACROS1 = round(TOTAL_MACROS*0.8)
bp1["pz_MeV"] = 250
bp1["eps_x_N"] = 1e-6
bp1["eps_y_N"] = 1e-6
bp1["sig_x_0"] = 50e-6
bp1["sig_y_0"] = 50e-6
SIGMA_LONG1 = 2e-12 * const.c
bp1["sig_z_0"] = SIGMA_LONG1
bp1["sig_pz_0"] = 0
bp1["lCorr_Fac"] = 0#1  # An arbitrary number, but want to convert to chirp at some point
bp1["x_0"] = 0
bp1["y_0"] = 0
bp1["xp_0"] = 0
bp1["yp_0"] = 0
bp1["z_0"] = 0#(1e-12 + 0.15e-12)*const.c
TOTAL_CHARGE1 = TOTAL_CHARGE*0.8 #Bunch charge in C
bp1["charge_per_macro"] = TOTAL_CHARGE1 / N_MACROS1 #Charge per macro
bp1["LongitudinalProfile"] = 'SkewGaussian'
bp1["skew"] = -8.0

bp2 = {}
N_MACROS2 = round(TOTAL_MACROS*0.2)
bp2["pz_MeV"] = 250
bp2["eps_x_N"] = 1e-6
bp2["eps_y_N"] = 1e-6
bp2["sig_x_0"] = 50e-6
bp2["sig_y_0"] = 50e-6
SIGMA_LONG2 = 0.2e-12 * const.c
bp2["sig_z_0"] = SIGMA_LONG1
bp2["sig_pz_0"] = 0
bp2["lCorr_Fac"] = 0#1  # An arbitrary number, but want to convert to chirp at some point
bp2["x_0"] = 0
bp2["y_0"] = 0
bp2["xp_0"] = 0
bp2["yp_0"] = 0
bp2["z_0"] = -0.9e-3
TOTAL_CHARGE2 = TOTAL_CHARGE*0.2 #Bunch charge in C
bp2["charge_per_macro"] = TOTAL_CHARGE2 / N_MACROS2 #Charge per macro
bp2["LongitudinalProfile"] = 'SkewGaussian'
bp2["skew"] = -8

beam1.setBeamParameters(bp1)
beam1.calcWaistOpticFunctions()
beam1.generate6DMacros(N_MACROS1)
beam2.setBeamParameters(bp2)
beam2.calcWaistOpticFunctions()
beam2.generate6DMacros(N_MACROS2)


beam1._beam['x'] = np.append(beam1._beam['x'], beam2._beam['x'])
beam1._beam['y'] = np.append(beam1._beam['y'], beam2._beam['y'])
beam1._beam['z'] = np.append(beam1._beam['z'], beam2._beam['z'])
beam1._beam['px'] = np.append(beam1._beam['px'], beam2._beam['px'])
beam1._beam['py'] = np.append(beam1._beam['py'], beam2._beam['py'])
beam1._beam['pz'] = np.append(beam1._beam['pz'], beam2._beam['pz'])
beam1._beam['t'] = np.append(beam1._beam['t'], beam2._beam['t'])
beam1._beam['charge'] = np.append(beam1._beam['charge'], beam2._beam['charge'])
beam1._beam['z'] = beam1.zn
beam1._beam['t'] = beam1.get_time_coord()

#Add the witness bunch
bp3 = {}
N_MACROS3 = round(TOTAL_MACROS*0.05)
bp3["pz_MeV"] = 250
bp3["eps_x_N"] = 1e-6
bp3["eps_y_N"] = 1e-6
bp3["sig_x_0"] = 50e-6
bp3["sig_y_0"] = 50e-6
SIGMA_LONG3 = 0.05e-12 * const.c
bp3["sig_z_0"] = SIGMA_LONG3
bp3["sig_pz_0"] = 0
bp3["lCorr_Fac"] = 0#1  # An arbitrary number, but want to convert to chirp at some point
bp3["x_0"] = 0
bp3["y_0"] = 0
bp3["xp_0"] = 0
bp3["yp_0"] = 0
bp3["z_0"] = 1e-3
TOTAL_CHARGE3 = TOTAL_CHARGE*0.05 #Bunch charge in C
bp3["charge_per_macro"] = TOTAL_CHARGE3 / N_MACROS3 #Charge per macro
bp3["LongitudinalProfile"] = 'Gaussian'

beam3.setBeamParameters(bp3)
beam3.calcWaistOpticFunctions()
beam3.generate6DMacros(N_MACROS3)

beam1._beam['x'] = np.append(beam1._beam['x'], beam3._beam['x'])
beam1._beam['y'] = np.append(beam1._beam['y'], beam3._beam['y'])
beam1._beam['z'] = np.append(beam1._beam['z'], beam3._beam['z'])
beam1._beam['px'] = np.append(beam1._beam['px'], beam3._beam['px'])
beam1._beam['py'] = np.append(beam1._beam['py'], beam3._beam['py'])
beam1._beam['pz'] = np.append(beam1._beam['pz'], beam3._beam['pz'])
beam1._beam['t'] = np.append(beam1._beam['t'], beam3._beam['t'])
beam1._beam['charge'] = np.append(beam1._beam['charge'], beam3._beam['charge'])
beam1._beam['z'] = beam1.zn
beam1._beam['t'] = beam1.get_time_coord()


Lx = 50e-6 * 5
Ly = 50e-6 * 5
Lz = 8e-12 * const.c
Dx = 20e-6
Dy = 20e-6
Dz = 15e-15 * const.c

binnedBeam = dbt.beamBinner(beam1)
'''
We want a custom mesh -> copy binTheBeam method but change maxZ (waste of time added loads of points before the beam itself)
'''
xMin = binnedBeam.beam.Mx - (Lx / 2)
xMax = binnedBeam.beam.Mx + (Lx / 2)
nBinsx = int(Lx / Dx)

yMin = binnedBeam.beam.My - (Ly / 2)
yMax = binnedBeam.beam.My + (Ly / 2)
nBinsy = int(Ly / Dy)

zMin = binnedBeam.beam.Mzn - 3*binnedBeam.beam.Sz
zMax = binnedBeam.beam.Mzn + Lz
nBinsz = int((zMax-zMin) / Dz)

# Sanity check mesh parameters
binnedBeam.mesh['cps_x'] = binnedBeam.beam.Sx / Dx
binnedBeam.mesh['cps_y'] = binnedBeam.beam.Sy / Dy
binnedBeam.mesh['cps_z'] = binnedBeam.beam.Sz / Dz

binnedBeam.mesh['xRange'] = (xMin, xMax)
binnedBeam.mesh['yRange'] = (yMin, yMax)
binnedBeam.mesh['zRange'] = (zMin, zMax)

# Make 3D hist to bin the beam
# Must use beam.zn here!
H, edges = np.histogramdd((binnedBeam.beam.x, binnedBeam.beam.y, binnedBeam.beam.zn),
                          bins=(nBinsx, nBinsy, nBinsz),
                          range=((xMin, xMax), (yMin, yMax), (zMin, zMax)))

binnedBeam.mesh['beamHist'] = H
binnedBeam.mesh['histEdges'] = edges

# Calculate the 'central' mesh points from the bin edges
binnedBeam.mesh['cenValsX'] = np.empty(len(edges[0][:-1]))  # empty array
binnedBeam.mesh['dx'] = edges[0][1] - edges[0][0]
for i, val in enumerate(edges[0][:-1]):
    binnedBeam.mesh['cenValsX'][i] = val + (binnedBeam.mesh['dx'] / 2)

binnedBeam.mesh['cenValsY'] = np.empty(len(edges[1][:-1]))  # empty array
binnedBeam.mesh['dy'] = edges[1][1] - edges[1][0]
for i, val in enumerate(edges[1][:-1]):
    binnedBeam.mesh['cenValsY'][i] = val + (binnedBeam.mesh['dy'] / 2)

binnedBeam.mesh['cenValsZ'] = np.empty(len(edges[2][:-1]))  # empty array
binnedBeam.mesh['dz'] = edges[2][1] - edges[2][0]
for i, val in enumerate(edges[2][:-1]):
    binnedBeam.mesh['cenValsZ'][i] = val + (binnedBeam.mesh['dz'] / 2)

binnedBeam.mesh['nMacroHist'] = np.sum(binnedBeam.mesh['beamHist'], axis=(0, 1, 2))


path = './../../../DWA_FEBE/DWA_Demonstrator_TwoBeams/DoubleTriangle_Total4ps_TwoBunch2.h5'
binnedBeam.write_hd5f_mesh_file(path,outputMacros=True, overwrite_existing_file=True,includeSmoothed=False,includeUnSmoothed=True)