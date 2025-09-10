'''
Example_DWA_Transport.py
T Overton - Sept25

Example tutorial to initialise a DWA stage using DiWaCAT functions, calculate the fields, and transport through each DLW
Beam saved after each DLW in the stage. Field recalculated each time
'''

#Some of these libraries might not be needed - just copied from another example script
import sys
import numpy as np
import math
import scipy.constants as const
#sys.path.append("PATH\\TO\\DiWaCAT")
from Python_Tools.Modules import diwacat_tools as DWA
from Cpp_Source import DiWaCAT_FieldCalc as DWACalc
from Python_Tools.Modules import beam_tools as dbt

SAVEFOLDER = 'BeamsOut/' #Don't forget the backslash at the end

#---------------------------------------------------------------#
#----- 1. Initialise the repeating pattern in the cell ---------#
#---------------------------------------------------------------#

#First we need to initialise a DiWaCAT object
#We will read in a beam file here - see other example scripts if you want to make the beam
DWA_Beam = DWA.DiWaCATOutput()
DWA_Beam.ReadFromFile("beamfile.h5")


#Next we define arrays of variables that will change from DLW to DLW and those that will be constant
#Remember values need to be in cm - I'll add a 1e2 to make the values clearer to read
a_gap = 1e-3 * 1e2
delta = 200e-6 * 1e2
dlw_length = 10e-2 #Not in cm
permitivity = 3.75
permeability = 1
width_cm = 1
max_z = np.max(DWA_Beam.beam.z) #Not in cm. If you want to see the field behind the bunch just extend this number to whatever
sN = 40
sI = 40
ModePrec = 0.01
ModeAcc = 0.01

#For the changing variables we need to know how many stages there will be and how many DLWs per stage
#For this example we'll have 3 stages and each stage will have 10 DLWs (HVVHHVVHHV)
N_STAGES = 3
N_DLW = 12

orientationhv = ['h', 'v', 'v', 'h', 'h', 'v', 'v', 'h', 'h', 'v', 'v', 'h']
orientationvh = ['v', 'h', 'h', 'v', 'v', 'h', 'h', 'v', 'v', 'h', 'h', 'v']
orientationlist = [orientationhv, orientationvh, orientationhv]
#Offsets are in the DLW geometry so y0 in a V DLW is actually a horizontal offset in the beam frame
#These are just example offsets so the for loops work later
y0_DWA1 = 10e-6
y0_DWA2 = -52e-6
y0_DWA3 = -34e-6
offsetlist = [y0_DWA1, y0_DWA2, y0_DWA3]

#---------------------------------------------------------#
#------- 2. Now loop through the stages and DLWs ---------#
#---------------------------------------------------------#

for m in range(N_STAGES):
    for n in range(N_DLW):
        offset = offsetlist[m]
        orientation = orientationlist[m][n]
        #Assume we don't want an offset for the vertical DLWs
        if orientation == 'v':
            offset = 0
		#I will add a cheat here - let's say we won't recalculate the field if the last DLW geometry is the same
		#This way we recalculate the fields when the orientation changes - saves a good few calculations
        recalculate = True
        if n > 0:
            if orientationlist[m][n] == orientationlist[m][n-1]:
                recalculate = False
		
        #--------------------------------------------#
        #---3. Calculate the fields if neccessary ---#
        #--------------------------------------------#
        if recalculate:
            #We need to have a beam mesh to calculate field
            #After a few field values just check the field looks sensible - this 
            Lx = 6 * DWA_Beam.beam.Sx
            Ly = 6 * DWA_Beam.beam.Sy
            totalLength = max_z - np.min(DWA_Beam.beam.z)
            Dx = 50e-6
            Dy = 50e-6
            Dz = totalLength * 1/25

            binnedBeam = dbt.beamBinner(DWA_Beam.beam)
            
            xMin = binnedBeam.beam.Mx - (Lx / 2)
            xMax = binnedBeam.beam.Mx + (Lx / 2)
            nBinsx = int(Lx / Dx)

            yMin = binnedBeam.beam.My - (Ly / 2)
            yMax = binnedBeam.beam.My + (Ly / 2)
            nBinsy = int(Ly / Dy)

            zMin = min(binnedBeam.beam.z) - Dz
            zMax = max_z
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
            H, edges = np.histogramdd((binnedBeam.beam.x, binnedBeam.beam.y, binnedBeam.beam.z),
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
            
            
			#The field calculation is a different class - we can call this and then calculate the field each time
            #We'll save the field calculated and then reload that file for the beam - a little convoluted but it works
            calculation = DWACalc.DiWaCAT_FieldCalc()
            #Setting macros for field calculations
            calculation._beam = DWA_Beam.beam._beam
            xReduced = []
            yReduced = []
            zReduced = []
            Charge3D = []
            #Use a flag to ignore deprecation warnings - caused by the h5py library so out of our control
            Charge3D = np.array(binnedBeam.beamHist)
            xReduced = np.array(binnedBeam.cenValsX) 
            yReduced = np.array(binnedBeam.cenValsY) 
            zReduced = np.array(binnedBeam.cenValsZ)
            calculation._macros['ChargePerMacro'] = np.array(binnedBeam.beam.charge_per_macro)
                    
            
            zListSize = zReduced.size
            xMacros = np.empty(xReduced.size * yReduced.size * zReduced.size)
            yMacros = np.empty(xReduced.size * yReduced.size * zReduced.size)
            zMacros = np.empty(xReduced.size * yReduced.size * zReduced.size)
            chargeMacro = np.empty(xReduced.size * yReduced.size * zReduced.size)
            ListElement = 0
            Extras = 0
            for i in range(yReduced.size):
                for j in range(xReduced.size):
                    for k in range(zReduced.size):
                        #Remember to convert distances to cm
                        zMacros[ListElement + Extras] = zReduced[k] * 100
                        xMacros[ListElement + Extras] = (xReduced[j] * 100)
                        yMacros[ListElement + Extras] = (yReduced[i] * 100)
                        if k < zListSize:
                            chargeMacro[ListElement + Extras] = Charge3D[j][i][k] * calculation.chargePerMacro
                            ListElement+=1
                        else:
                            chargeMacro[ListElement + Extras] = 0 
                            Extras+=1
                            
            calculation._macros['x'] = xMacros
            calculation._macros['y'] = yMacros
            calculation._macros['z'] = zMacros
            calculation._macros['charge'] = chargeMacro
            
            #The DLW values
            calculation._InputParameters['Geometry'] = orientation
            calculation._InputParameters['MaxZ'] = max_z
            calculation._InputParameters['x0'] = width_cm * 0.5
            calculation._InputParameters['y0'] = offset * 1e2
            calculation._InputParameters['z0'] = 0
            calculation._InputParameters['Permitivity'] = permitivity
            calculation._InputParameters['Permeability'] = permeability
            calculation._InputParameters['a'] = a_gap
            calculation._InputParameters['delta'] = a_gap + delta
            calculation._InputParameters['w'] = width_cm
            calculation._InputParameters['sN'] = sN
            calculation._InputParameters['sI'] = sI
            calculation._InputParameters['ModePrec'] = ModePrec
            calculation._InputParameters['ModeAcc'] = ModeAcc
            calculation._InputParameters['ConvergenceCalculate'] = False
            
            #Do the calculation
            calculation.FieldCalculation()
            #Save the field then read the field
            calculation.write_Field_HDF5(SAVEFOLDER + "fieldfile_stage" + str(m) + "_dlw"+str(n)+".h5")
            DWA_Beam.GetFieldValues(SAVEFOLDER + "fieldfile_stage" + str(m) + "_dlw"+str(n)+".h5")
            
        #----------------------------------------------#
        #-----------4. Transport the beam--------------#
        #----------------------------------------------#
        DWA_Beam._DielectricParameter['L'] = dlw_length
        DWA_Beam.TransportThroughDLW(nSteps = 3)
        DWA_Beam.ClearLostParticles()
		#Save the beam at each step
		#You could always remove this
        DWA_Beam.beam.write_HDF5_beam_file(SAVEFOLDER + "beamfile_stage" + str(m) + "_dlw"+str(n)+".h5", overwrite_existing_file=True)
		
		
