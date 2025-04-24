'''
Ocelot example beam tracking

In this example:
    *Read in a DiWaCAT beam and convert to an ocelot ParticleArray
    *Define a beamline as a quadrupole triplet followed by two DLWs
    *H+V DLW structures with an offset in the horizontal DLW - i.e. streaker
    *Saves the beam with DiWaCAT formatting
'''
import sys
#If ocelot is not installed - need to uncomment this line and change to direct to ocelot
#sys.path.append("PATH/PATH-TO-OCELOT")
from Python_Tools.Modules.diwacat_ocelot import DielectricWakefieldProc
from ocelot import * 
import numpy as np
import copy
import matplotlib.pyplot as plt
#This line assumes the base DiWaCAT library is added to system path - if not direct like ocelot (line 16)
import Python_Tools.Modules.beam_tools as dbt

#Read in a DiWaCAT beam
beam = dbt.BeamFromFile()
beam.read_WakeCode_beam_file(".\\..\\Example_Scripts\\Example_DiWaCAT_Beam.h5")
beamtrack = ParticleArray()
'''
First need to convert a DiWaCAT beam to a ParticleArray beam
'''
#The z positions may need flipping - DiWaCAT has z and t with the wrong sign. If it does multiply by -1
#Should be fine as long as other physics processes aren't involved - the DWA field calculation uses the DiWACAT convention
zFactor = 1

nMacros = beam.nMacros
beamtrack.__dict__['rparticles'] = np.zeros((6, nMacros))
beamtrack.rparticles[0] = beam.x  - np.mean(beam.x)
beamtrack.rparticles[1] = beam.xp
beamtrack.rparticles[2] = beam.y - np.mean(beam.y)
beamtrack.rparticles[3] = beam.yp
beamtrack.rparticles[4] = zFactor * beam.zn 
beamtrack.rparticles[5] = beam.pz / np.mean(beam.pz)
beamtrack.__dict__['E'] = beam.Mcpz * 1e-9
beamtrack.__dict__['s'] = 0
beamtrack.__dict__['q_array'] = beam.charge

'''
Now lets set up the beamline
For this example we'll do a quadrupole triplet, followed by a H DLW, V DLW and then drift to a screen
'''

# Quadrupoles
Quad1 = Quadrupole(l=0.24, k1=0.2237922517, eid='Quad1')
Quad2 = Quadrupole(l=0.24, k1=-0.2040914523, eid='Quad2')
Quad3 = Quadrupole(l=0.24, k1=0.24059635850000002, eid='Quad3')

#Drifts between the quadrupoles
drift1 = Drift(l = 0.1, eid="drift1")
drift2 = Drift(l = 0.5, eid="drift2")
drift3 = Drift(l = 0.5, eid="drift3")
drift4 = Drift(l = 0.1, eid="drift4")

#First define the positions within the DLWs and the screen
ws1_start = Marker(eid="ws1_start")
ws1_stop = Marker(eid="ws1_stop")
ws1_center = Marker(eid="ws1_center")
ws2_start = Marker(eid="ws2_start")
ws2_stop = Marker(eid="ws2_stop")
ws2_center = Marker(eid="ws2_center")
screen = Marker(eid="OTRB.2560.T3")
#The cell elements for the DLWs and then drift to the screen
wakefield_structureh = (ws1_start, Drift(l=0.1), ws1_center, Drift(l=0.1), ws1_stop)
wakefield_structurev = (ws2_start, Drift(l=0.1), ws2_center, Drift(l=0.1), ws2_stop)
drift_to_screen = (Drift(l=1, eid='DriftToScreen'), screen)


#Define the total beamline
cell = (drift1, Quad1, drift2, Quad2, drift3, Quad3, drift4, wakefield_structureh, wakefield_structurev, drift_to_screen)


#In this example the first DLW is streaking the beam (a = 2mm, offset = 1.7mm) and second is cancelling some of the quadrupole wakefield (a = 650um, offset = 0)
#In both the dielectric thickness is 200um (permitivity = 3.75). See the DielectricWakefieldProc class for how to define the other parameters
ws1 = DielectricWakefieldProc(2e-3, 200e-6, 1.7e-3)
ws2 = DielectricWakefieldProc(650e-6, 200e-6, 0, orientation = 'v') #Default orientation is 'h' so need to state if vertical


lattice = MagneticLattice(cell, stop=screen)
navi = Navigator(lattice)


navi.add_physics_proc(ws1, ws1_start, ws1_stop)
navi.add_physics_proc(ws2, ws2_start, ws2_stop)

_, trackBeam = track(lattice, p_array=beamtrack, navi=navi, calc_tws=False, print_progress=False)

#Here we can output our beam, plot things, whatever we like
#For this example lets get back to a DiWaCAT beam
outputbeam = dbt.GeneralBeam()
particles = trackBeam.rparticles
outputbeam._beam['x'] = particles[0]
outputbeam._beam['y'] = particles[2]
outputbeam._beam['z'] = zFactor * particles[4]
energycentral = trackBeam.E
outputbeam._beam['px'] = energycentral * (1/const.c) * particles[1]
outputbeam._beam['py'] = energycentral * (1/const.c) * particles[3]
outputbeam._beam['pz'] = energycentral * (1/const.c) * particles[5]
outputbeam._beam['t'] = outputbeam.zn * -1 * 1/const.c
outputbeam._beam['macro_type'] = 'variable_weight'
outputbeam._beam['code'] = 'ocelot'
outputbeam._beam['charge'] = beamtrack.q_array

outputbeam.write_HDF5_beam_file('testoutput.h5', overwrite_existing_file = True)