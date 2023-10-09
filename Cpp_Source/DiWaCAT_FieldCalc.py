# -*- coding: utf-8 -*-
"""
Field Calculator Functions
T Overton
4/10/23

Reads in a meshed beam file and calculates fields for circular and planar DLW
"""

import numpy as np
import h5py
import os, sys
import scipy.constants as const
import math
import warnings
sys.stdout.flush()

from Python_Tools.Modules import beam_tools as dbt
import Cpp_Source.DiWaCAT_library as DiWaCAT

def UnitConversion(prefix):
    if prefix == 'unit' or prefix == 'Unit':
        return 1
    else:
        return 1/getattr(const, prefix)
    

'''----------------------------------------


CREATE CLASS TO HOLD THE FIELD CALCULATION FUNCTIONS/INPUTS

-----------------------------------------'''

class DiWaCAT_FieldCalc(dbt.BeamFromFile):
    def __init__(self):
        super().__init__()
        self.resetVals()
    def resetVals(self):
        self.resetDicts()
        self._macros= {}
        self._Field = {};
        self._InputParameters = {}
        self._InputParameters['x0'] = 0
        self._InputParameters['y0'] = 0
        self._InputParameters['MaxZ'] = 0
    
    def BeamFile(self, FileName):
        self.BeamFile = FileName;
    @property
    def xMacros(self):
        return self._macros['x'] + self.x0
    @property
    def xMaxMacro(self):
        return np.max(np.abs(self.xMacros))
    @property
    def yMacros(self):
        return self._macros['y'] + self.y0
    @property
    def yMaxMacro(self):
        return np.max(np.abs(self.yMacros))
    @property
    def zMacros(self):
        return self._macros['z'] + self.z0
    @property
    def zMacroInterval(self):
        return self.zMacros[1] - self.zMacros[0]
    @property
    def rMacros(self):
        return np.sqrt(self.xMacros**2 + self.yMacros**2)   
    @property
    def thetaMacros(self):
        return np.arctan2(self.xMacros, self.yMacros)
    @property
    def MacroCharge(self):
        return self._macros['charge']
    @property
    def chargePerMacro(self):
        return self._macros['ChargePerMacro']
    @property
    def x0(self):
        return self._InputParameters['x0']
    @property
    def y0(self):
        return self._InputParameters['y0']
    @property
    def r0(self):
        return math.sqrt(self.x0**2 + self.y0**2)
    @property
    def theta0(self):
        return math.atan2(self.y0, self.x0)
    @property
    def z0(self):
        return self._InputParameters['z0']
    @property
    def Ep(self):
        return self._InputParameters['Permitivity']
    @property
    def Mu(self):
        return self._InputParameters['Permeability']
    @property
    def b(self):
        return self._InputParameters['a']
    @property
    def delta(self):
        return self._InputParameters['delta']
    @property
    def c(self):
        return self.b + self.delta
    @property
    def w(self):
        return self._InputParameters['w']
    @property
    def sN(self):
        return self._InputParameters['sN']
    @property
    def sI(self):
        return self._InputParameters['sI']
    @property
    def nR(self):
        return self._InputParameters['nR']
    @property
    def nTheta(self):
        return self._InputParameters['nT']
    @property
    def ModePrec(self):
        return self._InputParameters['ModePrec']
    @property
    def ModeAcc(self):
        return self._InputParameters['ModeAcc']
    @property
    def Convergence(self):
        return self._InputParameters['ConvergenceCalculate']
    @property
    def MaxZ(self):
        return self._InputParameters['MaxZ']
    @property
    def Geometry(self):
        return self._InputParameters['Geometry']
    @property
    def Ez(self):
        return self._Field['Ez']
    @property
    def Fx(self):
        return self._Field['Fx']
    @property
    def Fy(self):
        return self._Field['Fy']
    
    def ReadMeshValues(self, BeamFile):
        self.read_WakeCode_beam_file(BeamFile)
        with h5py.File(BeamFile, "r") as h5file:
            xReduced = []
            yReduced = []
            zReduced = []
            Charge3D = []
            #Use a flag to ignore deprecation warnings - caused by the h5py library so out of our control
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=DeprecationWarning)
                if h5file.get("MeshValues/macrosAtPoint") is not None:
                    Charge3D = np.array(h5file.get('MeshValues/macrosAtPoint'))
                if h5file.get("MeshValues/xPoints") is not None:
                    xReduced = np.array(h5file.get('MeshValues/xPoints')) 
                if h5file.get("MeshValues/yPoints") is not None:
                    yReduced = np.array(h5file.get('MeshValues/yPoints')) 
                if h5file.get("MeshValues/zPoints") is not None:
                    zReduced = np.array(h5file.get('MeshValues/zPoints'))
                if h5file.get("MeshParameters/chargePerMacro") is not None:
                    self._macros['ChargePerMacro'] = np.array(h5file.get("MeshParameters/chargePerMacro"))
                    
            if 'MaxZ' not in self._InputParameters:
                self._InputParameters['MaxZ'] = 0
            maxZValue = np.max(zReduced)
            zListSize = zReduced.size
            if maxZValue < self.MaxZ:
                zInterval = zReduced[1] - zReduced[0]
                ZAdditions = math.ceil((self.MaxZ - maxZValue)*1/zInterval)
                zReduced.resize(zListSize + ZAdditions)
                for i in range(ZAdditions):
                    zReduced[zListSize + i] = zReduced[zListSize + i - 1] + zInterval
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
                            chargeMacro[ListElement + Extras] = Charge3D[i][j][k] * self.chargePerMacro
                            ListElement+=1
                        else:
                            chargeMacro[ListElement + Extras] = 0 
                            Extras+=1
                            
            self._macros['x'] = xMacros
            self._macros['y'] = yMacros
            self._macros['z'] = zMacros
            self._macros['charge'] = chargeMacro
            
            
    def InputPlanarParameters(self, InputList):
        self._InputParameters['Geometry'] = InputList[0] #Either 'h' or 'v'
        self._InputParameters['MaxZ'] = InputList[12]
        
        #Read in beam and structure parameters
        #Functions are all in cm so make sure to convert all the inputs
        self._InputParameters['x0'] = InputList[1] * 100
        self._InputParameters['y0'] = InputList[2] * 100
        self._InputParameters['z0'] = 0
        self._InputParameters['Permitivity'] = InputList[3] # Dielectric Permitivity
        self._InputParameters['Permeability'] = InputList[4] # Dielectric Permeability
        self._InputParameters['a'] = InputList[5] * 100#DLW half-gap
        self._InputParameters['delta'] = InputList[6]* 100 #Half-gap + Dielectric Thickness
        self._InputParameters['w'] = InputList[7]* 100
        self._InputParameters['sN'] = InputList[8] #Number of x-modes
        self._InputParameters['sI'] = InputList[9] # Number of y-modes
        self._InputParameters['ModePrec'] = InputList[10] # Precision of mode root finder
        self._InputParameters['ModeAcc'] = InputList[11] * 1/100 #Amount of wake potential in the final 5 modes
        if (InputList[13] == 'f'): self._InputParameters['ConvergenceCalculate'] = False
        else: self._InputParameters['ConvergenceCalculate'] = True
        
        #Add half the width to the x offset and swap positions x/y offsets if vertical orientation
        #y0 from input is always towards the dielectric
        self._InputParameters['x0'] += self.w/2
        if (self.Geometry == 'v'):
            self._InputParameters['y0'] = self._InputParameters['x0']
            self._InputParameters['x0'] = InputList[2] * 100

    def InputCircularParameters(self,InputList):
        self._InputParameters['Geometry'] = 'c'
        self._InputParameters['MaxZ'] = InputList[10]
        #Read in beam/structure parameters
        self._InputParameters['x0'] = InputList[0] * 100
        self._InputParameters['y0'] = InputList[1] * 100
        self._InputParameters['z0'] = 0
        self._InputParameters['Permitivity'] = InputList[2] # Dielectric Permitivity
        self._InputParameters['Permeability'] = InputList[3] # Dielectric Permeability
        self._InputParameters['a'] = InputList[4] * 100#DLW half-gap
        self._InputParameters['delta'] = InputList[5]* 100 #Half-gap + Dielectric Thickness
        self._InputParameters['nR'] = InputList[6];     #Number of frequcy modes in r
        self._InputParameters['nT'] = InputList[7];     #Number of frequency modes in theta
        self._InputParameters['ModePrec'] = InputList[8] # Precision of mode root finder
        self._InputParameters['ModeAcc'] = InputList[9] * 1/100 #Amount of wake potential in the final 5 modes
        if (InputList[11] == 'f'): self._InputParameters['ConvergenceCalculate'] = False
        else: self._InputParameters['ConvergenceCalculate'] = True
        

    def FieldCalculation(self):
        #InputList is the list of variables for  calculation
        if self.Geometry == 'h' or self.Geometry=='v':
            print("Planar Field Calculation")
            
            #Calculate beta squared
            B2 = np.mean((self.Bz)**2)
            "---------------- Calculate Eigenvalues of Modes ----------------"
            #Calculate mode Convergence
            zeroesES = []
            zeroesEA = []
            zeroesHS =[]
            zeroesHA = []
            if self.Convergence == True:
                yModesAccurate = False
                xModesAccurate = False
                #Needs to be if either of the two are false since need convergence in both
                while yModesAccurate == False or xModesAccurate == False:
                    dFzFullModes = dFyFullModes = dFxFullModes = 0
                    dFzFewerModesX = dFyFewerModesX = dFxFewerModesX = 0
                    dFzFewerModesY = dFyFewerModesY = dFxFewerModesY = 0
                    zeroesES, zeroesEA, zeroesHS, zeroesHA = DiWaCAT.EigenvaluesCalculator(self.sI, self.sN, zeroesES, zeroesEA, zeroesHS, zeroesHA, self.Ep, self.Mu, B2, self.b, self.c, self.w, self.ModePrec)
                    if (np.isnan(zeroesES).any() or np.isnan(zeroesEA).any() or np.isnan(zeroesHS).any() or np.isnan(zeroesHA).any()):
                        zeroesES, zeroesEA, zeroesHS, zeroesHA = DiWaCAT.EigenvaluesCalculator(self.sI, self.sN, zeroesES, zeroesEA, zeroesHS, zeroesHA, self.Ep, self.Mu, B2, self.b, self.c, self.w, self.ModePrec)   
                    if self.Geometry == 'h':
                        ymax = self.yMaxMacro
                        if ymax > self.b: ymax = self.b - 0.0005
                        dFzFullModes, dFyFullModes, dFxFullModes = DiWaCAT.CalcWakeElement(dFzFullModes, dFxFullModes, dFyFullModes, self.xMaxMacro, self.x0, ymax, self.y0, self.zMacroInterval, 0, zeroesES, zeroesEA, zeroesHS, zeroesHA, self.sN, self.sI, self.Ep, self.Mu, B2, self.b, self.c, self.w);
                        dFzFewerModesX, dFyFewerModesX, dFxFewerModesX = DiWaCAT.CalcWakeElement(dFzFullModes, dFxFullModes, dFyFullModes, self.xMaxMacro, self.x0, ymax, self.y0, self.zMacroInterval, 0, zeroesES, zeroesEA, zeroesHS, zeroesHA, self.sN - 5, self.sI, self.Ep, self.Mu, B2, self.b, self.c, self.w);
                        dFzFewerModesY, dFyFewerModesY, dFxFewerModesY = DiWaCAT.CalcWakeElement(dFzFullModes, dFxFullModes, dFyFullModes, self.xMaxMacro, self.x0, ymax, self.y0, self.zMacroInterval, 0, zeroesES, zeroesEA, zeroesHS, zeroesHA, self.sN, self.sI - 5, self.Ep, self.Mu, B2, self.b, self.c, self.w);
                    elif self.Geometry == 'v':
                        xmax = self.xMaxMacro
                        if xmax > self.b: xmax = self.b - 0.0005
                        dFzFullModes, dFyFullModes, dFxFullModes = DiWaCAT.CalcWakeElement(dFzFullModes, dFxFullModes, dFyFullModes, self.yMaxMacro, self.y0, xmax, self.x0, self.zMacroInterval, 0, zeroesES, zeroesEA, zeroesHS, zeroesHA, self.sN, self.sI, self.Ep, self.Mu, B2, self.b, self.c, self.w);
                        dFzFewerModesX, dFyFewerModesX, dFxFewerModesX = DiWaCAT.CalcWakeElement(dFzFullModes, dFxFullModes, dFyFullModes, self.yMaxMacro, self.y0, xmax, self.x0, self.zMacroInterval, 0, zeroesES, zeroesEA, zeroesHS, zeroesHA, self.sN - 5, self.sI, self.Ep, self.Mu, B2, self.b, self.c, self.w);
                        dFzFewerModesY, dFyFewerModesY, dFxFewerModesY = DiWaCAT.CalcWakeElement(dFzFullModes, dFxFullModes, dFyFullModes, self.yMaxMacro, self.y0, xmax, self.x0, self.zMacroInterval, 0, zeroesES, zeroesEA, zeroesHS, zeroesHA, self.sN, self.sI - 5, self.Ep, self.Mu, B2, self.b, self.c, self.w);
                    
                    if (((dFzFullModes - dFzFewerModesX) * 1/dFzFullModes < self.ModeAcc) and ((dFyFullModes - dFyFewerModesX) * 1/dFyFullModes < self.ModeAcc) and ((dFxFullModes - dFxFewerModesX) * 1/dFxFullModes < self.ModeAcc)):
                        xModesAccurate = True
                    else: self._InputParameters['sN'] +=2
                    if (((dFzFullModes - dFzFewerModesY) * 1/dFzFullModes < self.ModeAcc) and ((dFyFullModes - dFyFewerModesY) * 1/dFyFullModes < self.ModeAcc) and ((dFxFullModes - dFxFewerModesY) * 1/dFxFullModes < self.ModeAcc)):
                        yModesAccurate = True
                    else: self._InputParameters['sI'] +=2
            
            print("Modes Used: Nx = ", self.sN, ", Ny = ", self.sI)
            zeroesES, zeroesEA, zeroesHS, zeroesHA = DiWaCAT.EigenvaluesCalculator(self.sI, self.sN, zeroesES, zeroesEA, zeroesHS, zeroesHA, self.Ep, self.Mu, B2, self.b, self.c, self.w, self.ModePrec)
            
            
            "------------------ SET UP 3D FORCE MATRICES  -------------------"
            nParticle = self.zMacros.size
            Fx = np.empty(nParticle)
            Fy = np.empty(nParticle)
            Ez = np.empty(nParticle)
            "-------------------- RUN CODE FOR FIELD VALUES ------------------"
            
            if self.Geometry == 'h':
                Fx,Fy,Ez = DiWaCAT.TotalForceMeshHDF5(Ez, Fx, Fy, nParticle, 0, self.x0, 0, self.y0, 0, np.min(self.zMacros), np.max(self.zMacros), self.z0, self.xMacros, self.yMacros, self.zMacros, zeroesES, zeroesEA, zeroesHS, zeroesHA, self.sN, self.sI, self.Ep, self.Mu, B2, self.b, self.c, self.w, self.MacroCharge)
            elif self.Geometry == 'v':
                Fx,Fy,Ez = DiWaCAT.TotalForceMeshHDF5VerticalPlate(Ez, Fx, Fy, nParticle, 0, self.x0, 0, self.y0, 0, np.min(self.zMacros), np.max(self.zMacros), self.z0, self.xMacros, self.yMacros, self.zMacros, zeroesES, zeroesEA, zeroesHS, zeroesHA, self.sN, self.sI, self.Ep, self.Mu, B2, self.b, self.c, self.w, self.MacroCharge)
            #Store the calculated fields
            self._Field['Fx'] = Fx
            self._Field['Fy'] = Fy
            self._Field['Ez'] = Ez
            
        elif self.Geometry == 'c':
            print('Circular Field Calculation')
            #REMEMBER CIRCULAR IS IN METRES NOT CM
            
            '''----------- Convert Macro Positions to Cylindrical Coordinates ---------------'''
            #Get maximum radial position
            rMax = np.max(self.rMacros)*0.01
            
            '''-------------- Set up the 3D Force Vectors -------------------'''
            nParticle = self.zMacros.size
            Fx = np.empty(nParticle)
            Fy = np.empty(nParticle)
            Ez = np.empty(nParticle)
            Fr = np.empty(nParticle)
            Ftheta = np.empty(nParticle)
            
            '''------------------- Calculate Mode Amplitudes and Convergence ---------------'''
            ModeAmplitude = []
            WaveVectors = []
            if self.Convergence == True:
                WaveVectors, ModeAmplitude = DiWaCAT.ModeConvergence(rMax, self.b * 0.01, self.c*0.01, self.Mu, self.Ep, self.ModePrec, self.ModeAcc)
            else:
                WaveVectors, ModeAmplitude = DiWaCAT.FindModes(self.b * 0.01, self.c * 0.01, self.Mu, self.Ep, self.nTheta, self.nR, self.ModePrec, False)
            '''---------------- Calculate Forces --------------'''
            Ez, Fx, Fy = DiWaCAT.TotalForceMeshCircular(Ez, Fr, Ftheta, self.rMacros, self.thetaMacros, self.zMacros, self.MacroCharge, ModeAmplitude, WaveVectors, self.b * 0.01, self.c * 0.01, self.Mu, self.Ep)
            
            '''-------------- Convert Cylidrical Forces back to Cartesian -----------------'''
            #Fy = (Fr * np.cos(self.thetaMacros)) - (Ftheta * np.sin(self.thetaMacros))
            #Fx = (Fr * np.sin(self.thetaMacros)) - (Ftheta * np.cos(self.thetaMacros))
            
            #Store the calculated fields
            self._Field['Fx'] = Fx
            self._Field['Fy'] = Fy
            self._Field['Ez'] = Ez
            
        else:
            print('No DLW Geometry Selected')
            
    def write_Field_HDF5(self, filename, overwrite_existing_file=False, centered=False, mass=const.m_e,
                             sourcefilename=None, pos=None, rotation=None, longitudinal_reference='t', xoffset=0,
                             yoffset=0, zoffset=0, toffset=0, ):
        #Save field to HDF5 file
        #Include the beam and mesh
        #Method copied from beam_tools for beam and macros
        with h5py.File(filename, "w") as f:
            '''--------------WRITE BEAM-------------------------'''
            beamgrp = f.create_group("beam")

            if 'reference_particle' in self._beam:
                beamgrp['reference_particle'] = self._beam['reference_particle']
            if 'status' in self._beam:
                beamgrp['status'] = self._beam['status']
            beamgrp['longitudinal_reference'] = longitudinal_reference

            # print('hdf5 write cathode', cathode, np.array(beamgrp['cathode']))
            if len(self._beam['charge']) == len(self.x):
                chargevector = self._beam['charge']
            else:
                chargevector = np.full(len(self.x), self.charge / len(self.x))
            array = np.array(
                [self.x + xoffset, self.y + yoffset, self.z + zoffset, self.cpx, self.cpy, self.cpz, self.t + toffset,
                 chargevector]).transpose()
            beamgrp['columns'] = np.array(['x', 'y', 'z', 'cpx', 'cpy', 'cpz', 't', 'q'], dtype='S')
            beamgrp['units'] = np.array(['m', 'm', 'm', 'eV', 'eV', 'eV', 's', 'e'], dtype='S')
            beamgrp.create_dataset("beam", data=array)
            
            '''---------------WRITE MACROS-----------------------'''
            #Need to work out exactly how to do
            
            
            '''------------- WRITE FIELD VALUES ----------------'''
            fieldgrp = f.create_group("DielectricForces")
            forcearray = np.array(
                [1e-2*(self.xMacros - self.x0), 1e-2*(self.yMacros - self.y0), 1e-2*(self.zMacros), self.Fx, self.Fy, self.Ez]).transpose()
            fieldgrp.create_dataset("ForceField", data = forcearray)
            fieldgrp['columns'] = np.array(['x [m]', 'y [m]', 'z [m]', 'Fx [eV]', 'Fy [eV]', 'Ez [eV]'], dtype='S')
            parametergrp = f.create_group("DielectricParameters")
            if self.Geometry == 'c':
                parameters = np.array([[self.b*1e4, self.delta*1e4, self.x0*0.01, self.y0*0.01, self.Ep, self.Mu, self.nR, self.nTheta]]).transpose()
                parametergrp.create_dataset("Parameters", data = parameters)
                parametergrp['columns'] = np.array(['a [micron]', 'delta [micron]', 'x0 [m]', 'y0 [m]', 'Permitivity', 'Permeability', 'nR Modes', 'nTheta Modes'], dtype='S')            
            else:
                parameters = np.array([[self.b*1e4, self.delta*1e4, self.w*1e4, self.x0*0.01, self.y0*0.01, self.Ep, self.Mu, self.sN, self.sI]]).transpose()
                parametergrp.create_dataset("Parameters", data = parameters)
                parametergrp['columns'] = np.array(['a [micron]', 'delta [micron]', 'w [micron]', 'x0 [m]', 'y0 [m]', 'Permitivity', 'Permeability', 'nX Modes', 'nY Modes'], dtype='S')
            

        return os.path.isfile(filename)  # because we wrote a file!
        
'''
            
a = DiWaCAT_FieldCalc()
a.InputPlanarParameters(['h', 0, 0, 4, 1, 1e-3,200e-6, 2e-2,5, 5, 0.001, 0.1, 0, 't'])
a.ReadMeshValues('./../../Streaker_Work/Slice_Emittance/Sim_DoubleGaus/DoubleGaus_In.h5')
a.read_WakeCode_beam_file('./../../Streaker_Work/Slice_Emittance/Sim_DoubleGaus/DoubleGaus_In.h5')
a.FieldCalculation()

a.write_Field_HDF5('test.h5')
'''