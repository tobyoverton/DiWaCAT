# -*- coding: utf-8 -*-
"""
Class and Functions for DiWaCAT Fields
"""

import os
import numpy as np
import h5py
import scipy
import matplotlib.pyplot as plt
import math
from collections import OrderedDict

from Python_Tools.Modules import beam_tools as dbt


class DiWaCATOutput(object):
    def _init_(self):
        #Initialise the variables we need
        self.beam = {}
        self._DielectricParameter = {}
        self._FieldPoints = {}
            
    def _init_(self, filename = 'NONE'):
        #Initialise the variables we need
        self.beam = {}
        self._DielectricParameter = {}
        self._FieldPoints = {}
        if filename != 'NONE':
            self.beam = dbt.BeamFromFile()
        else:
            self.beam = self.ReadFromFile(filename)
            
            
    def reset_dicts(self):
        self.beam = {}
        self._DielectricParameter = {}
        self._FieldPoints = {}
        
    def CompleteSimulation(self, filename):
        self.ReadFromFile(filename)
        self.TransportThroughDLW(nSteps = 5)
        return self.beam
            
    def ReadFromFile(self, filename):
        self.reset_dicts()
        self.beam = dbt.BeamFromFile();
        self.beam.read_WakeCode_beam_file(filename);
        with h5py.File(filename, "r") as h5file:
            if h5file.get('DielectricParameters/Parameters') is not None:
                parameter_array = np.array(h5file.get('DielectricParameters/Parameters'))
                self._DielectricParameter['a'] = parameter_array[0][0] * 1e-6;
                self._DielectricParameter['delta'] = parameter_array[1][0] * 1e-6;
                #Cheating way of checking if the width is the old format (i.e. in m) - no width will be on micron scale so assume if it is must be in m
                if parameter_array[2][0] < 1:
                    self._DielectricParameter['w'] = parameter_array[2][0]
                else:
                    self._DielectricParameter['w'] = parameter_array[2][0]*1e-6;
                #Check if there actually is a width. If not set to zero. This is how we check if the structure is circular
                self._DielectricParameter['L'] = parameter_array[3][0];
                self._DielectricParameter['y0'] = parameter_array[5][0];
                self._DielectricParameter['x0'] = parameter_array[4][0];
                if (len(parameter_array)<10):
                    #Also check the field calculated isn't in the old format
                    print(h5file['DielectricParameters/columns'][()][2][0])
                    if h5file['DielectricParameters/columns'][()][2] in [b'w [micron]',b'w [m]',b'w [cm]']:
                        print('Old Field Type')
                    else:
                        self._DielectricParameter['w'] = 0
                        self._DielectricParameter['y0'] = parameter_array[4][0];
                        self._DielectricParameter['x0'] = parameter_array[3][0];
            if h5file.get('DielectricForces/ForceField') is not None:
                x, y, z, Fx, Fy, Fz = np.array(h5file.get('DielectricForces/ForceField')).transpose()
                self._FieldPoints['x'] = x
                self._FieldPoints['y'] = y
                self._FieldPoints['z'] = z
                self._FieldPoints['Fx'] = Fx
                self._FieldPoints['Fy'] = Fy
                self._FieldPoints['Fz'] = Fz
                self._FieldPoints['ForceArray'] = [];
                for i in range(len(Fx)):
                    self._FieldPoints['ForceArray'].append([Fx[i], Fy[i], Fz[i]])
                    
    def SetDielectricLength(self, length):
        self._DielectricParameter['L'] = length
        
    def SetBeam(self, beam):
        self.beam = beam;
    
    def GetFieldValues(self, filename):
        self._DielectricParameter = {};
        self._FieldPoints = {}
        with h5py.File(filename, "r") as h5file:
            if h5file.get('DielectricParameters/Parameters') is not None:
                parameter_array = np.array(h5file.get('DielectricParameters/Parameters'))
                self._DielectricParameter['a'] = parameter_array[0][0] * 1e-6;
                self._DielectricParameter['delta'] = parameter_array[1][0] * 1e-6;
                #Cheating way of checking if the width is the old format (i.e. in m) - no width will be on micron scale so assume if it is must be in m
                if parameter_array[2][0] < 1:
                    self._DielectricParameter['w'] = parameter_array[2][0]
                else:
                    self._DielectricParameter['w'] = parameter_array[2][0]*1e-6;
                #Check if there actually is a width. If not set to zero. This is how we check if the structure is circular
                self._DielectricParameter['L'] = parameter_array[3][0];
                self._DielectricParameter['y0'] = parameter_array[5][0];
                self._DielectricParameter['x0'] = parameter_array[4][0];
                if (len(parameter_array)<10):
                    #Also check the field calculated isn't in the old format
                    print(h5file['DielectricParameters/columns'][()][2][0])
                    if h5file['DielectricParameters/columns'][()][2] in [b'w [micron]',b'w [m]',b'w [cm]']:
                        print('Old Field Type')
                    else:
                        self._DielectricParameter['w'] = 0
                        self._DielectricParameter['y0'] = parameter_array[4][0];
                        self._DielectricParameter['x0'] = parameter_array[3][0];
            if h5file.get('DielectricForces/ForceField') is not None:
                x, y, z, Fx, Fy, Fz = np.array(h5file.get('DielectricForces/ForceField')).transpose()
                self._FieldPoints['x'] = x
                self._FieldPoints['y'] = y
                self._FieldPoints['z'] = z
                self._FieldPoints['Fx'] = Fx
                self._FieldPoints['Fy'] = Fy
                self._FieldPoints['Fz'] = Fz
                self._FieldPoints['ForceArray'] = [];
                for i in range(len(Fx)):
                    self._FieldPoints['ForceArray'].append([Fx[i], Fy[i], Fz[i]])
    
    def ReturnFieldInterpolated(self,SamplePoint):
        xArray = self._FieldPoints['x']
        yArray = self._FieldPoints['y']
        zArray = self._FieldPoints['z']
        Fx_values = self._FieldPoints['Fx']
        Fy_values = self._FieldPoints['Fy']
        Ez_values = self._FieldPoints['Fz']            
        xs = list(OrderedDict.fromkeys(xArray))
        ys = list(OrderedDict.fromkeys(yArray))
        zs = list(OrderedDict.fromkeys(zArray))
        Fx_grid = np.array(Fx_values).reshape((len(ys), len(xs), len(zs)),order='C').transpose(1,0,2)
        Fy_grid = np.array(Fy_values).reshape((len(ys), len(xs), len(zs)),order='C').transpose(1,0,2)
        Ez_grid = np.array(Ez_values).reshape((len(ys), len(xs), len(zs)),order='C').transpose(1,0,2)
        
        interp_Fx = scipy.interpolate.RegularGridInterpolator((xs, ys, zs), Fx_grid, bounds_error=False, fill_value=None)
        interp_Fy = scipy.interpolate.RegularGridInterpolator((xs, ys, zs), Fy_grid, bounds_error=False, fill_value=None)
        interp_Ez = scipy.interpolate.RegularGridInterpolator((xs, ys, zs), Ez_grid, bounds_error=False, fill_value=None)
        
        return interp_Fx(SamplePoint),interp_Fy(SamplePoint),interp_Ez(SamplePoint)
    
    def ReturnFieldFunctions(self):
        #Returns the functions not actual values - for values use ReturnFieldInterpolated
        xArray = self._FieldPoints['x']
        yArray = self._FieldPoints['y']
        zArray = self._FieldPoints['z']
        Fx_values = self._FieldPoints['Fx']
        Fy_values = self._FieldPoints['Fy']
        Ez_values = self._FieldPoints['Fz']            
        xs = list(OrderedDict.fromkeys(xArray))
        ys = list(OrderedDict.fromkeys(yArray))
        zs = list(OrderedDict.fromkeys(zArray))
        Fx_grid = np.array(Fx_values).reshape((len(ys), len(xs), len(zs)),order='C').transpose(1,0,2)
        Fy_grid = np.array(Fy_values).reshape((len(ys), len(xs), len(zs)),order='C').transpose(1,0,2)
        Ez_grid = np.array(Ez_values).reshape((len(ys), len(xs), len(zs)),order='C').transpose(1,0,2)
        return [scipy.interpolate.RegularGridInterpolator((xs, ys, zs), Fx_grid, bounds_error=False, fill_value=None), scipy.interpolate.RegularGridInterpolator((xs, ys, zs), Fy_grid, bounds_error=False, fill_value=None), scipy.interpolate.RegularGridInterpolator((xs, ys, zs), Ez_grid, bounds_error=False, fill_value=None)]
    
    
    def TransportThroughDLW(self, nSteps=2, remove_losses = False):
        #Simple particle pusher through the DLW
        #Linear process with nSteps through the structure - applying the field at the end of each step
        #Can choose whether or not remove particles collamated by the DLW walls
        steplength = self._DielectricParameter['L'] * 1/nSteps
        for i in range(nSteps):
            #Apply in order - the kick has already been applied to momenta so just need to drift
            self.driftBeam(steplength * 0.5)
            #Save positions at the half-step point - this is where the field will be interpolated
            #The half-step keeps it a Boris pusher and means we're calculating field at the 'average' position in the step
            savedpositions = np.stack((self.beam._beam['x'],self.beam._beam['y'],self.beam._beam['z']),axis=-1)
            self.driftBeam(steplength*0.5);
            #Apply kick to the momenta of macroparticles
            KicksApplied = self.ReturnFieldInterpolated(savedpositions)
            for i in range(len(self.beam._beam['x'])):
                self.beam._beam['px'][i] += steplength * KicksApplied[0][i] * self.beam.q_over_c
                self.beam._beam['py'][i] += steplength * KicksApplied[1][i] * self.beam.q_over_c
                self.beam._beam['pz'][i] += steplength * KicksApplied[2][i] * self.beam.q_over_c
                
            #Remove particles collimated by the DLW
            if remove_losses == True:
                self.ClearLostParticles()
    
    def driftBeam(self, L):
        #Copy of the beam_tools driftBeam function
        #Don't move on the longitudinal coordinate - causes issues with the interpolation in our case
        self.beam._beam['x'] += L * self.beam.xp
        self.beam._beam['y'] += L * self.beam.yp
            
    def ClearLostParticles(self):
        LostParticlesY = [n for n,i in enumerate(self.beam._beam['y']) if abs(i+self._DielectricParameter['y0']) > self._DielectricParameter['a']]
        self.beam._beam['charge'][LostParticlesY] = 0; 
        LostParticlesX = [n for n,i in enumerate(self.beam._beam['x']) if abs(i+self._DielectricParameter['x0']) > self._DielectricParameter['w']]
        self.beam._beam['charge'][LostParticlesX] = 0; 
    
        charge = self.beam._beam['charge'];
        if np.any(charge == 0):
            macro_select = np.nonzero(charge)
            print("Beam losses detected. Charge Lost = ", (len(self.beam._beam['x']) - len(macro_select))*self.beam.charge_per_macro, " C")
        else:
            macro_select = ...

        self.beam._beam['x'] = self.beam._beam['x'][macro_select]
        self.beam._beam['y'] = self.beam._beam['y'][macro_select]
        self.beam._beam['z'] = self.beam._beam['z'][macro_select]
        self.beam._beam['px'] = self.beam._beam['px'][macro_select]
        self.beam._beam['py'] = self.beam._beam['py'][macro_select]
        self.beam._beam['pz'] = self.beam._beam['pz'][macro_select]
        self.beam._beam['t'] = self.beam._beam['t'][macro_select]
        self.beam._beam['charge'] = self.beam._beam['charge'][macro_select]
        
    def ChangeOrientation(self):
        #Note this assumes the fields can just be swapped (i.e. needs a symmetric beam)
        fx_temp = self._FieldPoints['Fx'];
        self._FieldPoints['Fx'] = self._FieldPoints['Fy'];
        self._FieldPoints['Fy'] = fx_temp;
        forcearray = [];
        for i in range(len(fx_temp)):
            forcearray.append([self._FieldPoints['Fx'][i], self._FieldPoints['Fy'][i], self._FieldPoints['Fx'][i]])
        self._FieldPoints['ForceArray'] = forcearray;
