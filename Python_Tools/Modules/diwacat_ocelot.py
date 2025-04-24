'''
Ocelot_StreakerWakefield

Implement DiWaCAT wakefields into Ocelot
Allowing for true generalisation of the phase space reconstruction using ocelot simulations
Use DiWaCAT to get the wake potential and then initialise a PhysProc object - this has function apply to do the tracking
'''


import sys, os
import math
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
#If ocelot is not installed - need to uncomment this line and change to direct to ocelot
#sys.path.append("PATH/PATH-TO-OCELOT")
from ocelot import *
from ocelot.gui import *
from ocelot.utils.dechirper_recon import *
from ocelot.cpbd.physics_proc import PhysProc
import scipy.constants as const
#This line assumes the base DiWaCAT library is added to system path - if not direct like ocelot (line 16)
import Cpp_Source.DiWaCAT_library as DiWaCAT


class DielectricWakefieldProc(PhysProc):
    '''
    Tracking through a DWA structures including
    '''
    def __init__(self, halfGap, dielectricThickness, offset, orientation = 'h', dielectricPermitivity = 3.75, dielectricPermeability = 1, meanEnergy_MeV = 250, length = 0.2, width = 1e-2, step = 1, maxZPosition = 5e-12*const.c):
        PhysProc.__init__(self)
        self.step = step
        self._parameters = {}
        self._calculation = {}
        self._parameters['a'] = halfGap
        self._parameters['delta'] = dielectricThickness
        self._parameters['width'] = width
        self._parameters['Permitivity'] = dielectricPermitivity
        self._parameters['Permeability'] = dielectricPermeability
        self._parameters['BeamEnergy'] = meanEnergy_MeV
        self._parameters['Offset'] = offset
        self._parameters['TotalLength'] = length
        if orientation == 'v' or orientation == 'V' or orientation == 'vertical' or orientation == 'Vertical':
            self._parameters['Horizontal'] = False
        else:
            self._parameters['Horizontal'] = True
        self.wake_calculated = False
        self._calculation['LongitudinalPositions'] = np.linspace(0, maxZPosition, 10000);
        
    '''
    Inputs for the parameters will be in m but wake potential requires cm so multiply by 1e2
    '''
    @property
    def b(self):
        return self._parameters['a'] * 1e2
    @property
    def c(self):
        return (self._parameters['a'] + self._parameters['delta'])*1e2
    @property
    def w(self):
        return self._parameters['width']*1e2
    @property
    def Energy(self):
        #Energy must be given in MeV
        return self._parameters['BeamEnergy']
    @property
    def B2(self):
        #return self._parameters['B2']
        return 1 - (0.5 * 1/self.Energy)**2
    @property
    def Ep(self):
        return self._parameters['Permitivity']
    @property
    def Mu(self):
        return self._parameters['Permeability']
    @property
    def y0(self):
        return self._parameters['Offset']*1e2
    @property
    def dz_step(self):
        return self._parameters['TotalLength'] * 1/step
        
    def calculateWakePotential(self):
        #Need the structure parameters and the bunch length range
        #No need for convergence as 1D calculation - just use a huge number of modes
        LongitudinalPositions = self._calculation['LongitudinalPositions']
        print('Calculating Eigenvalues')
        CalculationComplete = False
        CalculateFields = False
        while CalculationComplete == False:
            #Calculate the eigenvalues
            zeroesES = []
            zeroesEA= []
            zeroesHS = []
            zeroesHA = []
            zeroesES2 = []
            zeroesEA2= []
            zeroesHS2 = []
            zeroesHA2 = []
            Wx = []
            Wy = []
            WxDummy = []
            WyDummy = []
            Wz = []
            zeroesES, zeroesEA, zeroesHS, zeroesHA = DiWaCAT.EigenvaluesCalculator(50, 50, zeroesES, zeroesEA, zeroesHS, zeroesHA, self.Ep, self.Mu, self.B2, self.b, self.c, self.w, 1e-6)
            if np.isnan(zeroesES).any() or np.isnan(zeroesEA).any() or np.isnan(zeroesHS).any() or np.isnan(zeroesHA).any():
                print('Recalculating Eigenvalues')
                CalculateFields = False
            else:
                CalculateFields = True
            if CalculateFields == True:
                #Calculate the wake potential
                #Need to also define empty arrays for Fx and Ez even though we won't use them
                print('Calculating Fields')
                Wx = np.zeros(len(LongitudinalPositions)); Wy = np.zeros(len(LongitudinalPositions)); Wz = np.zeros(len(LongitudinalPositions));
                WxDummy = np.zeros(len(LongitudinalPositions)); WyDummy = np.zeros(len(LongitudinalPositions));
                for i in range(len(LongitudinalPositions)):
                    WxDummy[i], Wy[i], Wz[i] = DiWaCAT.CalcWakeElement(Wz[i], WxDummy[i], Wy[i], 0.5*self.w, 0.5*self.w, self.y0, self.y0, LongitudinalPositions[i] * 1e2, 0, zeroesES, zeroesEA, zeroesHS, zeroesHA, 50, 50, self.Ep, self.Mu, self.B2, self.b, self.c, self.w);
                    Wy[i] = Wy[i] * 1e2 *  1e6# * 1e8 * 1/const.c
                    Wx[i], WyDummy[i], Wz[i] = DiWaCAT.CalcWakeElement(Wz[i], Wx[i], WyDummy[i], 0.5*self.w + 5e-4, 0.5*self.w , 0, 0, LongitudinalPositions[i] * 1e2, 0, zeroesES, zeroesEA, zeroesHS, zeroesHA, 50, 50, self.Ep, self.Mu, self.B2, self.b, self.c, self.w);
                    Wx[i] = Wx[i] * 1e2 * 1e6 * 1/5e-6 #Add the normalisation to distance (we measure 5um away from the axis)
                    Wz[i] = Wz[i] * 1e2 *  1e6
            if np.isnan(Wy).any():
                print('Recalculating Eigenvalues')
                CalculationComplete = False
                CalculateFields == False
            else:
                CalculationComplete = True
        print('Fields Calculated')
        #These calculations are in MV/C/m so need to multiply by 1e6
        self._calculation['Dipole'] = Wy *1e6
        self._calculation['Quadrupole'] = Wx *1e6
        self._calculation['LongWake'] = Wz *1e6
        self.wake_calculated = True
    
    def calculateWakefields(self, longprofile):
        '''
        Calculate the wake potential, convolution with the wake potential
        NOTE: may need to flip the z positions - assumes negative z is the head. Will keep as is atm so consistent with DiWACAT rather than ocelot
        Return: longitudinal, dipole, and quadrupole component. All as a function of longitudinal position
        '''
        if self.wake_calculated == False:
            self.calculateWakePotential()
        LongPosition = self._calculation['LongitudinalPositions']
        DipoleWake = self._calculation['Dipole']
        QuadrupoleWake = self._calculation['Quadrupole']
        LongWake = self._calculation['LongWake']
        
        LongVariation = LongPosition[1] - LongPosition[0]
        totalcharge = np.sum(longprofile)
        DipoleArray = totalcharge*LongVariation * (1/np.trapz(longprofile,LongPosition[:-1])) * np.convolve(longprofile, DipoleWake, mode = 'full');
        quadWakeArray = totalcharge*LongVariation * (1/np.trapz(longprofile,LongPosition[:-1])) * np.convolve(longprofile, QuadrupoleWake, mode = 'full');
        EzArray = totalcharge*LongVariation * (1/np.trapz(longprofile,LongPosition[:-1])) * np.convolve(longprofile, LongWake, mode = 'full');
        return EzArray, DipoleArray, quadWakeArray
        
        
    def apply(self, p_array, dz):
        '''
        Get the current profile, calculate the wakefields, track
        '''
        ps = p_array.rparticles
        zpositions = ps[4] - np.min(ps[4])
        longprofile, _ = np.histogram(zpositions, bins = self._calculation['LongitudinalPositions'], weights = p_array.q_array)
        wakefield = self.calculateWakefields(longprofile)
        
        npointsinterpolate = len(self._calculation['LongitudinalPositions'])
        
        
        longitudinalinterp = interp1d(self._calculation['LongitudinalPositions'], wakefield[0][:npointsinterpolate])
        dipoleinterp = interp1d(self._calculation['LongitudinalPositions'], wakefield[1][:npointsinterpolate])
        quadinterp = interp1d(self._calculation['LongitudinalPositions'], wakefield[2][:npointsinterpolate])
        
        '''
        For the quadrupole wakefield we need a moving average position
        If the average position moves, the zero point for the quadrupole component moves too
        '''
        indices_zsorted = np.argsort(zpositions)
        x_sorted = p_array.rparticles[0][indices_zsorted]
        y_sorted = p_array.rparticles[2][indices_zsorted]
        z_sorted = p_array.rparticles[4][indices_zsorted]
        # Calculate cumulative sum of y_sorted and x_sorted
        cumulative_sum_y = np.cumsum(y_sorted)
        cumulative_sum_x = np.cumsum(x_sorted)
        # Calculate moving average
        moving_avg_y = cumulative_sum_y / np.arange(1, len(y_sorted) + 1)
        moving_avg_y_interpolate = interp1d(z_sorted, moving_avg_y, fill_value = "extrapolate")
        moving_avg_x = cumulative_sum_x / np.arange(1, len(x_sorted) + 1)
        moving_avg_x_interpolate = interp1d(z_sorted, moving_avg_x, fill_value = "extrapolate")

        p_array.rparticles[5] = p_array.rparticles[5] + longitudinalinterp(zpositions) * dz * 1/(p_array.E * 1e9)
        if self._parameters['Horizontal']:
            p_array.rparticles[3] = p_array.rparticles[3] + ((dipoleinterp(zpositions) + -1*quadinterp(zpositions)* (p_array.rparticles[2] - moving_avg_y_interpolate(zpositions))) * dz  * 1/(p_array.E * 1e9))
            p_array.rparticles[1] = p_array.rparticles[1] + (quadinterp(zpositions) * (p_array.rparticles[0] - moving_avg_x_interpolate(zpositions)) * dz  / (p_array.E * 1e9))
        else:
            p_array.rparticles[1] = p_array.rparticles[1] + ((dipoleinterp(zpositions) + -1*quadinterp(zpositions) * (p_array.rparticles[0] - moving_avg_x_interpolate(zpositions))) * dz * 1/(p_array.E * 1e9))
            p_array.rparticles[3] = p_array.rparticles[3] + (quadinterp(zpositions) * (p_array.rparticles[2] - moving_avg_y_interpolate(zpositions)) * dz  / (p_array.E * 1e9))
        