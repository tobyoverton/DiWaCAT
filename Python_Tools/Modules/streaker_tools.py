# -*- coding: utf-8 -*-
"""
Streaker Reconstruction - Header File
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import scipy
import scipy.interpolate, scipy.optimize
import time
from Python_Tools.Modules import beam_tools as dbt
from fastkde import fastKDE

from Python_Tools.Modules import diwacat_tools as DWA

from scipy import signal

cgstoSI = 1e8
c = 299792458

def Array2Dto1D(array):
    #Converts an array in the form [[x1,y1], [x2,y2],...] into [[x1,x2,...], [y1,y2, ...]] as needed
    ArrayXValues = [];
    ArrayYValues = [];
    for i in range(len(array)):
        ArrayXValues.append(array[i][0]);
        ArrayYValues.append(array[i][1]);
    return np.array([ArrayXValues,ArrayYValues])

def MeanVariance(array):
    #Returns the mean and variance from the arrays used here
    mean = np.average(array[0], weights = array[1])
    variance = np.average((array[0]-mean)**2, weights = array[1])
    return mean, variance

def changeArrayNPoints(array, nPoints, maintainIntegral=True):
    #Take the current array, make an interpolating function and create new array between two values
    #Needed to increase the sample size of a function (e.g. the longitudinal profile to make a streak)
    InitialIntegral = np.trapz(array[1],array[0])
    CurrentArrayInterpolate = scipy.interpolate.interp1d(array[0], array[1], bounds_error = False, fill_value = 0, kind = 'linear');
    #nPoints = round((max(array[0])-min(array[0]))/newDistance);
    newArrayValues = np.linspace(min(array[0]), max(array[0]), num=nPoints);
    newValues = CurrentArrayInterpolate(newArrayValues);
    NormalisingFactor = 1
    if maintainIntegral==1:
        NormalisingFactor = InitialIntegral* 1/np.trapz(newValues, newArrayValues)
    return [newArrayValues,NormalisingFactor*newValues];

def changeArraySeparation(array, newSeparation, maintainIntegral = False):
    #Change the separation between array x values
    #Convolution needs the x bins to have the same width -> need this to resample one of the arrays
    interpolatearray = scipy.interpolate.interp1d(array[0], array[1], bounds_error=False, fill_value = None)
    npoints = round((array[0][-1] - array[0][0])*1/newSeparation)+1
    xvalues = np.linspace(array[0][0], array[0][-1], num = npoints)
    yvalues = interpolatearray(xvalues)
    if maintainIntegral == True:
        yvalues = yvalues * 1/np.trapz(yvalues,xvalues)
    return [xvalues,yvalues]

def GaussianBunch(BunchLength, NumberOfPoints):
    #Gives a longitudinally Gaussian profile
    rmslength = BunchLength*c
    z = np.linspace(0, 6 * rmslength)
    charge = np.exp(-0.5*((z - 3*rmslength)/rmslength)**2)
    charge = charge * 1/np.trapz(charge, x = z)
    return [z, charge]


def GreenProfileConvolution(longprofile, green_function, Profile_Only = False):
    if (longprofile[0][1] - longprofile[0][0]) != (green_function[0][1]-green_function[0][0]):
        longprofile=changeArraySeparation(longprofile, green_function[0][1]-green_function[0][0], maintainIntegral=True) 
    FyArray = (green_function[0][1]-green_function[0][0]) * (1/np.trapz(longprofile[1],longprofile[0])) * np.convolve(longprofile[1], green_function[1], mode = 'full');
    FyArray = [np.linspace(0,len(FyArray)*(green_function[0][1]-green_function[0][0]), len(FyArray)), FyArray]
    if Profile_Only == True:
        EndOfProfileIndex = np.where(FyArray[0] < max(longprofile[0]))[0]
        FyArray = [FyArray[0][EndOfProfileIndex],FyArray[1][EndOfProfileIndex]]
    return FyArray

def BackwardPropagate(StreakProfile, fy_Array, NBins, NSample):
    #Take the Fy profile and backwards propagate the screen positions
    
    #First we want to make sure the profile is monotonic
    fy_monotonic = [fy_Array[0][0:np.argmax(fy_Array[1])], fy_Array[1][0:np.argmax(fy_Array[1])]]
    profileToPropagate = changeArrayNPoints(StreakProfile, NSample)
    CoordinateTransformFunction = scipy.interpolate.interp1d(fy_monotonic[1],fy_monotonic[0],fill_value=-1e3, bounds_error=False)
    times = CoordinateTransformFunction(profileToPropagate[0])
    t_histogram = np.histogram(times[np.nonzero(times+1e3)], weights=profileToPropagate[1][np.nonzero(times+1e3)], bins = NBins)
    t_histogram = [t_histogram[1][1:], scipy.signal.savgol_filter(t_histogram[0], 3, 1)]
    t_histogram[1] = t_histogram[1] * 1/np.trapz(t_histogram[1],t_histogram[0])
    t_histogram[0] = t_histogram[0] - t_histogram[0][0]
    return t_histogram

def ForwardPropagate(long_profile, fy_array, NSample):
    #Propagate the longitudinal profile using an Fy profile
    #Can always use GreenProfileConvolution as the second argument if need to calculate Fy
    
    FyFunction = scipy.interpolate.interp1d(fy_array[0], fy_array[1], bounds_error=False, fill_value = -2)
    
    #Increase the number of longitudinal points
    longprofile = changeArrayNPoints(long_profile, NSample)
    
    yPositions = FyFunction(longprofile[0])
    yPositions = [yPositions[np.nonzero(yPositions+2)], longprofile[1][np.nonzero(yPositions+2)]]
    return yPositions

def BackwardReconstruction(StreakProfile, GreensFunction, BunchLengths, nBins, numberofiterations, returnCompleteData = False):
    #Begin with a Gaussian profile with RMS length of initial
    
    measuredprofiles = []
    variations = []
    for j in range(len(BunchLengths)):
        profile_meas = GaussianBunch(BunchLengths[j],NBins)
        
        for i in range(numberofiterations):
            #First off calculate the Fy profile by convolving with the Green's function
            #We only want the Green's function in the range of the beam itself
            FyArray = GreenProfileConvolution(profile_meas, GreensFunction)
                    
            #Backwards propagate the measured streak profile using the Fy function calculated
            #This gives the new measured profile
            profile_meas = BackwardPropagate(StreakProfile, FyArray, NBins, 5000)
            #Add zero charge at the start or end of the distribution if needed - avoids histogram errors in some cases
            if profile_meas[1][0] != 0:
                profile_meas[1] = np.insert(profile_meas[1],0,0)[0:-1]
                profile_meas[0] = np.insert(profile_meas[0],0,-profile_meas[0][1])[0:-1]
                profile_meas[0] = profile_meas[0] - profile_meas[0][0]
            if profile_meas[1][-1] != 0:
                profile_meas[1] = np.append(profile_meas[1],0)
                profile_meas[0] = np.append(profile_meas[0],profile_meas[0][-1]+profile_meas[0][1])
                profile_meas[0] = profile_meas[0] - profile_meas[0][0]    
            
    
            MeasuredBeamStreak = ForwardPropagate(profile_meas, GreenProfileConvolution(profile_meas, GreensFunction), 5000)
            mv_sim = MeanVariance(MeasuredBeamStreak)
            mv_meas = MeanVariance(StreakProfile)
            
            CostFunction = (1 + abs(mv_sim[0] - mv_meas[0]))*(1 + abs(mv_sim[1] - mv_meas[1]))
            
            measuredprofiles.append(profile_meas)
            variations.append(CostFunction)
    
        #Provide the option to give the output of each backward propagation
    if returnCompleteData==True:
        return measuredprofiles, variations
    #If not return the profile that gave the minimum cost function
    else:
        return measuredprofiles[np.argmin(variations)], variations[np.argmin(variations)]

