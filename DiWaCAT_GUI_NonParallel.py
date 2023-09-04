# -*- coding: utf-8 -*-
"""
DiWaCAT GUI using PyQT
"""

# Only needed for access to command line arguments
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget
from PyQt5.QtCore import Qt, QSize, QThread, pyqtSignal, QObject
import PyQt5.QtWidgets as pyqt
import pyqtgraph as pg
import sys, os, subprocess
from pathlib import Path
#from time import sleep
import numpy as np
import copy
import scipy.constants as const
from pathos.helpers import mp
import multiprocessing

import threading



#import diwacat_tools as DWA
#import beam_tools as dbt
from Python_Tools.Modules import diwacat_tools as DWA
from Python_Tools.Modules import beam_tools as dbt

import matplotlib.pyplot as plt
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']


def UnitConversion(prefix):
    if prefix == 'unit' or prefix == 'Unit':
        return 1
    else:
        return 1/getattr(const, prefix)
'''
class OpenFile(QWidget):
    def _init_(self):
        super()._init_()
        self.setWindowTitle("Open File")
        self.FileDialog = pyqt.QFileDialog.getOpenFileName(None, 'Test Dialog', os.getcwd(), 'HDF5 Files (*.h5*)')
        self.FileDialog.clicked.connect(self.buttonClicked)
        self.filename = ''
    def buttonClicked(self):
        options = pyqt.QFileDialog.Options()
        options |= pyqt.QFileDialog.DontUseNativeDialog
        self.filename, _ = pyqt.QFileDialog.getOpenFileName(self,"Choose File", "","csv (*.csv)", options=options)
    def GetFileName(self):
        return self.filename
'''
class BeamMakeWindow(QWidget):
    """
    Make a beam and save the file
    """
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Beam File Maker")
        windowWidth = 1000
        windowHeight = 800
        self.setMinimumSize(windowWidth, windowHeight)
        
        self.TotalLayout = pyqt.QHBoxLayout()
        self.BeamParameters = pyqt.QGridLayout()
        self.BeamGraphs = pyqt.QGridLayout()
        
        self.scalelist = ['unit','kilo', 'mega', 'giga','milli', 'micro', 'nano', 'pico', 'femto']
        
        self.nMacros = pyqt.QSpinBox()
        self.nMacros.setRange(1000,2147483647)
        self.nMacros.setValue(32000)
        self.nMacros.setSingleStep(1000)
        self.BeamParameters.addWidget(pyqt.QLabel("Number of Macroparticles:"),0,0,1,2)
        self.BeamParameters.addWidget(self.nMacros,0,1)

        self.MomentumValue = pyqt.QDoubleSpinBox()
        self.MomentumValue.setMinimum(0.1)
        self.MomentumValue.setMaximum(10000)
        self.MomentumScale = pyqt.QComboBox()
        self.MomentumScale.addItems(self.scalelist)
        self.BeamParameters.addWidget(pyqt.QLabel('Beam Momentum [eV/c]:'),1,0)
        self.BeamParameters.addWidget(self.MomentumValue,1,1)
        self.BeamParameters.addWidget(self.MomentumScale,1,2)
        
        self.MomentumSigma = pyqt.QDoubleSpinBox()
        self.MomentumSigma.setMinimum(0.1)
        self.MomentumSigma.setMaximum(10000)
        self.MomentumSigmaScale = pyqt.QComboBox()
        self.MomentumSigmaScale.addItems(self.scalelist)
        self.BeamParameters.addWidget(pyqt.QLabel('Uncorrelated RMS Momentum Variation [eV/c]:'),2,0)
        self.BeamParameters.addWidget(self.MomentumSigma,2,1)
        self.BeamParameters.addWidget(self.MomentumSigmaScale,2,2)
        
        self.ChirpGradient = pyqt.QDoubleSpinBox()
        self.ChirpGradient.setRange(-1000,1000)
        self.ChirpGradient.setSingleStep(1)
        self.ChirpGradientScale = pyqt.QComboBox()
        self.ChirpGradientScale.addItems(self.scalelist)
        self.BeamParameters.addWidget(pyqt.QLabel('Correlated Momentum Gradient [eV/c/ps]:'),3,0)
        self.BeamParameters.addWidget(self.ChirpGradient,3,1)
        self.BeamParameters.addWidget(self.ChirpGradientScale,3,2)
        
        self.ChargeValue = pyqt.QDoubleSpinBox()
        self.ChargeValue.setMinimum(0.1)
        self.ChargeValue.setMaximum(10000)
        self.ChargeScale = pyqt.QComboBox()
        self.ChargeScale.addItems(self.scalelist)
        self.BeamParameters.addWidget(pyqt.QLabel('Total Charge [C]:'),4,0)
        self.BeamParameters.addWidget(self.ChargeValue,4,1)
        self.BeamParameters.addWidget(self.ChargeScale,4,2)
        
        
        self.epsx = pyqt.QDoubleSpinBox()
        self.epsx.setMinimum(0.1)
        self.epsx.setMaximum(10000)
        self.epsxScale = pyqt.QComboBox()
        self.epsxScale.addItems(self.scalelist)
        self.BeamParameters.addWidget(pyqt.QLabel('Normalised Emittance x [m rad]:'),5,0)
        self.BeamParameters.addWidget(self.epsx,5,1)
        self.BeamParameters.addWidget(self.epsxScale,5,2)
        
        self.epsy = pyqt.QDoubleSpinBox()
        self.epsy.setMinimum(0.1)
        self.epsy.setMaximum(10000)
        self.epsyScale = pyqt.QComboBox()
        self.epsyScale.addItems(self.scalelist)
        self.BeamParameters.addWidget(pyqt.QLabel('Normalised Emittance y [m rad]:'),6,0)
        self.BeamParameters.addWidget(self.epsy,6,1)
        self.BeamParameters.addWidget(self.epsyScale,6,2)
        
        self.sigx = pyqt.QDoubleSpinBox()
        self.sigx.setMinimum(0.1)
        self.sigx.setMaximum(10000)
        self.sigxScale = pyqt.QComboBox()
        self.sigxScale.addItems(self.scalelist)
        self.BeamParameters.addWidget(pyqt.QLabel('Horizontal RMS Waist Size [m]:'),7,0)
        self.BeamParameters.addWidget(self.sigx,7,1)
        self.BeamParameters.addWidget(self.sigxScale,7,2)
        
        self.sigy = pyqt.QDoubleSpinBox()
        self.sigy.setMinimum(0.1)
        self.sigy.setMaximum(10000)
        self.sigyScale = pyqt.QComboBox()
        self.sigyScale.addItems(self.scalelist)
        self.BeamParameters.addWidget(pyqt.QLabel('Vertical RMS Waist Size [m]:'),8,0)
        self.BeamParameters.addWidget(self.sigy,8,1)
        self.BeamParameters.addWidget(self.sigyScale,8,2)
        
        self.sigz = pyqt.QDoubleSpinBox()
        self.sigz.setMinimum(0.1)
        self.sigz.setMaximum(10000)
        self.sigzScale = pyqt.QComboBox()
        self.sigzScale.addItems(self.scalelist)
        self.BeamParameters.addWidget(pyqt.QLabel('Longitudinal RMS Bunch Length [m or s]:'),9,0)
        self.BeamParameters.addWidget(self.sigz,9,1)
        self.BeamParameters.addWidget(self.sigzScale,9,2,1,1)
        
        self.TorZ = pyqt.QComboBox()
        self.TorZ.addItems(['Distance', 'Time'])
        self.BeamParameters.addWidget(pyqt.QLabel('Bunch Length Scale:'),10,0)
        self.BeamParameters.addWidget(self.TorZ,10,1)
        
        self.ProfileShape = pyqt.QComboBox()
        self.ProfileShape.addItems(['Gaussian', 'SkewGaussian', 'Uniform', 'Plateau', 'DoubleGauss'])
        self.ProfileShape.currentIndexChanged.connect( self.ProfileShapeChange )
        self.BeamParameters.addWidget(pyqt.QLabel('Longitudinal Profile Shape:'),11,0)
        self.BeamParameters.addWidget(self.ProfileShape,11,1)
        self.LastProfile = self.ProfileShape.currentText()
        
        self.HorizontalWaistPosition = pyqt.QDoubleSpinBox()
        self.HorizontalWaistPosition.setMinimum(-5000)
        self.HorizontalWaistPosition.setSingleStep(0.1)
        self.HorizontalWaistPosition.setValue(0)
        self.BeamParameters.addWidget(pyqt.QLabel('Horizontal Waist Position [m]'), 15, 0)
        self.BeamParameters.addWidget(self.HorizontalWaistPosition,15,1)
        
        self.VerticalWaistPosition = pyqt.QDoubleSpinBox()
        self.VerticalWaistPosition.setMinimum(-5000)
        self.VerticalWaistPosition.setSingleStep(0.1)
        self.VerticalWaistPosition.setValue(0)
        self.BeamParameters.addWidget(pyqt.QLabel('Vertical Waist Position [m]'), 16, 0)
        self.BeamParameters.addWidget(self.VerticalWaistPosition,16,1)
        
        self.BeamParameters.addWidget(pyqt.QLabel('Mesh Variables:'),17,0)
        self.DefaultMesh = pyqt.QPushButton('Default Mesh Parameters')
        self.DefaultMesh.clicked.connect(self.DefaultMeshPressed)
        self.BeamParameters.addWidget(self.DefaultMesh, 17,1,1,2)
        
        self.CellsPerSigmaT = pyqt.QDoubleSpinBox()
        self.CellsPerSigmaT.setMinimum(1.5)
        self.CellsPerSigmaT.setSingleStep(0.5)
        self.CellsPerSigmaT.setValue(3)
        self.BeamParameters.addWidget(pyqt.QLabel('Cells per Sigma (Transverse):'),18,0)
        self.BeamParameters.addWidget(self.CellsPerSigmaT,18,1)
         
        self.CellsPerSigmaL = pyqt.QDoubleSpinBox()
        self.CellsPerSigmaL.setMinimum(1.5)
        self.CellsPerSigmaL.setSingleStep(0.5)
        self.CellsPerSigmaL.setValue(3.5)
        self.BeamParameters.addWidget(pyqt.QLabel('Cells per Sigma (Longitudinal):'),19,0)
        self.BeamParameters.addWidget(self.CellsPerSigmaL,19,1)
        
        self.MaxXCell = pyqt.QDoubleSpinBox()
        self.MaxXCell.setMinimum(1)
        self.MaxXCell.setSingleStep(0.25)
        self.MaxXCell.setValue(2.25)
        self.BeamParameters.addWidget(pyqt.QLabel('Maximum X Cell Position / sigmax:'),20,0)
        self.BeamParameters.addWidget(self.MaxXCell,20,1)
        
        self.MaxYCell = pyqt.QDoubleSpinBox()
        self.MaxYCell.setMinimum(1)
        self.MaxYCell.setSingleStep(0.25)
        self.MaxYCell.setValue(2.25)
        self.BeamParameters.addWidget(pyqt.QLabel('Maximum Y Cell Position / sigmay:'),21,0)
        self.BeamParameters.addWidget(self.MaxYCell,21,1)
        
        self.MaxZCell = pyqt.QDoubleSpinBox()
        self.MaxZCell.setMinimum(1)
        self.MaxZCell.setSingleStep(0.25)
        self.MaxZCell.setValue(3)
        self.BeamParameters.addWidget(pyqt.QLabel('Maximum Z Cell Position / sigmaz:'),22,0)
        self.BeamParameters.addWidget(self.MaxZCell,22,1)
               
        self.BeamMakeButton = pyqt.QPushButton("Plot Beam Profiles")
        self.BeamMakeButton.clicked.connect(self.BeamPlot)
        self.BeamParameters.addWidget(self.BeamMakeButton,23,0,1,3)
        
        self.BeamSaveButton = pyqt.QPushButton("Save Beam")
        self.BeamSaveButton.clicked.connect(self.SaveBeam)
        self.BeamParameters.addWidget(self.BeamSaveButton,24,0,1,3)
        
        self.BeamXData = []
        self.BeamYData = []
        self.BeamZData = []
        self.BeamXPData = []
        self.BeamYPData = []
        self.BeamcpzData = []
        
        self.XYPlot = pg.PlotWidget()
        self.XYPlot.setBackground('w')
        self.XYPlot.setLabel('bottom','x [m]')
        self.XYPlot.setLabel('left','y [m]')
        self.XYPlot.setMinimumSize(int(windowWidth*0.3), int(windowHeight*0.3))
        self.xyplot = pg.ScatterPlotItem()
        self.xyplot.setData(self.BeamXData,self.BeamYData)
        
        #self.XYPlot.plot(self.BeamXData, self.BeamYData, pen = 'k')
        self.XYPlot.addItem(self.xyplot)
        
        #Q vs Z histogram
        self.QPlot = pg.PlotWidget()
        self.QPlot.setBackground('w')
        self.QPlot.setLabel('bottom', 'z [m]')
        self.QPlot.setLabel('left','Intensity [arb.]')
        self.QPlot.setMinimumSize(int(windowWidth*0.3), int(windowHeight*0.3))
        self.ZHist = [[],[0]]
        self.QCurve = pg.PlotCurveItem(self.ZHist[1], self.ZHist[0], stepMode=True, fillLevel=0)
        self.QPlot.addItem(self.QCurve)
        
        
        self.XXPPlot = pg.PlotWidget()
        self.XXPPlot.setBackground('w')
        self.XXPPlot.setLabel('bottom','x [m]')
        self.XXPPlot.setLabel('left','x\' [rads]')
        self.XXPPlot.setMinimumSize(int(windowWidth*0.3), int(windowHeight*0.3))
        self.xxpplot =self.XXPPlot.plot(self.BeamXData, self.BeamXPData, pen = 'k')
        self.XXPPlot.addItem(self.xxpplot)
        
        self.YYPPlot = pg.PlotWidget()
        self.YYPPlot.setBackground('w')
        self.YYPPlot.setLabel('bottom','y [m]')
        self.YYPPlot.setLabel('left','y\' [rads]')
        self.YYPPlot.setMinimumSize(int(windowWidth*0.3), int(windowHeight*0.3))
        self.yypplot = self.YYPPlot.plot(self.BeamYData, self.BeamYPData, pen = 'k')
        self.YYPPlot.addItem(self.yypplot)
        
        self.ZcPZPlot = pg.PlotWidget()
        self.ZcPZPlot.setBackground('w')
        self.ZcPZPlot.setLabel('bottom','z [m]')
        self.ZcPZPlot.setLabel('left','cpz [eV/c]')
        self.ZcPZPlot.setMinimumSize(int(windowWidth*0.3), int(windowHeight*0.3))
        self.zcpzplot = self.ZcPZPlot.plot(self.BeamZData, self.BeamcpzData, pen = 'k')
        self.ZcPZPlot.addItem(self.zcpzplot)
        
        self.cpzPlot = pg.PlotWidget()
        self.cpzPlot.setBackground('w')
        self.cpzPlot.setLabel('bottom', 'cpz [eV/c]')
        self.cpzPlot.setLabel('left','Intensity [arb.]')
        self.cpzPlot.setMinimumSize(int(windowWidth*0.3), int(windowHeight*0.3))
        self.cpzHist = [[],[0]]
        self.cpzCurve = pg.PlotCurveItem(self.cpzHist[1], self.cpzHist[0], stepMode=True, fillLevel=0)
        self.cpzPlot.addItem(self.cpzCurve)
        
        self.BeamGraphs.addWidget(self.XYPlot, 0, 0)
        self.BeamGraphs.addWidget(self.QPlot, 0, 1)
        self.BeamGraphs.addWidget(self.XXPPlot, 1, 0)
        self.BeamGraphs.addWidget(self.YYPPlot, 1, 1)
        self.BeamGraphs.addWidget(self.ZcPZPlot, 2, 0)
        self.BeamGraphs.addWidget(self.cpzPlot, 2, 1)
        
               
        self.TotalLayout.addLayout(self.BeamParameters)
        self.TotalLayout.addLayout(self.BeamGraphs)
        self.setLayout(self.TotalLayout)
        
    def ProfileShapeChange(self):
        #BeamOptions for Specific Profile Shapes
        if self.LastProfile == 'SkewGaussian':
            self.SkewLabel.setVisible(False)
            self.Skew.setVisible(False)
        if self.LastProfile == 'Plateau':
            self.PlateauLabel.setVisible(False)
            self.PlateauTime.setVisible(False)
        if self.LastProfile == 'DoubleGauss':
            self.SecondGausSigmaLabel.setVisible(False)
            self.SecondGausSigma.setVisible(False)
            self.SecondGausSigmaScale.setVisible(False)
            self.SecondGausAmpLabel.setVisible(False)
            self.SecondGausAmp.setVisible(False)
            self.SecondGausOffsetLabel.setVisible(False)
            self.SecondGausOffset.setVisible(False)
            self.SecondGausOffsetScale.setVisible(False)
        
        self.Skew = pyqt.QDoubleSpinBox()
        self.SkewLabel = pyqt.QLabel('Skew Factor (alpha):')
        self.Skew.setSingleStep(0.5)
        if self.ProfileShape.currentText() == 'SkewGaussian':
            self.BeamParameters.addWidget(self.SkewLabel,12,0)
            self.BeamParameters.addWidget(self.Skew,12,1)
        
        self.PlateauTime = pyqt.QDoubleSpinBox()
        self.PlateauLabel = pyqt.QLabel('Plateau Rise Time / sigmaz')
        self.PlateauTime.setSingleStep(0.5)
        self.PlateauTime.setRange(0,3)
        if self.ProfileShape.currentText() == 'Plateau':
            self.BeamParameters.addWidget(self.PlateauLabel,12,0)
            self.BeamParameters.addWidget(self.PlateauTime,12,1)
            
        self.SecondGausSigmaLabel = pyqt.QLabel('Second Gaussian RMS Length [m or s]')
        self.SecondGausSigma = pyqt.QDoubleSpinBox()
        self.SecondGausSigma.setRange(0,999)
        self.SecondGausSigmaScale = pyqt.QComboBox()
        self.SecondGausSigmaScale.addItems(self.scalelist)
        self.SecondGausAmpLabel = pyqt.QLabel('Second Gaussian Relative Amplitude')
        self.SecondGausAmp = pyqt.QDoubleSpinBox()
        self.SecondGausAmp.setRange(0,1e6)
        self.SecondGausOffsetLabel = pyqt.QLabel('Gaussian Peak to Peak Distance [m or s]')
        self.SecondGausOffset = pyqt.QDoubleSpinBox()
        self.SecondGausOffset.setRange(0,999)
        self.SecondGausOffsetScale = pyqt.QComboBox()
        self.SecondGausOffsetScale.addItems(self.scalelist)
        if self.ProfileShape.currentText() == 'DoubleGauss':
            self.BeamParameters.addWidget(self.SecondGausSigmaLabel,12,0)
            self.BeamParameters.addWidget(self.SecondGausSigma,12,1)
            self.BeamParameters.addWidget(self.SecondGausSigmaScale,12,2)
            self.BeamParameters.addWidget(self.SecondGausAmpLabel,13,0)
            self.BeamParameters.addWidget(self.SecondGausAmp,13,1)
            self.BeamParameters.addWidget(self.SecondGausOffsetLabel,14,0)
            self.BeamParameters.addWidget(self.SecondGausOffset,14,1)
            self.BeamParameters.addWidget(self.SecondGausOffsetScale,14,2)
                       
        self.LastProfile = self.ProfileShape.currentText()
        
    def DefaultMeshPressed(self):
        self.CellsPerSigmaT.setValue(3)
        self.CellsPerSigmaL.setValue(3.5)
        self.MaxXCell.setValue(2.25)
        self.MaxYCell.setValue(2.25)
        self.MaxZCell.setValue(3)
        
    def MakeBeam(self):
        N_MACROS = self.nMacros.value()
        #Make empty dictionary for the beam parameters
        bp = {}
        bp["pz_MeV"] = self.MomentumValue.value() * 1/UnitConversion(self.MomentumScale.currentText()) * 1e-6
        bp["eps_x_N"] = self.epsx.value() * 1/UnitConversion(self.epsxScale.currentText())
        bp["eps_y_N"] = self.epsy.value() * 1/UnitConversion(self.epsyScale.currentText())
        bp["sig_x_0"] = self.sigx.value() * 1/UnitConversion(self.sigxScale.currentText())
        bp["sig_y_0"] = self.sigy.value() * 1/UnitConversion(self.sigyScale.currentText())
        SIGMA_LONG = self.sigz.value() * 1/UnitConversion(self.sigzScale.currentText())
        if self.TorZ.currentText() == 'Time':
            SIGMA_LONG *= const.c
        bp["sig_z_0"] = SIGMA_LONG
        bp["sig_pz_0"] = self.MomentumSigma.value() * 1/UnitConversion(self.MomentumSigmaScale.currentText()) * 1e-6
        bp["lCorr_Fac"] = 0#1  # An arbitrary number, but want to convert to chirp at some point
        bp["x_0"] = 0
        bp["y_0"] = 0
        bp["xp_0"] = 0
        bp["yp_0"] = 0
        bp["z_0"] = 0
        TOTAL_CHARGE = self.ChargeValue.value() * 1/UnitConversion(self.ChargeScale.currentText()) #Bunch charge in C
        bp["charge_per_macro"] = TOTAL_CHARGE / N_MACROS #Charge per macro
        
        bp["LongitudinalProfile"] = self.ProfileShape.currentText()
        if self.ProfileShape.currentText() == 'SkewGaussian':
            bp["skew"]= self.Skew.value();
        elif self.ProfileShape.currentText() == 'Plateau':
            bp["plat_rise"] =  self.PlateauTime.value() * SIGMA_LONG #For plateau type beam (rise and fall time)
        elif self.ProfileShape.currentText() == 'DoubleGauss':
            bp["sig_z_2"] = self.SecondGausSigma.value() * 1/UnitConversion(self.SecondGausSigmaScale.currentText())
            bp["offset"] = self.SecondGausOffset.value() * 1/UnitConversion(self.SecondGausOffsetScale.currentText())
            if (self.TorZ.currentText() == 'Time'):
                bp["sig_z_2"] = bp["sig_z_2"] * const.c
                bp["offset"] = bp["offset"] * const.c
            bp["rel_amp"] = self.SecondGausAmp.value()

        myBeam = dbt.GaussianBeam()
        myBeam.setBeamParameters(bp)
        myBeam.calcWaistOpticFunctions()
        myBeam.generate6DMacros(N_MACROS)
        
        myBeam._beam['x'] += self.HorizontalWaistPosition.value() * myBeam.xp
        myBeam._beam['y'] += self.VerticalWaistPosition.value() * myBeam.yp
        t = myBeam.z * 1/const.c * 1e12
        Chirp = self.ChirpGradient.value() * 1/UnitConversion(self.ChirpGradientScale.currentText()) * 1/const.c * myBeam.z * 1e12
        print(max(Chirp))
        myBeam._beam['pz'] += Chirp * myBeam.q_over_c
        #myBeam._beam['cpz'] += Chirp
        
        #Mesh the beam
        Lx = 2* self.MaxXCell.value() * myBeam.Sx
        Dx = Lx/(2* self.MaxXCell.value() *self.CellsPerSigmaT.value())
        Ly = 2* self.MaxYCell.value() * myBeam.Sy
        Dy = Ly/(2* self.MaxYCell.value() *self.CellsPerSigmaT.value())
        Lz = 2* self.MaxZCell.value() * myBeam.Sz
        Dz = Lz/(2* self.MaxZCell.value() *self.CellsPerSigmaL.value())
        
        binnedBeam = dbt.beamBinner(myBeam)
        binnedBeam.binTheBeam(Lx=Lx,Dx=Dx,Ly=Ly,Dy=Dy,Lz=Lz,Dz=Dz)
        print("Fraction of total charge captured on this mesh is",binnedBeam.nMacroHist/myBeam.nMacros)
        return binnedBeam      
        
    def BeamPlot(self):
        beam = self.MakeBeam()
        
        self.BeamXData = beam.beam.x
        self.BeamYData = beam.beam.y
        self.BeamZData = beam.beam.z
        self.BeamXPData = beam.beam.xp
        self.BeamYPData = beam.beam.yp
        self.BeamcpzData = beam.beam.cpz
        
        self.XYPlot.removeItem(self.xyplot)
        self.QPlot.removeItem(self.QCurve)
        self.XXPPlot.removeItem(self.xxpplot)
        self.YYPPlot.removeItem(self.yypplot)
        self.ZcPZPlot.removeItem(self.zcpzplot)
        self.cpzPlot.removeItem(self.cpzCurve)
        
        self.xyplot = pg.ScatterPlotItem(size = 1)
        self.xyplot.setData(self.BeamXData,self.BeamYData)
        self.xxpplot = pg.ScatterPlotItem(size = 1)
        self.xxpplot.setData(self.BeamXData,self.BeamXPData)
        self.yypplot = pg.ScatterPlotItem(size = 1)
        self.yypplot.setData(self.BeamYData,self.BeamYPData)
        self.zcpzplot = pg.ScatterPlotItem(size = 1)
        self.zcpzplot.setData(self.BeamZData,self.BeamcpzData)
        
        self.ZHist = np.histogram(self.BeamZData, bins = np.linspace(min(self.BeamZData), max(self.BeamZData), 100))
        self.cpzHist = np.histogram(self.BeamcpzData, bins = np.linspace(min(self.BeamcpzData), max(self.BeamcpzData), 100))
        #print(self.ZHist)
        self.QCurve = pg.PlotCurveItem(self.ZHist[1], self.ZHist[0], stepMode=True, fillLevel=0, pen = 'k')
        self.cpzCurve = pg.PlotCurveItem(self.cpzHist[1], self.cpzHist[0], stepMode=True, fillLevel=0, pen = 'k')
        
        self.XYPlot.addItem(self.xyplot)
        self.QPlot.addItem(self.QCurve)
        self.XXPPlot.addItem(self.xxpplot)
        self.YYPPlot.addItem(self.yypplot)
        self.ZcPZPlot.addItem(self.zcpzplot)
        self.cpzPlot.addItem(self.cpzCurve)
        
    def SaveBeam(self):
        beam = self.MakeBeam()
        filename, ok = pyqt.QFileDialog.getSaveFileName(self,"Select a File", "", "HDF5 Files (*.h5)")
        if filename:
            path = Path(filename)
        beam.write_hd5f_mesh_file(path,outputMacros=True, overwrite_existing_file=True,includeSmoothed=False,includeUnSmoothed=True)
        print('File Written')
        self.BeamPlot()
        
class BeamMeshWindow(QWidget):
    """
    Load in a beam file and mesh
    """        
    
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Mesh Existing Beam File")
        self.setMinimumWidth(400)
        
        self.TotalLayout = pyqt.QGridLayout()
        
        self.LoadButton = pyqt.QPushButton('Load Beam File')
        self.LoadButton.clicked.connect(self.OpenBeamFile)
        self.TotalLayout.addWidget(self.LoadButton, 0, 0, 1, 5)
        
        self.BeamPath = pyqt.QLineEdit()
        self.BeamPath.setText('No File Selected')
        self.TotalLayout.addWidget(pyqt.QLabel('File:'), 1,0)
        self.TotalLayout.addWidget(self.BeamPath, 1,1,1,4)
        
        self.DefaultMesh = pyqt.QPushButton('Default Mesh Parameters')
        self.DefaultMesh.clicked.connect(self.DefaultMeshPressed)
        self.TotalLayout.addWidget(self.DefaultMesh, 2,0,1,5)
        
        self.CellsPerSigmaT = pyqt.QDoubleSpinBox()
        self.CellsPerSigmaT.setMinimum(1.5)
        self.CellsPerSigmaT.setSingleStep(0.5)
        self.CellsPerSigmaT.setValue(3)
        self.TotalLayout.addWidget(pyqt.QLabel('Cells per Sigma (Transverse):'),3,0,1,2)
        self.TotalLayout.addWidget(self.CellsPerSigmaT,3,2,1,3)
         
        self.CellsPerSigmaL = pyqt.QDoubleSpinBox()
        self.CellsPerSigmaL.setMinimum(1.5)
        self.CellsPerSigmaL.setSingleStep(0.5)
        self.CellsPerSigmaL.setValue(3.5)
        self.TotalLayout.addWidget(pyqt.QLabel('Cells per Sigma (Longitudinal):'),4,0,1,2)
        self.TotalLayout.addWidget(self.CellsPerSigmaL,4,2,1,3)
        
        self.MaxXCell = pyqt.QDoubleSpinBox()
        self.MaxXCell.setMinimum(1)
        self.MaxXCell.setSingleStep(0.25)
        self.MaxXCell.setValue(2.25)
        self.TotalLayout.addWidget(pyqt.QLabel('Maximum X Cell Position / sigmax:'),5,0,1,2)
        self.TotalLayout.addWidget(self.MaxXCell,5,2,1,3)
        
        self.MaxYCell = pyqt.QDoubleSpinBox()
        self.MaxYCell.setMinimum(1)
        self.MaxYCell.setSingleStep(0.25)
        self.MaxYCell.setValue(2.25)
        self.TotalLayout.addWidget(pyqt.QLabel('Maximum Y Cell Position / sigmay:'),6,0,1,2)
        self.TotalLayout.addWidget(self.MaxYCell,6,2,1,3)
        
        self.MaxZCell = pyqt.QDoubleSpinBox()
        self.MaxZCell.setMinimum(1)
        self.MaxZCell.setSingleStep(0.25)
        self.MaxZCell.setValue(3)
        self.TotalLayout.addWidget(pyqt.QLabel('Maximum Z Cell Position / sigmaz:'),7,0,1,2)
        self.TotalLayout.addWidget(self.MaxZCell,7,2,1,3)
        
        self.MeshAndSaveButton = pyqt.QPushButton('Save Meshed Beam')
        self.MeshAndSaveButton.clicked.connect(self.MeshAndSave)
        self.TotalLayout.addWidget(self.MeshAndSaveButton,8,0,1,5)
        
        self.setLayout(self.TotalLayout)

    def OpenBeamFile(self):
        filename, ok = pyqt.QFileDialog.getOpenFileName(self,"Select a File", "", "HDF5 Files (*.h5)")
        if filename:
            self.BeamPath.setText(filename)
            
    def DefaultMeshPressed(self):
        self.CellsPerSigmaT.setValue(3)
        self.CellsPerSigmaL.setValue(3.5)
        self.MaxXCell.setValue(2.25)
        self.MaxYCell.setValue(2.25)
        self.MaxZCell.setValue(3)
        
    def MeshAndSave(self):
        BeamIn = dbt.BeamFromFile()
        try:
            BeamIn.read_WakeCode_beam_file(self.BeamPath.text())
            #Mesh the beam
            Lx = 2* self.MaxXCell.value() * BeamIn.Sx
            Dx = Lx/(2* self.MaxXCell.value() *self.CellsPerSigmaT.value())
            Ly = 2* self.MaxYCell.value() * BeamIn.Sy
            Dy = Ly/(2* self.MaxYCell.value() *self.CellsPerSigmaT.value())
            Lz = 2* self.MaxZCell.value() * BeamIn.Sz
            Dz = Lz/(2* self.MaxZCell.value() *self.CellsPerSigmaL.value())
            
            binnedBeam = dbt.beamBinner(BeamIn)
            binnedBeam.binTheBeam(Lx=Lx,Dx=Dx,Ly=Ly,Dy=Dy,Lz=Lz,Dz=Dz)
            print("Fraction of total charge captured on this mesh is",binnedBeam.nMacroHist/BeamIn.nMacros)
            
            filename, ok = pyqt.QFileDialog.getSaveFileName(self,"Select a File", "", "HDF5 Files (*.h5)")
            if filename:
                path = Path(filename)
            binnedBeam.write_hd5f_mesh_file(path,outputMacros=True, overwrite_existing_file=True,includeSmoothed=False,includeUnSmoothed=True)
            print('File Written')
        except:
            print('File Unable to Load')
        

class FieldPlotWindow(QWidget):
    """
    Definitions for the field plotter window. Includes loading a file
    """
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Field Profile Plotter")
        
        self.TotalLayout = pyqt.QHBoxLayout()
        self.OptionsLayout = pyqt.QGridLayout();
        
        self.beams = []
        self.beamlabels = []
        
        #Dialog to open a .h5 file
        self.LoadBeam = pyqt.QPushButton('Load Beam')
        self.LoadBeam.clicked.connect(self.BeamLoad)
        #self.LoadBeam.setFixedWidth(0.25*windowWidth)
        self.OptionsLayout.addWidget(self.LoadBeam,0,0,1,2)
        
        self.BeamList = pyqt.QListWidget()
        self.BeamList.addItems(self.beams)
        #self.BeamList.setFixedWidth(0.25*windowWidth)
        #self.BeamList.setMinimumHeight(0.2*windowHeight)
        self.OptionsLayout.addWidget(self.BeamList,1,0,1,2)
        
        self.BeamLabelText = pyqt.QLineEdit()
        self.BeamLabelText.setPlaceholderText('Beam Label')
        #self.BeamLabelText.setFixedWidth(0.12*windowWidth)
        self.OptionsLayout.addWidget(self.BeamLabelText,2,0)
        
        self.ChangeLabelButton = pyqt.QPushButton('Change Label')
        self.ChangeLabelButton.clicked.connect(self.ChangeLabel)
        #self.ChangeLabelButton.setFixedWidth(0.12*windowWidth)
        self.OptionsLayout.addWidget(self.ChangeLabelButton,2,1)
        
        self.RemoveBeamButton = pyqt.QPushButton('Remove Beam')
        self.RemoveBeamButton.clicked.connect(self.RemoveBeam)
        #self.RemoveBeamButton.setFixedWidth(0.25*windowWidth)
        self.OptionsLayout.addWidget(self.RemoveBeamButton,3,0,1,2)
        
        self.ForceSelect = pyqt.QComboBox();
        self.ForceSelect.addItems(["Fx", "Fy", "Ez"])
        self.ForceSelect.currentIndexChanged.connect( self.PlotData )
        self.OptionsLayout.addWidget(pyqt.QLabel("Field"),4,0)
        self.OptionsLayout.addWidget(self.ForceSelect,4,1)
        
        self.yScaleSelect = pyqt.QComboBox();
        self.yScaleSelect.addItems(['unit', 'kilo', 'mega', 'giga'])
        self.ScaleConstX = self.ScaleConstY = 1
        self.yScaleSelect.currentIndexChanged.connect( self.UpdatePlot )
        self.OptionsLayout.addWidget(pyqt.QLabel("Force Scale"),5,0)
        self.OptionsLayout.addWidget(self.yScaleSelect,5,1)
        
        self.AxisSelect = pyqt.QComboBox();
        self.AxisSelect.addItems(['x', 'y', 'z', 't'])
        self.AxisSelect.currentIndexChanged.connect( self.PlotData )
        self.OptionsLayout.addWidget(pyqt.QLabel("Plot Axis"),6,0)
        self.OptionsLayout.addWidget(self.AxisSelect,6,1)
        
        self.xScaleSelect = pyqt.QComboBox();
        self.xScaleSelect.addItems(['unit','kilo', 'mega', 'giga','milli', 'micro', 'nano', 'pico', 'femto'])
        self.xScaleSelect.currentIndexChanged.connect( self.UpdatePlot )
        self.OptionsLayout.addWidget(pyqt.QLabel("X-Axis Scale"),7,0)
        self.OptionsLayout.addWidget(self.xScaleSelect,7,1)
        
        self.NormaliseXAxis = pyqt.QCheckBox()
        self.NormaliseXAxis.stateChanged.connect(self.PlotData)
        self.OptionsLayout.addWidget(pyqt.QLabel("Normalise X-Axis to RMS"),8,0)
        self.OptionsLayout.addWidget(self.NormaliseXAxis,8,1)
        
        #Set the central position for the field plotting
        self.XCentralPosition = pyqt.QDoubleSpinBox()
        self.YCentralPosition = pyqt.QDoubleSpinBox()
        self.ZCentralPosition = pyqt.QDoubleSpinBox()
        self.XCentralPosition.setRange(-3,3)
        self.XCentralPosition.setSingleStep(0.5)
        self.YCentralPosition.setRange(-3,3)
        self.YCentralPosition.setSingleStep(0.5)
        self.ZCentralPosition.setRange(-3,3)
        self.ZCentralPosition.setSingleStep(0.5)
        self.XCentralPosition.valueChanged.connect(self.PlotData)
        self.YCentralPosition.valueChanged.connect(self.PlotData)
        self.ZCentralPosition.valueChanged.connect(self.PlotData)
        self.OptionsLayout.addWidget(pyqt.QLabel("Central X Position/sigmax:"),9,0)
        self.OptionsLayout.addWidget(pyqt.QLabel("Central Y Position/sigmay:"),10,0)
        self.OptionsLayout.addWidget(pyqt.QLabel("Central Z Position/sigmaz:"),11,0)
        self.OptionsLayout.addWidget(self.XCentralPosition,9,1)
        self.OptionsLayout.addWidget(self.YCentralPosition,10,1)
        self.OptionsLayout.addWidget(self.ZCentralPosition,11,1)
        
        
        
        
        self.XLabel = pyqt.QLineEdit()
        self.XLabel.setPlaceholderText('X Axis Label')
        self.XLabel.textChanged.connect(self.XAxisChanged)
        self.OptionsLayout.addWidget(pyqt.QLabel("X-Axis Label"),12,0)
        self.OptionsLayout.addWidget(self.XLabel,12,1)
        
        self.YLabel = pyqt.QLineEdit()
        self.YLabel.setPlaceholderText('Y Axis Label')
        self.YLabel.textChanged.connect(self.YAxisChanged)
        self.OptionsLayout.addWidget(pyqt.QLabel("Y-Axis Label"),13,0)
        self.OptionsLayout.addWidget(self.YLabel,13,1)
    
        self.graphwidget = pg.PlotWidget();
        self.xdata = []
        self.ydata = []
        self.graphwidget.setBackground('w')
        self.graphwidget.getAxis('bottom').setPen('k')
        self.graphwidget.getAxis('left').setPen('k')
        self.graphwidget.getAxis('bottom').setTextPen('k')
        self.graphwidget.getAxis('left').setTextPen('k')
        self.graphlegend = self.graphwidget.addLegend(labelTextColor='k')
        self.plot = self.graphwidget.plot(self.xdata,self.ydata);
        self.graphwidget.addItem(self.plot)
        
        #toolbar = NavigationToolbar(self.graphwidget, self)
        
        
        self.TotalLayout.addLayout(self.OptionsLayout)
        self.TotalLayout.addWidget(self.graphwidget)
        self.setLayout(self.TotalLayout)
        
    def BeamLoad(self):
        filename, ok = pyqt.QFileDialog.getOpenFileName(self,"Select a File", "", "HDF5 Files (*.h5)")
        if filename:
            self.beams.append(filename)
            self.beamlabels.append(filename)
            self.BeamList.addItem(filename)
        self.PlotData()
    def ChangeLabel(self):
        SelectedIndex = self.BeamList.currentRow()
        if SelectedIndex == -1:
            print('None Selected')
        else:
            if self.BeamLabelText.text() == '':
                self.beamlabels[SelectedIndex] = 'Beam'
            else:
                self.beamlabels[SelectedIndex] = self.BeamLabelText.text()
                self.BeamList.currentItem().setText(self.BeamLabelText.text())
    def RemoveBeam(self):
        SelectedIndex = self.BeamList.currentRow()
        if SelectedIndex == -1:
            print('None Selected')
        else:
            self.BeamList.takeItem(SelectedIndex)
            self.beamlabels.pop(SelectedIndex)
            self.beams.pop(SelectedIndex)
        self.PlotData()
    def XAxisChanged(self,text):
        self.graphwidget.setLabel('bottom',self.XLabel.text())
    def YAxisChanged(self,text):
        self.graphwidget.setLabel('left',self.YLabel.text())
        
        
    #Update the plot using new data
    def PlotData(self):
        self.xdata =[[]]
        self.ydata = [[]]
        self.graphwidget.clear()
        self.graphlegend.clear()
        for k in range(len(self.beams)):
            DWA_Beam = DWA.DiWaCATOutput()
            try:
                DWA_Beam.ReadFromFile(self.beams[k])
            except ValueError:
                print("ERROR: File Error. Check file path/extensions.")
                return
            XAxis = self.AxisSelect.currentText()
            YAxis = self.ForceSelect.currentIndex()
            FieldFunctions = DWA_Beam.ReturnFieldFunctions()
            if XAxis == 't':
                PlotValues = np.unique(DWA_Beam._FieldPoints['z'])
            else:
                PlotValues = np.unique(DWA_Beam._FieldPoints[XAxis])
                
            #Change this to be for each axis
            XValue = self.XCentralPosition.value() * DWA_Beam.beam.Sx
            YValue = self.YCentralPosition.value() * DWA_Beam.beam.Sy
            ZValue = self.ZCentralPosition.value() * DWA_Beam.beam.Sz
            
            InterpolationPoints = [];
            if XAxis == 'z' or XAxis == 't':
                for i in range(len(PlotValues)):
                    InterpolationPoints.append([XValue,YValue,PlotValues[i]])
            if XAxis == 'y':
                for i in range(len(PlotValues)):
                    InterpolationPoints.append([XValue,PlotValues[i],ZValue])
            if XAxis == 'x':
                for i in range(len(PlotValues)):
                    InterpolationPoints.append([PlotValues[i],YValue,ZValue])
                    
            self.xdata.append(np.array(PlotValues))
            self.ydata.append(np.array(FieldFunctions[YAxis](InterpolationPoints)))
            #self.graphwidget.removeItem(self.plot)
            if self.NormaliseXAxis.isChecked() == False and XAxis == 't':
                self.xdata[k+1] *= 1/const.c
            if self.NormaliseXAxis.isChecked() == True:
                if XAxis == 'x':
                    self.xdata[k+1] *= 1/DWA_Beam.beam.Sx
                if XAxis == 'y':
                    self.xdata[k+1] *= 1/DWA_Beam.beam.Sy
                if XAxis == 'z' or XAxis == 't':
                    self.xdata[k+1] *= 1/DWA_Beam.beam.Sz
            self.plot = self.graphwidget.plot(self.xdata[k+1] * self.ScaleConstX,self.ydata[k+1] * self.ScaleConstY, pen = colors[k]);
            self.graphwidget.addItem(self.plot)
            if len(self.beams) > 1:
                self.graphlegend.addItem(self.plot,self.beamlabels[k])
            
    def UpdatePlot(self):
        XAxisScale = self.xScaleSelect.currentText()
        YAxisScale = self.yScaleSelect.currentText()
        self.ScaleConstX = UnitConversion(XAxisScale)
        self.ScaleConstY = UnitConversion(YAxisScale)
        self.graphwidget.clear()
        for i in range(len(self.beams)):
            self.plot = self.graphwidget.plot(self.xdata[i+1] * self.ScaleConstX,self.ydata[i+1]* self.ScaleConstY, pen = colors[i]);
            self.graphwidget.addItem(self.plot)
            if len(self.beams) > 1:
                self.graphlegend.addItem(self.plot,self.beamlabels[i])
        
class BeamPlotWindow(QWidget):
    """
    Load in a beam - plot or print variables, drift beam
    """
    
    def __init__(self):
        super().__init__()
        
        self.setWindowTitle("Plot Beam Properties")
        windowWidth = 1000
        windowHeight = 600
        self.setFixedSize(windowWidth, windowHeight)
        
        self.TotalLayout = pyqt.QGridLayout()
        self.BeamsLayout = pyqt.QGridLayout()
        self.PropertyPrintLayout = pyqt.QGridLayout()
        self.BeamManipulateLayout = pyqt.QGridLayout()
        self.BeamPlotLayout = pyqt.QGridLayout()
        
        self.beams = [];
        self.beamlabel = [];
        #Loading the beams, changing labels and removing beams
        BeamsLayoutBackground = pyqt.QWidget(self)
        BeamsLayoutBackground.setStyleSheet('background-color : white')
        self.BeamsLayout.addWidget(BeamsLayoutBackground,0,0,4,2)
        
        self.LoadBeam = pyqt.QPushButton('Load Beam')
        self.LoadBeam.clicked.connect(self.BeamLoad)
        self.LoadBeam.setFixedWidth(0.25*windowWidth)
        self.BeamsLayout.addWidget(self.LoadBeam,0,0,1,2)
        
        self.BeamList = pyqt.QListWidget()
        self.BeamList.addItems(self.beams)
        self.BeamList.setFixedWidth(0.25*windowWidth)
        self.BeamList.setMinimumHeight(0.2*windowHeight)
        self.BeamsLayout.addWidget(self.BeamList,1,0,1,2)
        
        self.BeamLabelText = pyqt.QLineEdit()
        self.BeamLabelText.setPlaceholderText('Beam Label')
        self.BeamLabelText.setFixedWidth(0.12*windowWidth)
        self.BeamsLayout.addWidget(self.BeamLabelText,2,0)
        
        self.ChangeLabelButton = pyqt.QPushButton('Change Label')
        self.ChangeLabelButton.clicked.connect(self.ChangeLabel)
        self.ChangeLabelButton.setFixedWidth(0.12*windowWidth)
        self.BeamsLayout.addWidget(self.ChangeLabelButton,2,1)
        
        self.RemoveBeamButton = pyqt.QPushButton('Remove Beam')
        self.RemoveBeamButton.clicked.connect(self.RemoveBeam)
        self.RemoveBeamButton.setFixedWidth(0.25*windowWidth)
        self.BeamsLayout.addWidget(self.RemoveBeamButton,3,0,1,2)
        
        #Beam Manipulation: Drift beams a given distance
        self.DriftDistance = pyqt.QDoubleSpinBox()
        self.DriftDistance.setRange(-1000,1000)
        self.DriftScale = pyqt.QComboBox()
        self.DriftScale.addItems(['unit','kilo', 'mega', 'giga','milli', 'micro', 'nano', 'pico', 'femto'])
        self.BeamManipulateLayout.addWidget(pyqt.QLabel('Drift Distance [m]: '),0,0)
        self.BeamManipulateLayout.addWidget(self.DriftDistance, 0,1)
        self.BeamManipulateLayout.addWidget(self.DriftScale,0,2)
        
        
        #Print Beam Property
        BeamPropertyList = ['Total Charge','RMS Size - X', 'RMS Size - Y','RMS Size - Z', 'RMS Size - T', 'RMS Momentum Spread', 'Emittance - X', 'Emittance - Y', 
                            'Mean Position - X','Mean Position - Y','Mean Position - Z', 'Alpha - X', 'Alpha - Y', 'Beta - X', 'Beta - Y',
                            'Gamma - X', 'Gamma - Y', '95% RMS Size - X', '95% RMS Size - Y', '95% RMS Size - Z','95% RMS Size - T',
                            '90% Emittance - X', '90% Emittance - Y', 'Number of Macros']
        self.BeamProperties = ['total_charge', 'Sx', 'Sy', 'Sz', 'Sz * 1/const.c', 'Spz', 'normalized_horizontal_emittance', 'normalized_vertical_emittance',
                               'Mx', 'My', 'Mz', 'alpha_x', 'alpha_y', 'beta_x', 'beta_y', 'gamma_x', 'gamma_y', 'Sx_95', 'Sy_95',
                               'Sz_95','Sz_95 * 1/const.c', 'normalized_horizontal_emittance_90', 'normalized_vertical_emittance_90',
                               'nMacros']
        self.BeamPropertySelect = pyqt.QListWidget()
        self.BeamPropertySelect.addItems(BeamPropertyList)
        self.BeamPropertySelect.setFixedSize(0.15*windowWidth, 0.2*windowHeight)
        self.PropertyPrintLayout.addWidget(self.BeamPropertySelect, 0, 0, 2, 1)
        self.PrintPropertyScale = pyqt.QComboBox()
        self.PrintPropertyScale.addItems(['unit','kilo', 'mega', 'giga','milli', 'micro', 'nano', 'pico', 'femto'])
        self.PrintPropertyScale.setFixedWidth(0.09*windowWidth)
        self.PropertyPrintLayout.addWidget(self.PrintPropertyScale,0,1,1,1)
        
        self.PrintPropertyButton = pyqt.QPushButton('Print Property')
        self.PrintPropertyButton.setFixedWidth(0.09*windowWidth)
        self.PrintPropertyButton.clicked.connect(self.PrintProperty)
        self.PropertyPrintLayout.addWidget(self.PrintPropertyButton,1,1,1,1)
        self.PropertyPrinted = pyqt.QListView()
        
        self.PropertyPrintBox = pyqt.QListWidget()
        self.PropertyPrintBox.setFixedWidth(0.25*windowWidth)
        self.PropertyPrintBox.setMinimumHeight(0.2*windowHeight)
        self.PropertyPrintLayout.addWidget(self.PropertyPrintBox, 2,0,2,2)
        
        #Beam Property Plot
        self.VariableList = ['x','y','z','t','px','py','pz','cpx','cpy','cpz','xp','yp']
        #First set up the axes
        self.XAxisSelect = pyqt.QComboBox()
        self.XAxisSelect.addItems(self.VariableList)
        self.XAxisScale = pyqt.QComboBox()
        self.XAxisScale.addItems(['unit','kilo', 'mega', 'giga','milli', 'micro', 'nano', 'pico', 'femto'])
        self.BeamPlotLayout.addWidget(pyqt.QLabel('X Axis:'),0,0,1,1)
        self.BeamPlotLayout.addWidget(self.XAxisSelect,0,1,1,1)
        self.BeamPlotLayout.addWidget(self.XAxisScale,0,2,1,1)
        self.XAxisLabel = pyqt.QLineEdit()
        self.BeamPlotLayout.addWidget(pyqt.QLabel('X Axis Label:'),1,0,1,1)
        self.BeamPlotLayout.addWidget(self.XAxisLabel,1,1,1,2)
        
        self.YAxisSelect = pyqt.QComboBox()
        self.YAxisSelect.addItems(self.VariableList)
        self.YAxisScale = pyqt.QComboBox()
        self.YAxisScale.addItems(['unit','kilo', 'mega', 'giga','milli', 'micro', 'nano', 'pico', 'femto'])
        self.BeamPlotLayout.addWidget(pyqt.QLabel('Y Axis:'),2,0,1,1)
        self.BeamPlotLayout.addWidget(self.YAxisSelect,2,1,1,1)
        self.BeamPlotLayout.addWidget(self.YAxisScale,2,2,1,1)
        self.YAxisLabel = pyqt.QLineEdit()
        self.BeamPlotLayout.addWidget(pyqt.QLabel('Y Axis Label:'),3,0,1,1)
        self.BeamPlotLayout.addWidget(self.YAxisLabel,3,1,1,2)
        
        self.nHistBins = pyqt.QSpinBox()
        self.nHistBins.setRange(10,1000)
        self.nHistBins.setSingleStep(10)
        self.nHistBins.setValue(100)
        self.BeamPlotLayout.addWidget(pyqt.QLabel('Number of Histogram Bins:'),4,0,1,2)
        self.BeamPlotLayout.addWidget(self.nHistBins,4,2,1,1)
        
        self.PlotOpacity = pyqt.QDoubleSpinBox()
        self.PlotOpacity.setRange(0,100)
        self.PlotOpacity.setValue(80)
        self.PlotOpacity.setSingleStep(10)
        self.BeamPlotLayout.addWidget(pyqt.QLabel('Scatter Point Opacity [%]'),5,0,1,2)
        self.BeamPlotLayout.addWidget(self.PlotOpacity,5,2,1,1)
        
        self.UpdatePlotButton = pyqt.QPushButton('Update Plots')
        self.UpdatePlotButton.clicked.connect(self.UpdatePlots)
        self.BeamPlotLayout.addWidget(self.UpdatePlotButton,6,0,1,3)
        
        #Create the plots
        
        self.XAxisHistogram = pg.PlotWidget()
        self.XAxisHistogram.setBackground('w')
        self.XAxisHistogram.setLabel('left','Intensity [arb.]')
        self.XAxisHistogram.setFixedSize(windowWidth*0.48, windowHeight*0.3)
        self.XHistValues = [[[],[0]]]
        #Remember this one is rotated!
        self.XAxisHist = pg.PlotCurveItem( self.XHistValues[0][1],self.XHistValues[0][0], stepMode=True, fillLevel=0)
        self.XAxisHistogram.addItem(self.XAxisHist)
        
        self.YAxisHistogram = pg.PlotWidget()
        self.YAxisHistogram.setBackground('w')
        self.YAxisHistogram.setLabel('bottom','Intensity [arb.]')
        self.YAxisHistogram.setFixedSize(windowWidth*0.24, windowHeight*0.65)
        self.YHistValues = [[[],[0]]]
        #Remember this one is rotated!
        self.YAxisHist = pg.PlotCurveItem( self.YHistValues[0][1],self.YHistValues[0][0], stepMode=True, fillLevel=0)
        self.YAxisHistogram.addItem(self.YAxisHist)
        
        self.ScatterPlot = pg.PlotWidget()
        self.ScatterPlot.setBackground('w')
        self.ScatterPlot.setFixedSize(windowWidth*0.48, windowHeight*0.65)
        self.XData = []
        self.YData = []
        self.XYScatter = pg.ScatterPlotItem()
        self.XYScatter.setData(self.XData,self.YData)
        self.ScatterPlot.addItem(self.XYScatter)
        self.legend = self.ScatterPlot.addLegend(frame=True, brush = 'w',labelTextColor='k')
        
        
        self.XAxisHistogram.sigRangeChanged.connect(self.AxisScroll)
        self.YAxisHistogram.sigRangeChanged.connect(self.AxisScroll)
        self.ScatterPlot.sigRangeChanged.connect(self.AxisScroll)
        
        self.AllPlots = [self.ScatterPlot, self.XAxisHistogram, self.YAxisHistogram]
        for plot_item in self.AllPlots:
            plot_item.getAxis('bottom').setPen('k')
            plot_item.getAxis('left').setPen('k')
            plot_item.getAxis('bottom').setTextPen('k')
            plot_item.getAxis('left').setTextPen('k')
            
        
        self.BeamPlotLayout.addWidget(self.XAxisHistogram,0,4,7,6)
        self.BeamPlotLayout.addWidget(self.YAxisHistogram,7,0,14,3)
        self.BeamPlotLayout.addWidget(self.ScatterPlot,7,4,14,6)
        
        
        self.TotalLayout.addLayout(self.BeamsLayout,0,0)
        self.TotalLayout.addLayout(self.BeamManipulateLayout,1,0)
        self.TotalLayout.addLayout(self.PropertyPrintLayout,2,0)
        self.TotalLayout.addLayout(self.BeamPlotLayout,0,1,3,1)
        
        self.setLayout(self.TotalLayout)
        
        
        
    def BeamLoad(self):
        filename, ok = pyqt.QFileDialog.getOpenFileName(self,"Select a File", "", "HDF5 Files (*.h5)")
        if filename:
            self.beams.append(filename)
            self.beamlabel.append(filename)
            self.BeamList.addItem(filename)
    def ChangeLabel(self):
        SelectedIndex = self.BeamList.currentRow()
        if SelectedIndex == -1:
            print('None Selected')
        else:
            if self.BeamLabelText.text() == '':
                self.beamlabel[SelectedIndex] = 'Beam'
            else:
                self.beamlabel[SelectedIndex] = self.BeamLabelText.text()
                self.BeamList.currentItem().setText(self.BeamLabelText.text())
    def RemoveBeam(self):
        SelectedIndex = self.BeamList.currentRow()
        if SelectedIndex == -1:
            print('None Selected')
        else:
            self.BeamList.takeItem(SelectedIndex)
            self.beamlabel.pop(SelectedIndex)
            self.beams.pop(SelectedIndex)
            
    def PrintProperty(self):
        self.PropertyPrintBox.clear()
        SelectedProperty = self.BeamProperties[self.BeamPropertySelect.currentRow()]
        print(SelectedProperty)
        if self.BeamPropertySelect.currentRow() == -1:
            self.PropertyPrintBox.addItem('No Property Selected')
        elif len(self.beams)==0:
            self.PropertyPrintBox.addItem('No Beams Loaded')
        else:
            for i in range(len(self.beams)):
                loadedBeam = dbt.BeamFromFile()
                loadedBeam.read_WakeCode_beam_file(self.beams[i])
                loadedBeam.driftBeam(self.DriftDistance.value() * 1/UnitConversion(self.DriftScale.currentText()))
                if SelectedProperty == 'Spz':
                    value = np.sqrt(loadedBeam.covariance(loadedBeam.cpz, loadedBeam.cpz)) * UnitConversion(self.PrintPropertyScale.currentText())
                else:
                    value = eval("loadedBeam."+SelectedProperty) * UnitConversion(self.PrintPropertyScale.currentText())
                self.PropertyPrintBox.addItem( self.beamlabel[i] + ":  " + str(value))
        
    def AxisScroll(self,r):
        self.XAxisHistogram.blockSignals(True)
        self.YAxisHistogram.blockSignals(True)
        self.ScatterPlot.blockSignals(True)
        
        if r == self.XAxisHistogram:
            self.ScatterPlot.setRange(xRange=r.getAxis('bottom').range)
        if r == self.YAxisHistogram:
            self.ScatterPlot.setRange(yRange=r.getAxis('left').range)
        if r == self.ScatterPlot:
            self.XAxisHistogram.setRange(xRange=self.ScatterPlot.getAxis('bottom').range)
            self.YAxisHistogram.setRange(yRange=self.ScatterPlot.getAxis('left').range)
        self.XAxisHistogram.blockSignals(False)
        self.YAxisHistogram.blockSignals(False)
        self.ScatterPlot.blockSignals(False)
                
    def UpdatePlots(self):
        if len(self.beams)==0:
            print('No Beams Loaded')
        else:
            self.XData = []
            self.YData = []
            self.XHistValues = []
            self.YHistValues = []
            for plot_item in self.AllPlots:
                plot_item.clear()
            self.legend.clear()
            for i in range(len(self.beams)):
                loadedBeam = dbt.BeamFromFile()
                loadedBeam.read_WakeCode_beam_file(self.beams[i])
                loadedBeam.driftBeam(self.DriftDistance.value() * 1/UnitConversion(self.DriftScale.currentText()))
                self.XData = eval("loadedBeam." + self.XAxisSelect.currentText()) * UnitConversion(self.XAxisScale.currentText())
                self.YData = eval("loadedBeam." + self.YAxisSelect.currentText()) * UnitConversion(self.YAxisScale.currentText())
                xyscatter = pg.ScatterPlotItem(size = 1, pen = colors[i], brush = colors[i])
                #Make the same scatter plot but with bigger points for the legend
                LargePointScatter = pg.ScatterPlotItem(size = 5, pen = colors[i], brush = colors[i])
                xyscatter.setData(self.XData,self.YData)
                xyscatter.setOpacity(self.PlotOpacity.value()*0.01)
                self.ScatterPlot.addItem(xyscatter, name = self.beamlabel[i])
                if len(self.beams) > 1:
                    self.legend.addItem(LargePointScatter,self.beamlabel[i])
                self.XHistValues = np.histogram(self.XData, bins = np.linspace(min(self.XData), max(self.XData), self.nHistBins.value()))
                self.YHistValues = np.histogram(self.YData, bins = np.linspace(min(self.YData), max(self.YData), self.nHistBins.value()))
                self.XAxisHist = pg.PlotCurveItem(self.XHistValues[1], self.XHistValues[0]*1/max(self.XHistValues[0]), stepMode=True, fillLevel=0, pen = colors[i])
                self.YAxisHist = pg.PlotCurveItem(self.YHistValues[1], self.YHistValues[0]*1/max(self.YHistValues[0]), stepMode=True, fillLevel=0, pen = colors[i])
                self.YAxisHist.setRotation(90)
                self.XAxisHistogram.addItem(self.XAxisHist)
                self.YAxisHistogram.addItem(self.YAxisHist)
                #self.YAxisHistogram.getPlotItem(self.YAxisHist).invertX(True)
                
                
            self.ScatterPlot.addLegend()
            self.ScatterPlot.setLabel('bottom', self.XAxisLabel.text())
            self.ScatterPlot.setLabel('left', self.YAxisLabel.text())
            self.XAxisHistogram.setLabel('bottom', self.XAxisLabel.text())
            self.YAxisHistogram.setLabel('left', self.YAxisLabel.text())


def WorkerFunction_Interpolation(beam, position):
    return beam.ReturnFieldInterpolated(position)


class BeamTrackFunction(QThread):
    progress = pyqtSignal(float)
    stepdone = pyqtSignal()
    finished = pyqtSignal()
    newbeam = pyqtSignal(object)
    
    def __init__(self, BeamFile, FieldFile, DLWLength, nSteps):
        QThread.__init__(self)
        self.BeamFile = BeamFile
        self.FieldFile = FieldFile
        self.DLWLength = DLWLength
        self.nSteps = nSteps
    
    def run(self):
        #self.ProgressBar.setRange(0, self.nSteps.value())
        #self.ProgressBar.setValue(0)
        self.progress.emit(0)
        InitialBeam = DWA.DiWaCATOutput()
        InitialBeam.ReadFromFile(self.BeamFile)
        InitialBeam.GetFieldValues(self.FieldFile)
        StepLength = self.DLWLength * 1/self.nSteps
        InitialBeam.SetDielectricLength(StepLength)   
        CurrentBeam = InitialBeam
        self.newbeam.emit(CurrentBeam)
        for i in range(self.nSteps):
            #First Drift the Beam by 1/2 a step
            CurrentBeam.driftBeam(0.5 * StepLength)
            MacroPositions = np.stack((CurrentBeam.beam.x,CurrentBeam.beam.y,CurrentBeam.beam.z),axis = -1)
            #Interpolate the field at the new positions
            KicksApplied = np.empty((len(MacroPositions),3))
            #ForceFunctions = CurrentBeam.ReturnFieldFunctions()
            for j in range(len(KicksApplied)):
                #KicksApplied[i] = np.empty(3)
                KicksApplied[j] = CurrentBeam.ReturnFieldInterpolated(MacroPositions[j])
                if (j%1000 == 0):
                    self.progress.emit((i + j/len(KicksApplied))*1/self.nSteps)
                    
            #Drift the beam to the end of the step and then apply the kicks
            CurrentBeam.driftBeam(0.5 * StepLength)
            CurrentBeam.beam._beam['px'] += StepLength * KicksApplied[:,0] * CurrentBeam.beam.q_over_c
            CurrentBeam.beam._beam['py'] += StepLength * KicksApplied[:,1] * CurrentBeam.beam.q_over_c
            CurrentBeam.beam._beam['pz'] += StepLength * KicksApplied[:,2] * CurrentBeam.beam.q_over_c
            
            #Clear lost particles
            if CurrentBeam._DielectricParameter['w'] == 0:
                Boundary = CurrentBeam._DielectricParameter['a']
                RadialPositions = np.sqrt((CurrentBeam.beam.x+CurrentBeam._DielectricParameter['x0'])**2 + (CurrentBeam.beam.y+CurrentBeam._DielectricParameter['y0'])**2)
                LostParticles = [n for n,i in enumerate(RadialPositions) if i>Boundary]
                CurrentBeam.beam._beam['charge'][LostParticles] = 0;
            else:
                LostParticlesY = [n for n,i in enumerate(CurrentBeam.beam._beam['y']) if abs(i+CurrentBeam._DielectricParameter['y0']) > CurrentBeam._DielectricParameter['a']]
                LostParticlesX = [n for n,i in enumerate(CurrentBeam.beam._beam['x']) if abs(i+CurrentBeam._DielectricParameter['x0']) > CurrentBeam._DielectricParameter['w']]
                CurrentBeam.beam._beam['charge'][LostParticlesY] = 0;
                CurrentBeam.beam._beam['charge'][LostParticlesX] = 0;
            macro_select = np.nonzero(CurrentBeam.beam._beam['charge'])
            
            CurrentBeam.beam._beam['x'] = CurrentBeam.beam._beam['x'][macro_select]
            CurrentBeam.beam._beam['y'] = CurrentBeam.beam._beam['y'][macro_select]
            CurrentBeam.beam._beam['z'] = CurrentBeam.beam._beam['z'][macro_select]
            CurrentBeam.beam._beam['px'] = CurrentBeam.beam._beam['px'][macro_select]
            CurrentBeam.beam._beam['py'] = CurrentBeam.beam._beam['py'][macro_select]
            CurrentBeam.beam._beam['pz'] = CurrentBeam.beam._beam['pz'][macro_select]
            CurrentBeam.beam._beam['t'] = CurrentBeam.beam._beam['t'][macro_select]
            CurrentBeam.beam._beam['charge'] = CurrentBeam.beam._beam['charge'][macro_select]
            
            self.progress.emit((i+1)*1/self.nSteps)
            self.newbeam.emit(CurrentBeam)
            
        self.finished.emit()
            
class BeamTrackWindow(QWidget):
    '''
    Track the beam through a DLW - ability to mesh and save beam too
    '''
    def __init__(self):
        super().__init__()
        
        self.setWindowTitle("Tracking within DLW")
        self.OptionsLayout = pyqt.QGridLayout()
        self.PlotsLayout = pyqt.QGridLayout()
        self.TotalLayout = pyqt.QGridLayout()
        
        self.LoadBeamButton = pyqt.QPushButton('Load Beam File')
        self.LoadBeamButton.clicked.connect(self.LoadBeamFile)
        self.BeamFile = pyqt.QLineEdit()
        self.BeamFile.setPlaceholderText('Beam File Name')
        self.BeamFieldSame = pyqt.QPushButton('Beam File = Field File')
        self.BeamFieldSame.clicked.connect(self.BeamFieldFileSame)
        self.LoadFieldButton = pyqt.QPushButton('Load Field File')
        self.LoadFieldButton.clicked.connect(self.LoadFieldFile)
        self.FieldFile = pyqt.QLineEdit()
        self.FieldFile.setPlaceholderText('Field File Name')
        self.OptionsLayout.addWidget(self.LoadBeamButton,0,0,1,3)
        self.OptionsLayout.addWidget(self.BeamFile,1,0,1,3)
        self.OptionsLayout.addWidget(self.BeamFieldSame,2,0,1,3)
        self.OptionsLayout.addWidget(self.LoadFieldButton,3,0,1,3)
        self.OptionsLayout.addWidget(self.FieldFile,4,0,1,3)
        
        self.DLWLength = pyqt.QDoubleSpinBox()
        self.DLWLength.setRange(0,1000)
        self.DLWLengthScale = pyqt.QComboBox()
        self.DLWLengthScale.addItems(['unit','kilo', 'mega', 'giga','milli', 'micro', 'nano', 'pico', 'femto'])
        self.OptionsLayout.addWidget(pyqt.QLabel('DLW Length [m]: '), 5,0,1,1)
        self.OptionsLayout.addWidget(self.DLWLength, 5,1,1,1)
        self.OptionsLayout.addWidget(self.DLWLengthScale, 5,2,1,1)
        
        self.nSteps = pyqt.QSpinBox()
        self.nSteps.setRange(1,100)
        self.nSteps.setValue(3)
        self.OptionsLayout.addWidget(pyqt.QLabel('Number of Steps: '),6,0,1,1)
        self.OptionsLayout.addWidget(self.nSteps,6,1,1,1)
        
        self.DriftDistance = pyqt.QDoubleSpinBox()
        self.DriftDistance.setRange(0,1000)
        self.DriftDistanceScale = pyqt.QComboBox()
        self.DriftDistanceScale.addItems(['unit','kilo', 'mega', 'giga','milli', 'micro', 'nano', 'pico', 'femto'])
        self.OptionsLayout.addWidget(pyqt.QLabel('Post DLW Drift Distance [m]'), 7, 0)
        self.OptionsLayout.addWidget(self.DriftDistance, 7, 1)
        self.OptionsLayout.addWidget(self.DriftDistanceScale, 7, 2)
        
        self.SimulateButton = pyqt.QPushButton('Run Simulation')
        #self.SimulateButton.clicked.connect(self.TrackBeam)
        self.SimulateButton.clicked.connect(self.RunSim)
        self.OptionsLayout.addWidget(self.SimulateButton,8,0,1,3)
        self.ProgressBar = pyqt.QProgressBar()
        self.ProgressBar.setRange(0, 100)
        self.ProgressBar.setValue(0)
        self.OptionsLayout.addWidget(self.ProgressBar,9,0,1,3)
        
        self.BeamList = pyqt.QListWidget()
        self.TrackedBeams = []
        self.OptionsLayout.addWidget(pyqt.QLabel('Select Tracked Beam:'),10,0,1,1)
        self.OptionsLayout.addWidget(self.BeamList,11,0,1,3)
        
        self.MeshBeam = pyqt.QCheckBox('Mesh Simulated Beam')
        self.MeshBeam.stateChanged.connect(self.MeshVariables)
        self.OptionsLayout.addWidget(self.MeshBeam,12,0)
        
        self.SaveBeam = pyqt.QPushButton('Save Selected Beam')
        self.SaveBeam.clicked.connect(self.SaveSingleBeam)
        self.OptionsLayout.addWidget(self.SaveBeam,19,0,1,3)
        
        self.SaveAllBeamsButton = pyqt.QPushButton('Save All Beams')
        self.SaveAllBeamsButton.clicked.connect(self.SaveAllBeams)
        self.OptionsLayout.addWidget(self.SaveAllBeamsButton,20,0,1,3)
        
        #Plot a handful of variables - set up the plots for when tracking
        self.LengthArray = []
        self.SigmaX = []
        self.SigmaY = []
        self.Mx = []
        self.My = []
        self.epsx = []
        self.epsy = []
        self.Charge = []
        self.SigmaT = []
        
        self.SigXPlot = pg.PlotWidget()
        self.SigYPlot = pg.PlotWidget()
        self.SigTPlot = pg.PlotWidget()
        self.MXPlot = pg.PlotWidget()
        self.MYPlot = pg.PlotWidget()
        self.epsXPlot = pg.PlotWidget()
        self.epsYPlot = pg.PlotWidget()
        self.ChargePlot = pg.PlotWidget()
        
        self.Plots = [self.SigXPlot,self.SigYPlot,self.SigTPlot,self.MXPlot,self.MYPlot,self.epsXPlot,self.epsYPlot,self.ChargePlot]
        PlotAxisLabel = ['sigmax [m]','sigmay [m]','sigmat [s]','Mean Position - X [m]','Mean Position - Y [m]','Emittance - X [m rad]',
                         'Emittance - Y [m rad]','Charge Transported [%]']
        for i in range(len(self.Plots)):
            self.Plots[i].setBackground('w')
            self.Plots[i].getAxis('bottom').setPen('k')
            self.Plots[i].getAxis('left').setPen('k')
            self.Plots[i].getAxis('bottom').setTextPen('k')
            self.Plots[i].getAxis('left').setTextPen('k')
            self.Plots[i].setLabel('bottom', 'L [m]')
            self.Plots[i].setLabel('left', PlotAxisLabel[i])
        
        self.SigXPlot.plot(self.LengthArray,self.SigmaX, pen = 'k')
        self.SigYPlot.plot(self.LengthArray,self.SigmaY, pen = 'k')
        self.SigTPlot.plot(self.LengthArray,self.SigmaT, pen = 'k')
        self.MXPlot.plot(self.LengthArray,self.Mx, pen = 'k')
        self.MYPlot.plot(self.LengthArray,self.My, pen = 'k')
        self.epsXPlot.plot(self.LengthArray,self.epsx, pen = 'k')
        self.epsYPlot.plot(self.LengthArray,self.epsy, pen = 'k')
        self.ChargePlot.plot(self.LengthArray,self.Charge, pen = 'k')
        
        self.PlotsLayout.addWidget(self.SigXPlot,0,0)
        self.PlotsLayout.addWidget(self.SigYPlot,0,1)
        self.PlotsLayout.addWidget(self.MXPlot,1,0)
        self.PlotsLayout.addWidget(self.MYPlot,1,1)
        self.PlotsLayout.addWidget(self.epsXPlot,2,0)
        self.PlotsLayout.addWidget(self.epsYPlot,2,1)
        self.PlotsLayout.addWidget(self.SigTPlot,3,0)
        self.PlotsLayout.addWidget(self.ChargePlot,3,1)
        
        self.TotalLayout.addLayout(self.OptionsLayout,0,0,1,1)
        self.TotalLayout.addLayout(self.PlotsLayout,0,1,1,2)
        self.setLayout(self.TotalLayout)
        
    def LoadBeamFile(self):
        filename, ok = pyqt.QFileDialog.getOpenFileName(self,"Select a File", "", "HDF5 Files (*.h5)")
        self.BeamFile.setText(filename)
    def LoadFieldFile(self):
        filename, ok = pyqt.QFileDialog.getOpenFileName(self,"Select a File", "", "HDF5 Files (*.h5)")
        self.FieldFile.setText(filename)
    def BeamFieldFileSame(self):
        self.FieldFile.setText(self.BeamFile.text())
    def MeshVariables(self):
        if self.MeshBeam.isChecked():
            self.DefaultMesh = pyqt.QPushButton('Default Mesh Parameters')
            self.DefaultMesh.clicked.connect(self.DefaultMeshPressed)
            self.OptionsLayout.addWidget(self.DefaultMesh, 13,0,1,3)
            
            self.CellsPerSigmaT = pyqt.QDoubleSpinBox()
            self.CellsPerSigmaTLabel =pyqt.QLabel('Cells per Sigma (Transverse):')
            self.CellsPerSigmaT.setMinimum(1.5)
            self.CellsPerSigmaT.setSingleStep(0.5)
            self.CellsPerSigmaT.setValue(3)
            self.OptionsLayout.addWidget(self.CellsPerSigmaTLabel,14,0,1,1)
            self.OptionsLayout.addWidget(self.CellsPerSigmaT,14,1,1,2)
             
            self.CellsPerSigmaL = pyqt.QDoubleSpinBox()
            self.CellsPerSigmaLLabel =pyqt.QLabel('Cells per Sigma (Longitudinal):')
            self.CellsPerSigmaL.setMinimum(1.5)
            self.CellsPerSigmaL.setSingleStep(0.5)
            self.CellsPerSigmaL.setValue(3.5)
            self.OptionsLayout.addWidget(self.CellsPerSigmaLLabel,15,0,1,1)
            self.OptionsLayout.addWidget(self.CellsPerSigmaL,15,1,1,2)
            
            self.MaxXCell = pyqt.QDoubleSpinBox()
            self.MaxXCellLabel = pyqt.QLabel('Maximum X Cell Position / sigmax:')
            self.MaxXCell.setMinimum(1)
            self.MaxXCell.setSingleStep(0.25)
            self.MaxXCell.setValue(2.25)
            self.OptionsLayout.addWidget(self.MaxXCellLabel,16,0,1,1)
            self.OptionsLayout.addWidget(self.MaxXCell,16,1,1,2)
            
            self.MaxYCell = pyqt.QDoubleSpinBox()
            self.MaxYCellLabel = pyqt.QLabel('Maximum Y Cell Position / sigmay:')
            self.MaxYCell.setMinimum(1)
            self.MaxYCell.setSingleStep(0.25)
            self.MaxYCell.setValue(2.25)
            self.OptionsLayout.addWidget(self.MaxYCellLabel,17,0,1,1)
            self.OptionsLayout.addWidget(self.MaxYCell,17,1,1,2)
            
            self.MaxZCell = pyqt.QDoubleSpinBox()
            self.MaxZCellLabel = pyqt.QLabel('Maximum Z Cell Position / sigmaz:')
            self.MaxZCell.setMinimum(1)
            self.MaxZCell.setSingleStep(0.25)
            self.MaxZCell.setValue(3)
            self.OptionsLayout.addWidget(self.MaxZCellLabel,18,0,1,1)
            self.OptionsLayout.addWidget(self.MaxZCell,18,1,1,2)
        else:
            self.DefaultMesh.setVisible(False)
            self.CellsPerSigmaT.setVisible(False)
            self.CellsPerSigmaL.setVisible(False)
            self.MaxXCell.setVisible(False)
            self.MaxYCell.setVisible(False)
            self.MaxZCell.setVisible(False)
            self.CellsPerSigmaTLabel.setVisible(False)
            self.CellsPerSigmaLLabel.setVisible(False)
            self.MaxXCellLabel.setVisible(False)
            self.MaxYCellLabel.setVisible(False)
            self.MaxZCellLabel.setVisible(False)
    def DefaultMeshPressed(self):
        self.CellsPerSigmaT.setValue(3)
        self.CellsPerSigmaL.setValue(3.5)
        self.MaxXCell.setValue(2.25)
        self.MaxYCell.setValue(2.25)
        self.MaxZCell.setValue(3)
        
    def SaveSingleBeam(self):
        if len(self.TrackedBeams)>0 and self.BeamList.currentRow()!=-1:
            filename, ok = pyqt.QFileDialog.getSaveFileName(self,"Select a File", "", "HDF5 Files (*.h5)")
            if filename:
                path = Path(filename)
            beamSelected = self.TrackedBeams[self.BeamList.currentRow()]
            if self.MeshBeam.isChecked():
                Lx = 2* self.MaxXCell.value() * beamSelected.Sx
                Dx = Lx/(2* self.MaxXCell.value() *self.CellsPerSigmaT.value())
                Ly = 2* self.MaxYCell.value() * beamSelected.Sy
                Dy = Ly/(2* self.MaxYCell.value() *self.CellsPerSigmaT.value())
                Lz = 2* self.MaxZCell.value() * beamSelected.Sz
                Dz = Lz/(2* self.MaxZCell.value() *self.CellsPerSigmaL.value())
                
                binnedBeam = dbt.beamBinner(beamSelected)
                binnedBeam.binTheBeam(Lx=Lx,Dx=Dx,Ly=Ly,Dy=Dy,Lz=Lz,Dz=Dz)
                binnedBeam.write_hd5f_mesh_file(path,outputMacros=True, overwrite_existing_file=True,includeSmoothed=False,includeUnSmoothed=True)
                print('File Written')
            else:
                beamSelected.write_HDF5_beam_file(path,overwrite_existing_file=True)
                print('File Written')
        else:
            print('No Beams Tracked or No Beam Selected')
    def SaveAllBeams(self):
        folderpath = pyqt.QFileDialog.getExistingDirectory(self, "Select a Directory")
        if len(self.TrackedBeams)>0:
            for i in range(len(self.TrackedBeams)):
                beamSelected = self.TrackedBeams[i]
                if i==0:
                    path = Path(folderpath + '/BeamIn.h5')
                else:
                    path = Path(folderpath + '/Beam' + str(i) + '.h5')
                if self.MeshBeam.isChecked():
                    Lx = 2* self.MaxXCell.value() * beamSelected.Sx
                    Dx = Lx/(2* self.MaxXCell.value() *self.CellsPerSigmaT.value())
                    Ly = 2* self.MaxYCell.value() * beamSelected.Sy
                    Dy = Ly/(2* self.MaxYCell.value() *self.CellsPerSigmaT.value())
                    Lz = 2* self.MaxZCell.value() * beamSelected.Sz
                    Dz = Lz/(2* self.MaxZCell.value() *self.CellsPerSigmaL.value())
                    
                    binnedBeam = dbt.beamBinner(beamSelected)
                    binnedBeam.binTheBeam(Lx=Lx,Dx=Dx,Ly=Ly,Dy=Dy,Lz=Lz,Dz=Dz)
                    
                    binnedBeam.write_hd5f_mesh_file(path,outputMacros=True, overwrite_existing_file=True,includeSmoothed=False,includeUnSmoothed=True)
                    print('File Written')
                else:
                    beamSelected.write_HDF5_beam_file(path,overwrite_existing_file=True)
                    print('File Written')
        else:
            print('No Beams Tracked or No Beam Selected')
            

    def RunSim(self):
        if self.DLWLength.value() == 0:
            print('DLW Length = 0. No Sim Performed')
        else:
            #Set the files and run button as non-editable
            self.BeamFile.setEnabled(False)
            self.LoadBeamButton.setEnabled(False)
            self.BeamFieldSame.setEnabled(False)
            self.LoadFieldButton.setEnabled(False)
            self.FieldFile.setEnabled(False)
            self.DLWLength.setEnabled(False)
            self.DLWLengthScale.setEnabled(False)
            self.nSteps.setEnabled(False)
            self.DriftDistance.setEnabled(False)
            self.DriftDistanceScale.setEnabled(False)
            self.SimulateButton.setText('Simulation Running')
            self.SimulateButton.setEnabled(False)
            
            clearlist = [self.TrackedBeams, self.BeamList, self.LengthArray, self.SigmaX, self.SigmaY, self.SigmaT, self.Charge, self.Mx, self.My, self.epsx, self.epsy]
            for item in clearlist:
                item.clear()
            
            #Load up the function
            #This class (BeamTrackFunction) starts a new thread on which to run the simulation - then feeds progress and the beam at each step back
            self.function = BeamTrackFunction(self.BeamFile.text(), self.FieldFile.text(), self.DLWLength.value() * 1/UnitConversion(self.DLWLengthScale.currentText()), self.nSteps.value())
            self.function.progress.connect(self.UpdateProgress)
            self.function.newbeam.connect(self.AddNewBeam)
            self.function.finished.connect(self.function.quit)
            self.function.finished.connect(self.SimFinished)
            self.function.start()
        
        
    def UpdateProgress(self,num):
        self.ProgressBar.setValue(round(100*num))
    def AddNewBeam(self,newBeam):
        BeamIndex = len(self.TrackedBeams)
        beamToAdd = dbt.GeneralBeam()
        beamToAdd = newBeam.beam
        self.TrackedBeams.append(copy.deepcopy(beamToAdd))
        if BeamIndex == 0:
            self.BeamList.addItem('Beam In')
        elif BeamIndex == self.nSteps.value():
            self.BeamList.addItem('DLW End')
        else:
            self.BeamList.addItem('Beam ' + str(BeamIndex))
        
        self.LengthArray.append(self.DLWLength.value() * 1/UnitConversion(self.DLWLengthScale.currentText()) * BeamIndex/self.nSteps.value())
        self.SigmaX.append(newBeam.beam.Sx)
        self.SigmaY.append(newBeam.beam.Sy)
        self.Mx.append(newBeam.beam.Mx)
        self.My.append(newBeam.beam.My)
        self.epsx.append(newBeam.beam.normalized_horizontal_emittance)
        self.epsy.append(newBeam.beam.normalized_vertical_emittance)
        if BeamIndex == 0:
            self.Charge.append(100)
        else:    
            self.Charge.append(100 * newBeam.beam.total_charge * 1/self.TrackedBeams[0].total_charge)
        self.SigmaT.append(newBeam.beam.Sz * 1/const.c)
        for plot_item in self.Plots:
            plot_item.clear()
        self.SigXPlot.addItem(self.SigXPlot.plot(self.LengthArray,self.SigmaX, pen = 'k'))
        self.SigYPlot.plot(self.LengthArray,self.SigmaY, pen = 'k')
        self.SigTPlot.plot(self.LengthArray,self.SigmaT, pen = 'k')
        self.MXPlot.plot(self.LengthArray,self.Mx, pen = 'k')
        self.MYPlot.plot(self.LengthArray,self.My, pen = 'k')
        self.epsXPlot.plot(self.LengthArray,self.epsx, pen = 'k')
        self.epsYPlot.plot(self.LengthArray,self.epsy, pen = 'k')
        self.ChargePlot.plot(self.LengthArray,self.Charge, pen = 'k')
        
        
    def SimFinished(self):
        
        self.BeamFile.setEnabled(True)
        self.LoadBeamButton.setEnabled(True)
        self.BeamFieldSame.setEnabled(True)
        self.LoadFieldButton.setEnabled(True)
        self.FieldFile.setEnabled(True)
        self.DLWLength.setEnabled(True)
        self.DLWLengthScale.setEnabled(True)
        self.nSteps.setEnabled(True)
        self.DriftDistance.setEnabled(True)
        self.DriftDistanceScale.setEnabled(True)
        self.SimulateButton.setText('Run Simulation')
        self.SimulateButton.setEnabled(True)
        if self.DriftDistance.value() > 0:
            self.DriftBeamEnd()
    def DriftBeamEnd(self):
        beamToDrift = DWA.DiWaCATOutput()
        beamToDrift.beam = copy.deepcopy(self.TrackedBeams[-1])
        driftDistance = self.DriftDistance.value() * 1/UnitConversion(self.DriftDistanceScale.currentText())
        beamToDrift.driftBeam(driftDistance)
        self.TrackedBeams.append(beamToDrift.beam)
        self.BeamList.addItem('Post Drift')
        
        self.LengthArray.append(self.DLWLength.value() * 1/UnitConversion(self.DLWLengthScale.currentText()) + driftDistance)
        self.SigmaX.append(beamToDrift.Sx)
        self.SigmaY.append(beamToDrift.Sy)
        self.Mx.append(beamToDrift.Mx)
        self.My.append(beamToDrift.My)
        self.epsx.append(beamToDrift.normalized_horizontal_emittance)
        self.epsy.append(beamToDrift.normalized_vertical_emittance)
        self.Charge.append(100 * beamToDrift.total_charge * 1/self.TrackedBeams[0].total_charge)
        self.SigmaT.append(beamToDrift.Sz * 1/const.c)
        for plot_item in self.Plots:
            plot_item.clear()
        self.SigXPlot.addItem(self.SigXPlot.plot(self.LengthArray,self.SigmaX, pen = 'k'))
        self.SigYPlot.plot(self.LengthArray,self.SigmaY, pen = 'k')
        self.SigTPlot.plot(self.LengthArray,self.SigmaT, pen = 'k')
        self.MXPlot.plot(self.LengthArray,self.Mx, pen = 'k')
        self.MYPlot.plot(self.LengthArray,self.My, pen = 'k')
        self.epsXPlot.plot(self.LengthArray,self.epsx, pen = 'k')
        self.epsYPlot.plot(self.LengthArray,self.epsy, pen = 'k')
        self.ChargePlot.plot(self.LengthArray,self.Charge, pen = 'k')
            
class MainWindow(QMainWindow):

    def __init__(self):
        super().__init__()
        
        self.setWindowTitle("DiWaCAT Beam Functions")
        self.setMinimumWidth(350)
        MainBoxLayout = pyqt.QGridLayout();
        
        DiWaCATLoader = pyqt.QPushButton("Open DiWaCAT Field Calculator")
        DiWaCATLoader.setStyleSheet("background-color: rgb(230,85,85)")
        MainBoxLayout.addWidget(DiWaCATLoader,0,0,1,2)
        DiWaCATLoader.clicked.connect(self.DiWaCATOpen)
        
        BeamMakeButton = pyqt.QPushButton("Beam Maker");
        MainBoxLayout.addWidget(BeamMakeButton, 1, 0)
        BeamMakeButton.clicked.connect(self.BeamMakeOpen)
        
        BeamMeshButton = pyqt.QPushButton("Mesh Existing Beam");
        MainBoxLayout.addWidget(BeamMeshButton, 1, 1)
        BeamMeshButton.clicked.connect(self.BeamMeshOpen)
        
        ForcePlotButton = pyqt.QPushButton("Field Plotter");
        MainBoxLayout.addWidget(ForcePlotButton, 2, 0)
        ForcePlotButton.clicked.connect(self.FieldPlotOpen)
        
        BeamPlotButton = pyqt.QPushButton("Beam Plotter");
        MainBoxLayout.addWidget(BeamPlotButton, 2, 1)
        BeamPlotButton.clicked.connect(self.BeamPlotOpen)
        
        SimulateButton = pyqt.QPushButton("Simulate Beam within DLW");
        SimulateButton.setStyleSheet("background-color: rgb(90,205,215)")
        MainBoxLayout.addWidget(SimulateButton, 3, 0,1,2)
        SimulateButton.clicked.connect(self.SimulateBeamWindow)
        
        
        widget = pyqt.QWidget()
        widget.setLayout(MainBoxLayout)
        self.setCentralWidget(widget)

    def BeamMakeOpen(self,checked):
        self.b = BeamMakeWindow()
        self.b.show()
        
    def BeamMeshOpen(self,checked):
        self.c = BeamMeshWindow()
        self.c.show()

    def FieldPlotOpen(self, checked):
        self.d = FieldPlotWindow()
        self.d.show()
        
    def BeamPlotOpen(self,checked):
        self.e = BeamPlotWindow()
        self.e.show()
        
    def SimulateBeamWindow(self,checked):
        self.f = BeamTrackWindow()
        self.f.show()
        
    def DiWaCATOpen(self,checked):
        popen = subprocess.Popen(["./FieldSolver_Executable/DiWaCAT_FieldSolverUI.exe"])
        '''
        popen = subprocess.Popen(["./DWA_Code/DiWaCAT_Executable/DiWaCAT_FieldSolverUI.exe"], stdout=subprocess.PIPE,bufsize=1)
        lines_iterator = iter(popen.stdout.readline, b"")
        while popen.poll() is None:
            for line in lines_iterator:
                nline = line.rstrip()
                print(nline.decode("latin"), end = "\r\n",flush =True)
        '''

if __name__ == "__main__":
    multiprocessing.freeze_support()
    app = QApplication(sys.argv)
    w = MainWindow()
    w.show()
    app.exec()