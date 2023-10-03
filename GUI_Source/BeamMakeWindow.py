# -*- coding: utf-8 -*-
"""
BeamMakeWindow.py
T Overton
03/10/2023

Class containing the beam making GUI functions
Uses functions from beam_tools to make beam from the 6D distribution
"""

from PyQt5.QtWidgets import QWidget
import PyQt5.QtWidgets as pyqt
import pyqtgraph as pg
from pathlib import Path
import numpy as np
import scipy.constants as const

from Python_Tools.Modules import beam_tools as dbt

import matplotlib.pyplot as plt
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']


def UnitConversion(prefix):
    if prefix == 'unit' or prefix == 'Unit':
        return 1
    else:
        return 1/getattr(const, prefix)


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
        
        #Set up widgets for all input parameters
        #Use of a drop box for unit where appropriate
        
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
        
        #For correlated energy spread: positive correlation = higher energy at the tail
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
        
        #For waist position a positive number is a waist downstream
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
        #General sensible mesh density values for planar DLW
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
        
        #The beam is made at the waist
        #Want a positive number to be a downstream waist position so we're pushing the beam backwards
        myBeam._beam['x'] += (-1) * self.HorizontalWaistPosition.value() * myBeam.xp
        myBeam._beam['y'] += (-1) * self.VerticalWaistPosition.value() * myBeam.yp
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