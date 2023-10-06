# -*- coding: utf-8 -*-
"""
Field Calculator GUI
T Overton
5/10/23


Using the Eigenfunction method or Circular DLW method
Fields calculated using C++ functions for speed
"""

from PyQt5.QtWidgets import QWidget
from PyQt5.QtCore import QThread, pyqtSignal,Qt
import PyQt5.QtWidgets as pyqt
import pyqtgraph as pg
from pathlib import Path
import numpy as np
import copy
import scipy.constants as const
from pathos.helpers import mp
import os

from Cpp_Source import DiWaCAT_FieldCalc as DWA

def UnitConversion(prefix):
    if prefix == 'unit' or prefix == 'Unit':
        return 1
    else:
        return 1/getattr(const, prefix)
    
    
class CalculationWorkerFunction(QThread):
    finished = pyqtSignal(object)
    nosim = pyqtSignal()
    def __init__(self,InputList, BeamFile):
        QThread.__init__(self)
        self.InputList=InputList
        self.BeamFile=BeamFile

    def run(self):
        if os.path.isfile(self.BeamFile):
            FieldCalcObject = DWA.DiWaCAT_FieldCalc()
            FieldCalcObject.ReadMeshValues(self.BeamFile)
            if len(self.InputList) == 14:
                FieldCalcObject.InputPlanarParameters(self.InputList)
            else:
                FieldCalcObject.InputCircularParameters(self.InputList)
            FieldCalcObject.FieldCalculation()
            self.finished.emit(FieldCalcObject)
        else: self.nosim.emit()
    def stop(self):
        self._isRunning = False


def clearLayout(layout):
    if layout is not None:
        while layout.count():
            child = layout.takeAt(0)
            if child.widget() is not None:
                child.widget().deleteLater()
            elif child.layout() is not None:
                clearLayout(child.layout())
    

class CalculationWindow(QWidget):
    '''
    Read in beam and perform field calculation
    '''
    def __init__(self):
        super().__init__()
        self.setWindowTitle("DiWaCAT Field Calculation")
        self.scalelist = ['unit', 'centi', 'milli', 'micro', 'nano', 'pico', 'femto']
        
        windowWidth = 700
        windowHeight = 500
        self.setMinimumSize(windowWidth, windowHeight)
        
        self.layout = pyqt.QGridLayout()
        
        self.LoadButton = pyqt.QPushButton('Load Beam File')
        self.BeamPath = pyqt.QLineEdit()
        self.BeamPath.setPlaceholderText('No File Loaded')
        self.LoadButton.clicked.connect(self.OpenBeamFile)
        self.layout.addWidget(self.LoadButton, 0, 0, 1, 1)
        self.layout.addWidget(self.BeamPath, 0, 1, 1, 2)
        
        self.runButton = pyqt.QPushButton('Run Field Calculation')
        self.runButton.clicked.connect(self.runCalculation)
        self.layout.addWidget(self.runButton, 1, 1, 1, 1)
        
        #List of variables for all simulations
        #Beam Position
        self.CentreBeam = pyqt.QPushButton('Center Beam')
        self.CentreBeam.clicked.connect(self.centerbeam)
        self.x0Input = pyqt.QDoubleSpinBox()
        self.x0Input.setRange(-10000000000000,10000000000000)
        self.x0Scale = pyqt.QComboBox()
        self.x0Scale.addItems(self.scalelist)
        self.y0Input = pyqt.QDoubleSpinBox()
        self.y0Input.setRange(-10000000000000,10000000000000)
        self.y0Scale = pyqt.QComboBox()
        self.y0Scale.addItems(self.scalelist)
        #Max Z Position for Sim
        self.MaxZInput = pyqt.QDoubleSpinBox()
        self.MaxZInput.setRange(0,10000000000000)
        self.MaxZInput.setSingleStep(0.001)
        self.MaxZInput.setValue(0)
        self.MaxZScale = pyqt.QComboBox()
        self.MaxZScale.addItems(['m', 'mm', 'um', 'nm', 's', 'ns', 'ps', 'fs', 'as'])
        self.MaxZConversion = [1, 1e-3, 1e-6, 1e-9, 1/const.c,1e-9 * 1/const.c,1e-12*1/const.c,1e-15*1/const.c,1e-18*1/const.c]
        #Structure
        self.aInput = pyqt.QDoubleSpinBox()
        self.aInput.setRange(0,10000000000000)
        self.aInput.setSingleStep(0.01)
        self.aScale = pyqt.QComboBox()
        self.aScale.addItems(self.scalelist)
        self.deltaInput = pyqt.QDoubleSpinBox()
        self.deltaInput.setRange(0,10000000000000)
        self.deltaInput.setSingleStep(0.01)
        self.deltaScale = pyqt.QComboBox()
        self.deltaScale.addItems(self.scalelist)
        self.PermitivityInput = pyqt.QDoubleSpinBox()
        self.PermitivityInput.setRange(1,100)
        self.PermitivityInput.setValue(3.75)
        self.PermitivityInput.setSingleStep(0.01)
        self.PermeabilityInput = pyqt.QDoubleSpinBox()
        self.PermeabilityInput.setRange(1,100)
        self.PermeabilityInput.setRange(1,100)
        self.PermeabilityInput.setValue(1)
        self.PermeabilityInput.setSingleStep(0.01)
        #Mode Precision and Accuracy
        self.PrecisionInput = pyqt.QDoubleSpinBox()
        self.PrecisionInput.setRange(0,1)
        self.PrecisionInput.setSingleStep(0.0001)
        self.PrecisionInput.setValue(0.01)
        self.AccuracyInput = pyqt.QDoubleSpinBox()
        self.AccuracyInput.setRange(0,100)
        self.AccuracyInput.setSingleStep(0.001)
        self.AccuracyInput.setValue(0.1)
        self.ConvergenceTF = pyqt.QCheckBox()
        self.ConvergenceTF.setCheckState(Qt.Checked)
        
        # Initialize tab screen
        self.tabs = pyqt.QTabWidget()
        self.PlanarTab = QWidget()
        self.CircularTab = QWidget()
        self.PlanarTab.layout = pyqt.QGridLayout(self)
        self.CircularTab.layout = pyqt.QGridLayout(self)
        
        # Add tabs
        self.tabs.addTab(self.PlanarTab,"Planar DLW")
        self.tabs.addTab(self.CircularTab,"Circular DLW")
        
        self.tabs.currentChanged.connect(self.setupTabs)
        self.setupTabs()
        # Add tabs to widget
        self.layout.addWidget(self.tabs, 2, 0, 1, 3)
        self.setLayout(self.layout)
        
    def setupTabs(self):
        if self.tabs.currentIndex() == 0:
            # Create planar tab
            #Drop box for orientation
            self.PlanarOrientation = pyqt.QComboBox()
            self.PlanarOrientation.addItems(['Horizontal', 'Vertical'])
            self.PlanarTab.layout.addWidget(pyqt.QLabel('Planar DLW Orientation'), 0, 0, 1, 1)
            self.PlanarTab.layout.addWidget(self.PlanarOrientation, 0, 1, 1, 2)
            #Add info
            self.PlanarTab.layout.addWidget(pyqt.QLabel('x0 [m]'), 1, 0, 1, 1)
            self.PlanarTab.layout.addWidget(self.x0Input, 1, 1, 1, 1)
            self.PlanarTab.layout.addWidget(self.x0Scale, 1, 2, 1, 1)
            self.PlanarTab.layout.addWidget(pyqt.QLabel('y0 [m]'), 2, 0, 1, 1)
            self.PlanarTab.layout.addWidget(self.y0Input, 2, 1, 1, 1)
            self.PlanarTab.layout.addWidget(self.y0Scale, 2, 2, 1, 1)
            self.PlanarTab.layout.addWidget(pyqt.QLabel('DLW Half-Gap [m]'), 3, 0, 1, 1)
            self.PlanarTab.layout.addWidget(self.aInput, 3, 1, 1, 1)
            self.PlanarTab.layout.addWidget(self.aScale, 3, 2, 1, 1)
            self.PlanarTab.layout.addWidget(pyqt.QLabel('Dielectric Thickness [m]'), 4, 0, 1, 1)
            self.PlanarTab.layout.addWidget(self.deltaInput, 4, 1, 1, 1)
            self.PlanarTab.layout.addWidget(self.deltaScale, 4, 2, 1, 1)
            self.PlanarTab.layout.addWidget(pyqt.QLabel('Dielectric Permitivity'), 5, 0, 1, 1)
            self.PlanarTab.layout.addWidget(self.PermitivityInput, 5, 1, 1, 1)
            self.PlanarTab.layout.addWidget(pyqt.QLabel('Dielectric Permeability'), 6, 0, 1, 1)
            self.PlanarTab.layout.addWidget(self.PermeabilityInput, 6, 1, 1, 1)
            self.WidthInput = pyqt.QDoubleSpinBox()
            self.WidthInput.setRange(0,1000000000)
            self.WidthInput.setValue(1)
            self.WidthScale = pyqt.QComboBox()
            self.WidthScale.addItems(self.scalelist)
            self.PlanarTab.layout.addWidget(pyqt.QLabel('DLW Width [m]'), 7, 0, 1, 1)
            self.PlanarTab.layout.addWidget(self.WidthInput, 7, 1, 1, 1)
            self.PlanarTab.layout.addWidget(self.WidthScale, 7, 2, 1, 1)
            #Simulation Parameters - Include Default Button
            self.DefaultPlanar = pyqt.QPushButton('Default Simulation Parameters')
            self.DefaultPlanar.clicked.connect(self.DefaultParametersP)
            self.PlanarTab.layout.addWidget(self.DefaultPlanar, 8, 0, 1, 2)
            self.PlanarTab.layout.addWidget(pyqt.QLabel('Mode Precision'), 9, 0, 1, 1)
            self.PlanarTab.layout.addWidget(self.PrecisionInput, 9, 1, 1, 1)
            self.PlanarTab.layout.addWidget(pyqt.QLabel('Mode Convergence Accuracy [%]'), 10, 0, 1, 1)
            self.PlanarTab.layout.addWidget(self.AccuracyInput, 10, 1, 1, 1)
            self.nXMin = pyqt.QSpinBox()
            self.nXMin.setRange(5,1000)
            self.nXMin.setValue(15)
            self.PlanarTab.layout.addWidget(pyqt.QLabel('Minimum nX Modes'), 11, 0, 1, 1)
            self.PlanarTab.layout.addWidget(self.nXMin, 11, 1, 1, 1)
            self.nYMin = pyqt.QSpinBox()
            self.nYMin.setRange(5,1000)
            self.nYMin.setValue(15)
            self.PlanarTab.layout.addWidget(pyqt.QLabel('Minimum nY Modes'), 12, 0, 1, 1)
            self.PlanarTab.layout.addWidget(self.nYMin, 12, 1, 1, 1)
            self.PlanarTab.layout.addWidget(pyqt.QLabel('Convergence [True/False]'), 13, 0, 1, 1)
            self.PlanarTab.layout.addWidget(self.ConvergenceTF, 13, 1, 1, 1)
            self.PlanarTab.layout.addWidget(pyqt.QLabel('Maximum Longitudinal Position'), 14, 0, 1, 1)
            self.PlanarTab.layout.addWidget(self.MaxZInput, 14, 1, 1, 1)
            self.PlanarTab.layout.addWidget(self.MaxZScale, 14, 2, 1, 1)
            self.DefaultParametersP()
            self.PlanarTab.setLayout(self.PlanarTab.layout)
        elif self.tabs.currentIndex() == 1:
            #Create circular tab
            #Add info
            self.CircularTab.layout.addWidget(pyqt.QLabel('x0 [m]'), 1, 0, 1, 1)
            self.CircularTab.layout.addWidget(self.x0Input, 1, 1, 1, 1)
            self.CircularTab.layout.addWidget(self.x0Scale, 1, 2, 1, 1)
            self.CircularTab.layout.addWidget(pyqt.QLabel('y0 [m]'), 2, 0, 1, 1)
            self.CircularTab.layout.addWidget(self.y0Input, 2, 1, 1, 1)
            self.CircularTab.layout.addWidget(self.y0Scale, 2, 2, 1, 1)
            self.CircularTab.layout.addWidget(pyqt.QLabel('DLW Half-Gap [m]'), 3, 0, 1, 1)
            self.CircularTab.layout.addWidget(self.aInput, 3, 1, 1, 1)
            self.CircularTab.layout.addWidget(self.aScale, 3, 2, 1, 1)
            self.CircularTab.layout.addWidget(pyqt.QLabel('Dielectric Thickness [m]'), 4, 0, 1, 1)
            self.CircularTab.layout.addWidget(self.deltaInput, 4, 1, 1, 1)
            self.CircularTab.layout.addWidget(self.deltaScale, 4, 2, 1, 1)
            self.CircularTab.layout.addWidget(pyqt.QLabel('Dielectric Permitivity'), 5, 0, 1, 1)
            self.CircularTab.layout.addWidget(self.PermitivityInput, 5, 1, 1, 1)
            self.CircularTab.layout.addWidget(pyqt.QLabel('Dielectric Permeability'), 6, 0, 1, 1)
            self.CircularTab.layout.addWidget(self.PermeabilityInput, 6, 1, 1, 1)
            #Simulation Parameters - Include Default Button
            self.DefaultCircular = pyqt.QPushButton('Default Simulation Parameters')
            self.DefaultCircular.clicked.connect(self.DefaultParametersC)
            self.CircularTab.layout.addWidget(self.DefaultCircular, 7, 0, 1, 2)
            self.CircularTab.layout.addWidget(pyqt.QLabel('Mode Precision'), 8, 0, 1, 1)
            self.CircularTab.layout.addWidget(self.PrecisionInput, 8, 1, 1, 1)
            self.CircularTab.layout.addWidget(pyqt.QLabel('Mode Convergence Accuracy [%]'), 9, 0, 1, 1)
            self.CircularTab.layout.addWidget(self.AccuracyInput, 9, 1, 1, 1)
            self.nRMin = pyqt.QSpinBox()
            self.nRMin.setRange(3,1000)
            self.nRMin.setValue(10)
            self.CircularTab.layout.addWidget(pyqt.QLabel('Minimum Radial Modes'), 10, 0, 1, 1)
            self.CircularTab.layout.addWidget(self.nRMin, 10, 1, 1, 1)
            self.nThetaMin = pyqt.QSpinBox()
            self.nThetaMin.setRange(3,1000)
            self.nThetaMin.setValue(10)
            self.CircularTab.layout.addWidget(pyqt.QLabel('Minimum Azimuthal Modes'), 11, 0, 1, 1)
            self.CircularTab.layout.addWidget(self.nThetaMin, 11, 1, 1, 1)
            self.CircularTab.layout.addWidget(pyqt.QLabel('Convergence [True/False]'), 12, 0, 1, 1)
            self.CircularTab.layout.addWidget(self.ConvergenceTF, 12, 1, 1, 1)
            self.CircularTab.layout.addWidget(pyqt.QLabel('Maximum Longitudinal Position'), 13, 0, 1, 1)
            self.CircularTab.layout.addWidget(self.MaxZInput, 13, 1, 1, 1)
            self.CircularTab.layout.addWidget(self.MaxZScale, 13, 2, 1, 1)
            
            self.DefaultParametersC()
            self.CircularTab.setLayout(self.CircularTab.layout)
        
    def centerbeam(self):   
        self.x0Input.setValue(0)
        self.y0Input.setValue(0)
         
    def DefaultParametersP(self):
        self.PrecisionInput.setValue(0.01)
        self.AccuracyInput.setValue(0.1)
        self.nXMin.setValue(15)
        self.nYMin.setValue(15)
        self.ConvergenceTF.setCheckState(Qt.Checked)
    def DefaultParametersC(self):
        self.PrecisionInput.setValue(0.1)
        self.AccuracyInput.setValue(0.1)
        self.nRMin.setValue(10)
        self.nThetaMin.setValue(10)
        self.ConvergenceTF.setCheckState(Qt.Checked)
        
    def OpenBeamFile(self):
        filename, ok = pyqt.QFileDialog.getOpenFileName(self,"Select a File", "", "HDF5 Files (*.h5)")
        if filename:
            self.BeamPath.setText(filename)
            
    def runCalculation(self):
        Inputs = []
        convergence = 'f'
        if self.ConvergenceTF.isChecked(): convergence = 't'
        #Lock all variables
        if self.tabs.currentIndex() == 0:
            #First set all 
            for i in reversed(range(self.PlanarTab.layout.count())): 
                self.PlanarTab.layout.itemAt(i).widget().setEnabled(False)
            #Create the input list
            orientation = ['h','v'][self.PlanarOrientation.currentIndex()]
            Inputs = [orientation, self.x0Input.value() * 1/UnitConversion(self.x0Scale.currentText()), self.y0Input.value() * 1/UnitConversion(self.y0Scale.currentText()),
                      self.PermitivityInput.value(), self.PermeabilityInput.value(), self.aInput.value() * 1/UnitConversion(self.aScale.currentText()),
                      self.deltaInput.value() * 1/UnitConversion(self.deltaScale.currentText()), self.WidthInput.value() * 1/UnitConversion(self.WidthScale.currentText()),
                      self.nXMin.value(), self.nYMin.value(), self.PrecisionInput.value(), self.AccuracyInput.value(),
                      self.MaxZInput.value() * self.MaxZConversion[self.MaxZScale.currentIndex()], convergence]
        elif self.tabs.currentIndex() == 1:
            for i in reversed(range(self.CircularTab.layout.count())): 
                self.CircularTab.layout.itemAt(i).widget().setEnabled(False)
            Inputs = [self.x0Input.value() * 1/UnitConversion(self.x0Scale.currentText()), self.y0Input.value() * 1/UnitConversion(self.y0Scale.currentText()),
                      self.PermitivityInput.value(), self.PermeabilityInput.value(), self.aInput.value() * 1/UnitConversion(self.aScale.currentText()),
                      self.deltaInput.value() * 1/UnitConversion(self.deltaScale.currentText()),
                      self.nRMin.value(), self.nThetaMin.value(), self.PrecisionInput.value(), self.AccuracyInput.value(),
                      self.MaxZInput.value() * self.MaxZConversion[self.MaxZScale.currentIndex()], convergence]
        self.tabs.setEnabled(False)
        self.BeamPath.setEnabled(False)
        self.LoadButton.setEnabled(False)
        self.runButton.setText('Simulation Running')
        self.runButton.setEnabled(False)
        #Send calculation to an external thread to stop wheel of death - same as with tracking
        self.runFunction = CalculationWorkerFunction(Inputs, self.BeamPath.text())
        self.runFunction.finished.connect(self.SimulationFinished)
        self.runFunction.nosim.connect(self.NoSimRun)
        self.runFunction.start()
        
        
    def SimulationFinished(self, FieldOutput):
        #Unlock all variables
        self.stop_thread()
        if self.tabs.currentIndex() == 0:
            for i in reversed(range(self.PlanarTab.layout.count())): 
                self.PlanarTab.layout.itemAt(i).widget().setEnabled(True)
        elif self.tabs.currentIndex() == 1:
            for i in reversed(range(self.CircularTab.layout.count())): 
                self.CircularTab.layout.itemAt(i).widget().setEnabled(True)
        self.tabs.setEnabled(True)
        self.BeamPath.setEnabled(True)
        self.LoadButton.setEnabled(True)
        self.runButton.setText('Run Field Calculation')
        self.runButton.setEnabled(True)
        filename, ok = pyqt.QFileDialog.getSaveFileName(self,"Save File Name", "", "HDF5 Files (*.h5)")
        if filename:
            path = Path(filename)
            FieldOutput.write_Field_HDF5(path)
        else:
            print('Field Not Saved')
    def NoSimRun(self):
        print('File not Found')
        self.stop_thread()
        #Unlock variables
        if self.tabs.currentIndex() == 0:
            for i in reversed(range(self.PlanarTab.layout.count())): 
                self.PlanarTab.layout.itemAt(i).widget().setEnabled(True)
        elif self.tabs.currentIndex() == 1:
            for i in reversed(range(self.CircularTab.layout.count())): 
                self.CircularTab.layout.itemAt(i).widget().setEnabled(True)
        self.tabs.setEnabled(True)
        self.BeamPath.setEnabled(True)
        self.LoadButton.setEnabled(True)
        self.runButton.setText('Run Field Calculation')
        self.runButton.setEnabled(True)
        
    def stop_thread(self):
        #Make sure our thread is closed - attempt to stop memory leaks
        self.runFunction.stop()
        self.runFunction.quit()
        self.runFunction.wait()
        
    