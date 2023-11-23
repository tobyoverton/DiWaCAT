# -*- coding: utf-8 -*-
"""
Beam Tools GUI
18/10/2023
T Overton

Bringing together our beam manipulation tools into a GUI
Non-DWA functions so kept together
Using the functions from beam_tools and beam_to_screen
3 operations: slicing, transport, screen simulation
"""

from PyQt5.QtWidgets import QWidget
from PyQt5.QtGui import QFont
import PyQt5.QtWidgets as pyqt
import scipy.constants as const
import numpy as np
import pyqtgraph as pg
from pathlib import Path

from Python_Tools.Modules import beam_tools as dbt
from Python_Tools.Modules import beam_to_screen as bts

def UnitConversion(prefix):
    if prefix == 'unit' or prefix == 'Unit':
        return 1
    else:
        return getattr(const, prefix)

#Main set of windows - i.e. selecting the function we want. Panel in the centre will then be the selected function
class BeamToolsWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Beam Manipulation and Tools")
        self.setMinimumSize(1000,400)
        
        self.layout = pyqt.QGridLayout()
        
        #Constants
        self.ScaleList = ['unit', 'centi', 'milli', 'micro', 'nano', 'pico', 'femto', 'kilo', 'mega', 'giga']
        self.scaleunit = ['', 'c', 'm', 'u', 'n', 'p', 'f', 'k', 'M', 'G']
        self.ParameterList = ['x','y','z','t','r','px','py','pz','cpx','cpy','cpz','xp','yp']
        self.ParameterListSmaller = ['x','y','z','t','r','px','py','pz']
        
        #Input Beam File
        self.LoadButton = pyqt.QPushButton('Load Beam File')
        self.filename = ''
        self.BeamPath = pyqt.QLineEdit()
        self.BeamPath.setPlaceholderText('No File Loaded')
        self.LoadButton.clicked.connect(self.OpenBeamFile)
        self.layout.addWidget(self.LoadButton, 0, 0, 1, 1)
        self.layout.addWidget(self.BeamPath, 0, 1, 1, 2)
        self.beamIn = dbt.BeamFromFile()
        self.beamList = [];
        self.beamListLabels = [];
        self.CurrentBeam = pyqt.QLabel('')
        self.BeamListWidget = pyqt.QListWidget()
        
        
        # Initialize tab screen
        self.tabs = pyqt.QTabWidget()
        self.SliceTab = QWidget()
        self.SliceTabSetup()
        self.SliceTab.setLayout(self.SliceTab.layout)
        self.ManipulationTab = QWidget()
        self.ManipulateTabSetup()
        self.ManipulationTab.setLayout(self.ManipulationTab.layout)
        self.ScreenSimTab = QWidget()
        self.ScreenTabSetup()
        self.ScreenSimTab.setLayout(self.ScreenSimTab.layout)
        
        # Add tabs
        self.tabs.addTab(self.SliceTab,"Beam Slices")
        self.tabs.addTab(self.ManipulationTab,"Beam Manipulation/Transport")
        self.tabs.addTab(self.ScreenSimTab,"Screen Simulation")
        
        # Add tabs to widget
        self.layout.addWidget(self.tabs, 1, 0, 1, 3)
        
        self.setLayout(self.layout)
        
    ### --------------------------------------------------------------------------- ###    
    ### ---------------------- Slice Function Definitions ------------------------- ###
    ### --------------------------------------------------------------------------- ###
    def SliceTabSetup(self):
        self.SliceTab.layout = pyqt.QGridLayout()
        self.SliceParameter = pyqt.QComboBox()
        self.SliceParameter.addItems(self.ParameterListSmaller)
        self.SliceTab.layout.addWidget(pyqt.QLabel('Slice Parameter: '), 0, 0, 1, 1)
        self.SliceTab.layout.addWidget(self.SliceParameter, 0, 1, 1, 1)
        
        self.SliceScale = pyqt.QComboBox()
        self.SliceScale.addItems(self.ScaleList)
        self.SliceTab.layout.addWidget(pyqt.QLabel('Slice Parameter Scale: '),1,0,1,1)
        self.SliceTab.layout.addWidget(self.SliceScale,1,1,1,1)
        
        self.SliceMin = pyqt.QDoubleSpinBox()
        self.SliceMin.setRange(-999999,999999)
        self.SliceMax = pyqt.QDoubleSpinBox()
        self.SliceMax.setRange(-999999,999999)
        self.SliceTab.layout.addWidget(pyqt.QLabel('Slice Minimum Value'), 2, 0, 1, 1)
        self.SliceTab.layout.addWidget(self.SliceMin, 2, 1, 1, 1)
        self.SliceTab.layout.addWidget(pyqt.QLabel('Slice Maximum Value'), 3, 0, 1, 1)
        self.SliceTab.layout.addWidget(self.SliceMax, 3, 1, 1, 1)
        
        self.NumberOfSlices = pyqt.QSpinBox()
        self.NumberOfSlices.setRange(2,100)
        self.NumberOfSlices.setValue(5)
        self.SliceTab.layout.addWidget(pyqt.QLabel('Number of Slices'), 4, 0, 1, 1)
        self.SliceTab.layout.addWidget(self.NumberOfSlices, 4, 1, 1, 1)
        
        self.SliceTestButton = pyqt.QPushButton('Plot Slices')
        self.SliceTestButton.clicked.connect(self.SliceTest)
        self.SliceSaveButton = pyqt.QPushButton('Slicing and Save to Folder')
        self.SliceSaveButton.clicked.connect(self.SaveSlice)
        self.SliceTab.layout.addWidget(self.SliceTestButton, 5, 0, 1, 2)
        self.SliceTab.layout.addWidget(self.SliceSaveButton, 6, 0, 1, 2)
        
        self.SliceHistogram = pg.PlotWidget()
        self.SliceHistogram.setBackground('w')
        self.SliceHistogram.getAxis('bottom').setPen('k')
        self.SliceHistogram.getAxis('left').setPen('k')
        self.SliceHistogram.getAxis('bottom').setTextPen('k')
        self.SliceHistogram.getAxis('left').setTextPen('k')
        self.SliceHistogram.setLabel('left','No. of Macroparticles')
        self.HistValues = [[[0],[]]]
        self.Histogram = pg.PlotCurveItem( self.HistValues[0][0],self.HistValues[0][1], stepMode=True, fillLevel=0)
        self.SliceHistogram.addItem(self.Histogram)
        
        
        self.SliceTab.layout.addWidget(self.SliceHistogram,0,3,6,3)
    
    def SaveSlice(self):
        self.PerformSlicing(True)
    def SliceTest(self):
        self.PerformSlicing(False)
    def PerformSlicing(self, SaveTrueFalse):
        if self.filename != '':
            slicer = dbt.BeamSlicer(self.beamIn, self.SliceParameter.currentText())
            NSlices = self.NumberOfSlices.value();
            Scaling = UnitConversion(self.SliceScale.currentText())
            slice_w = Scaling * (self.SliceMax.value() - self.SliceMin.value()) * 1/NSlices
            minimum = (self.SliceMin.value() * Scaling) + 0.5*slice_w
            maximum = (self.SliceMax.value() * Scaling) - 0.5*slice_w
            slicer.set_slice_parameters(min_p = minimum,max_p = maximum, N_slices=NSlices, slice_width=slice_w)
            slicer.perform_slicing()
            
            slicepositions = (slicer.slice_centres - 0.5*slice_w) * 1/Scaling
            slicepositions = np.append(slicepositions,slicepositions[-1] + slice_w*1/Scaling)
            slicemacros = [beam.nMacros for beam in slicer.sliced_beams]
            
            self.HistValues[0][0] = slicepositions
            self.HistValues[0][1] = slicemacros
            Histogram = pg.PlotCurveItem( self.HistValues[0][0],self.HistValues[0][1], stepMode=True, fillLevel=0, pen = 'k')
            self.SliceHistogram.clear()
            self.SliceHistogram.addItem(Histogram)
            self.SliceHistogram.setLabel('bottom',self.SliceParameter.currentText() +' [' + self.scaleunit[self.SliceScale.currentIndex()] + self.beamIn.parameter_units[self.SliceParameter.currentText()] + ']')
            
            if SaveTrueFalse == True:
                print('True')
                #Open dialog to find folder
                folderpath = pyqt.QFileDialog.getExistingDirectory(self, "Select a Directory")
                if folderpath:
                #If the folder is found save beams
                    for i in range(len(slicer.sliced_beams)):
                        beamSelected = slicer.sliced_beams[i]
                        path = Path(folderpath + '/Slice' + str(i+1) + '.h5')
                        beamSelected.write_HDF5_beam_file(path,overwrite_existing_file=True)
                    print('Files Written')
                else:
                    print('Sliced Beams not Saved')
            else:
                print('Sliced Beams not Saved')
        else:
            print('No Beam Loaded')
    
    ### --------------------------------------------------------------------------- ###
    ### ------------------------ Beam Manipulation Functions ---------------------- ###
    ### --------------------------------------------------------------------------- ###
    def ManipulateTabSetup(self):
        self.ManipulationTab.layout = pyqt.QGridLayout()
        
        self.OpticsFunctions = ["Thin-Lens Quadrupole", "Thick-Lens Quadrupole", "Dipole", "Drift Space", "Collimate"]
        self.OptionSelected = pyqt.QComboBox()
        self.OptionSelected.addItems(self.OpticsFunctions)
        
        self.ManipulationTab.layout.addWidget(pyqt.QLabel('Current Beam: '), 0, 0, 1, 1)
        self.ManipulationTab.layout.addWidget(self.CurrentBeam, 0, 1, 1, 1)
        
        # Initialize tab screen
        self.ManipulationFunctionTabs = pyqt.QTabWidget()
        self.tablist = [pyqt.QWidget() for _ in range(len(self.OpticsFunctions))]
        self.ManipulationFunctionSetup()
        # Add tabs
        for i, tab in enumerate(self.tablist):
            tab.setLayout(tab.layout)
            self.ManipulationFunctionTabs.addTab(tab, self.OpticsFunctions[i])
        # Add tabs to widget
        self.ManipulationFunctionTabs.setFixedHeight(200)
        self.ManipulationTab.layout.addWidget(self.ManipulationFunctionTabs, 1, 0, 1, 2)
        
        
        #Add List of Beams
        self.BeamListWidget.setFixedHeight(250)
        self.ManipulationTab.layout.addWidget(self.BeamListWidget, 0, 2, 2, 2)
        #Reset Buttons
        self.ResetToSelected = pyqt.QPushButton('Reset to Selected Beam')
        self.ResetToSelected.clicked.connect(self.ResetSelect)
        self.ResetToBeamIn = pyqt.QPushButton('Reset to Beam In')
        self.ResetToBeamIn.clicked.connect(self.ResetAll)
        self.ManipulationTab.layout.addWidget(self.ResetToSelected, 2,2, 1, 1)
        self.ManipulationTab.layout.addWidget(self.ResetToBeamIn, 2,3, 1, 1)
        #Save Button
        self.SaveSelectButton = pyqt.QPushButton('Save Selected Beam')
        self.SaveSelectButton.clicked.connect(self.SaveSelect)
        self.SaveAllButton = pyqt.QPushButton('Save All Beams')
        self.SaveAllButton.clicked.connect(self.SaveAll)
        self.ManipulationTab.layout.addWidget(self.SaveSelectButton, 3,2, 1, 1)
        self.ManipulationTab.layout.addWidget(self.SaveAllButton, 3,3, 1, 1)
        
        self.ApplyFunctionButton = pyqt.QPushButton('Apply Function')
        self.ApplyFunctionButton.clicked.connect(self.ApplyFunction)
        self.ManipulationTab.layout.addWidget(self.ApplyFunctionButton, 2, 0, 1, 2)
    def ManipulationFunctionSetup(self):
        #First tab is the thin lens
        self.tablist[0].layout = pyqt.QGridLayout()
        self.FocalLength = pyqt.QDoubleSpinBox()
        self.FocalLength.setRange(-999999, 999999)
        self.FocalLengthScale = pyqt.QComboBox()
        self.FocalLengthScale.addItems(self.ScaleList)
        self.tablist[0].layout.addWidget(pyqt.QLabel('Focal Length (Horizontal) [m]:'), 0, 0, 1, 1)
        self.tablist[0].layout.addWidget(self.FocalLength, 0, 1, 1, 1)
        self.tablist[0].layout.addWidget(self.FocalLengthScale, 0, 2, 1, 1)
        
        #Second tab is the thick lens
        self.tablist[1].layout = pyqt.QGridLayout()
        self.ThickQuadStrength = pyqt.QDoubleSpinBox()
        self.ThickQuadStrength.setRange(-999999,999999)
        self.ThickQuadStrengthScale = pyqt.QComboBox()
        self.ThickQuadStrengthScale.addItems(self.ScaleList)
        self.tablist[1].layout.addWidget(pyqt.QLabel('Quadrupole Strength [m^-2]'), 0, 0, 1, 1)
        self.tablist[1].layout.addWidget(self.ThickQuadStrength, 0, 1, 1, 1)
        self.tablist[1].layout.addWidget(self.ThickQuadStrengthScale, 0, 2, 1, 1)
        self.ThickQuadLength= pyqt.QDoubleSpinBox()
        self.ThickQuadLength.setRange(-999999,999999)
        self.ThickQuadLengthScale = pyqt.QComboBox()
        self.ThickQuadLengthScale.addItems(self.ScaleList)
        self.tablist[1].layout.addWidget(pyqt.QLabel('Quadrupole Length [m]'), 1, 0, 1, 1)
        self.tablist[1].layout.addWidget(self.ThickQuadLength, 1, 1, 1, 1)
        self.tablist[1].layout.addWidget(self.ThickQuadLengthScale, 1, 2, 1, 1)
        
        #Third tab is dipole
        self.tablist[2].layout = pyqt.QGridLayout()
        self.tablist[2].layout.addWidget(pyqt.QLabel('Dipole Function Needs Writing'))
        
        #Fourth tab is only drift
        self.tablist[3].layout = pyqt.QGridLayout()
        self.DriftDistance = pyqt.QDoubleSpinBox()
        self.DriftDistance.setRange(-999999,999999)
        self.DriftScale = pyqt.QComboBox()
        self.DriftScale.addItems(self.ScaleList)
        self.tablist[3].layout.addWidget(pyqt.QLabel('Drift Distance [m]'), 0, 0, 1, 1)
        self.tablist[3].layout.addWidget(self.DriftDistance, 0, 1, 1, 1)
        self.tablist[3].layout.addWidget(self.DriftScale, 0, 2, 1, 1)
        
        
        #Fifth tab is collimation
        self.tablist[4].layout = pyqt.QGridLayout()
        self.CollimationVariable = pyqt.QComboBox()
        self.CollimationVariable.addItems(self.ParameterListSmaller)
        self.tablist[4].layout.addWidget(pyqt.QLabel('Collimation Variable'), 0, 0, 1, 1)
        self.tablist[4].layout.addWidget(self.CollimationVariable, 0, 1, 1, 1)
        self.CollimationScale = pyqt.QComboBox()
        self.CollimationScale.addItems(self.ScaleList)
        self.tablist[4].layout.addWidget(pyqt.QLabel('Collimation Limit Scale'), 1, 0, 1, 1)
        self.tablist[4].layout.addWidget(self.CollimationScale, 1, 1, 1, 1)
        self.CollimationMin = pyqt.QDoubleSpinBox()
        self.CollimationMin.setRange(-999999,999999)
        self.CollimationMax = pyqt.QDoubleSpinBox()
        self.CollimationMax.setRange(-999999,999999)
        self.tablist[4].layout.addWidget(pyqt.QLabel('Collimation Minimum Limit'), 2, 0, 1, 1)
        self.tablist[4].layout.addWidget(self.CollimationMin, 2, 1, 1, 1)
        self.tablist[4].layout.addWidget(pyqt.QLabel('Collimation Maximum Limit'), 3, 0, 1, 1)
        self.tablist[4].layout.addWidget(self.CollimationMax, 3, 1, 1, 1)
        
        
        
    def ResetSelect(self):
        IndexToResetTo = self.BeamListWidget.currentRow()
        self.beamList = self.beamList[:IndexToResetTo+1]
        self.beamListLabels = self.beamListLabels[:IndexToResetTo+1]
        self.CurrentBeam.setText(self.beamListLabels[-1])
        self.BeamListWidget.clear()
        self.BeamListWidget.addItems(self.beamListLabels)
    def ResetAll(self):
        if self.filename != '':
            self.beamList = [self.beamIn]
            self.beamListLabels = ["BeamIn"]
            self.CurrentBeam.setText("BeamIn")
            self.BeamListWidget.clear()
            self.BeamListWidget.addItem("BeamIn")
    def SaveSelect(self):
        beamSelected = self.beamList[self.BeamListWidget.currentRow()]
        filename, ok = pyqt.QFileDialog.getSaveFileName(self,"Save File Name", "", "HDF5 Files (*.h5)")
        if filename:
            path = Path(filename)
            beamSelected.write_HDF5_beam_file(path,overwrite_existing_file=True)
        else:
            print('Beam Not Saved')
    def SaveAll(self):
        folderpath = pyqt.QFileDialog.getExistingDirectory(self, "Select a Directory")
        if folderpath:
        #If the folder is found save beams
            for i in range(len(self.beamList)):
                beamSelected = self.beamList[i]
                path = Path(folderpath + '//' + self.beamListLabels[i] + '.h5')
                beamSelected.write_HDF5_beam_file(path,overwrite_existing_file=True)
            print('Files Written')
        else:
            print('Sliced Beams not Saved')
    def ApplyFunction(self):
        if len(self.beamList)>0:
            beamnew = self.beamList[-1].returnBeam()
            beamlabelpreppend = 'Beam' + str(len(self.beamList) + 1) + '_'
            if self.ManipulationFunctionTabs.currentIndex() == 0:
                beamnew.thinLensFQuadBeam(self.FocalLength.value() * UnitConversion(self.FocalLengthScale.currentText()))
                self.beamList = np.append(self.beamList,beamnew)
                self.beamListLabels = np.append(self.beamListLabels, beamlabelpreppend+'thinQuad_f_'+str(self.FocalLength.value())+'_'+self.FocalLengthScale.currentText())
            if self.ManipulationFunctionTabs.currentIndex() == 1:
                beamnew.thickQuadBeam(self.ThickQuadStrength.value() * UnitConversion(self.ThickQuadStrengthScale.currentText()),self.ThickQuadLength.value() * UnitConversion(self.ThickQuadLengthScale.currentText()))
                beamnew._beam['z'] = beamnew._beam['z'] - self.ThickQuadLength.value() * UnitConversion(self.ThickQuadLengthScale.currentText())
                self.beamList = np.append(self.beamList,beamnew)
                self.beamListLabels = np.append(self.beamListLabels,beamlabelpreppend+ 'thickQuad_K_'+str(self.ThickQuadStrength.value())+'_'+self.ThickQuadStrengthScale.currentText()+'L_'+str(self.ThickQuadLength.value())+'_'+self.ThickQuadLengthScale.currentText())
            if self.ManipulationFunctionTabs.currentIndex() == 2:
                print('Dipole Function not Written')
                self.beamList = np.append(self.beamList,beamnew)
                self.beamListLabels = np.append(self.beamListLabels,beamlabelpreppend+ 'Drift_0')
            if self.ManipulationFunctionTabs.currentIndex() == 3:
                beamnew.driftBeam(self.DriftDistance.value() * UnitConversion(self.DriftScale.currentText()))
                beamnew._beam['z'] = beamnew._beam['z'] - self.DriftDistance.value() * UnitConversion(self.DriftScale.currentText())
                self.beamList = np.append(self.beamList,beamnew)
                self.beamListLabels = np.append(self.beamListLabels, beamlabelpreppend+'Drift_'+str(self.DriftDistance.value())+'_'+self.DriftScale.currentText())
            if self.ManipulationFunctionTabs.currentIndex() == 4:
                beamnew.collimateBeam(self.CollimationVariable.currentText(), self.CollimationMin.value() * UnitConversion(self.CollimationScale.currentText()), self.CollimationMax.value() * UnitConversion(self.CollimationScale.currentText()))
                self.beamList = np.append(self.beamList,beamnew)
                self.beamListLabels = np.append(self.beamListLabels,beamlabelpreppend+ 'Collimation_'+self.CollimationVariable.currentText() + '_'+self.CollimationScale.currentText() + '_' + str(self.CollimationMin.value())+ '_' + str(self.CollimationMax.value()) + self.CollimationScale.currentText())

            self.BeamListWidget.addItem(self.beamListLabels[-1])
            self.CurrentBeam.setText(self.beamListLabels[-1])
            print('Function Applied')
        else:
            print('No Beam Loaded')
        
        

    ### --------------------------------------------------------------------------- ###
    ### ------------------------- Screen Simulation ------------------------------- ###
    ### --------------------------------------------------------------------------- ###
    def ScreenTabSetup(self):
        self.ScreenSimTab.layout = pyqt.QGridLayout()
        
        self.ScreenScale = pyqt.QComboBox()
        self.ScreenScale.addItems(self.ScaleList)
        
        
        self.ScreenSimTab.layout.addWidget(pyqt.QLabel('Screen Unit Scale: '), 0, 0, 1, 1)
        self.ScreenSimTab.layout.addWidget(self.ScreenScale, 0, 1, 1, 1)
        
        self.PixSizeX = pyqt.QDoubleSpinBox()
        self.PixSizeX.setRange(0.01,1000000)
        self.PixSizeX.setValue(20)
        self.PixSizeY = pyqt.QDoubleSpinBox()
        self.PixSizeY.setRange(0.01,1000000)
        self.PixSizeY.setValue(20)
        self.ScreenSimTab.layout.addWidget(pyqt.QLabel('Horizontal Pixel Size [Scale]: '), 1, 0, 1, 1)
        self.ScreenSimTab.layout.addWidget(self.PixSizeX, 1, 1, 1, 1)
        self.ScreenSimTab.layout.addWidget(pyqt.QLabel('Vertical Pixel Size [Scale]: '), 2, 0, 1, 1)
        self.ScreenSimTab.layout.addWidget(self.PixSizeY, 2, 1, 1, 1)
        
        self.nPixX = pyqt.QSpinBox()
        self.nPixX.setRange(10,1000000)
        self.nPixX.setValue(2000)
        self.CentralPixX = pyqt.QSpinBox()
        self.CentralPixX.setRange(10,1000000)
        self.CentralPixX.setValue(1000)
        self.nPixY = pyqt.QSpinBox()
        self.nPixY.setRange(10,1000000)
        self.nPixY.setValue(2000)
        self.CentralPixY = pyqt.QSpinBox()
        self.CentralPixY.setRange(10,1000000)
        self.CentralPixY.setValue(1000)
        
        self.ScreenSimTab.layout.addWidget(pyqt.QLabel('Number of Horizontal Pixels: '), 3, 0, 1, 1)
        self.ScreenSimTab.layout.addWidget(self.nPixX, 3, 1, 1, 1)
        self.ScreenSimTab.layout.addWidget(pyqt.QLabel('Central Horizontal Pixel: '), 4, 0, 1, 1)
        self.ScreenSimTab.layout.addWidget(self.CentralPixX, 4, 1, 1, 1)
        self.ScreenSimTab.layout.addWidget(pyqt.QLabel('Number of Vertical Pixels: '), 5, 0, 1, 1)
        self.ScreenSimTab.layout.addWidget(self.nPixY, 5, 1, 1, 1)
        self.ScreenSimTab.layout.addWidget(pyqt.QLabel('Central Vertical Pixel: '), 6, 0, 1, 1)
        self.ScreenSimTab.layout.addWidget(self.CentralPixY, 6, 1, 1, 1)
        
        self.PhysicalResolution = pyqt.QDoubleSpinBox()
        self.PhysicalResolution.setRange(0.01, 100000)
        self.PhysicalResolution.setValue(50)
        self.ScreenSimTab.layout.addWidget(pyqt.QLabel('Physical Resolution [Scale]: '), 7, 0, 1, 1)
        self.ScreenSimTab.layout.addWidget(self.PhysicalResolution, 7, 1, 1, 1)
        
        
        self.ScreenDriftDistance = pyqt.QDoubleSpinBox()
        self.ScreenDriftDistance.setRange(-999999,999999)
        self.ScreenDriftScale = pyqt.QComboBox()
        self.ScreenDriftScale.addItems(self.ScaleList)
        self.ScreenSimTab.layout.addWidget(pyqt.QLabel('Drift Before Screen [m]: '), 8, 0, 1, 1)
        self.ScreenSimTab.layout.addWidget(self.ScreenDriftDistance, 8, 1, 1, 1)
        self.ScreenSimTab.layout.addWidget(self.ScreenDriftScale, 8, 2, 1, 1)
        
        self.ScreenColourScheme = pyqt.QComboBox()
        self.ScreenColourScheme.addItems(['viridis', 'jet', 'turbo','plasma', 'rainbow','binary_r','Greys', 'Reds','Oranges'])
        self.ScreenSimTab.layout.addWidget(pyqt.QLabel('Screen Colour Scheme: '), 9, 0, 1, 1)
        self.ScreenSimTab.layout.addWidget(self.ScreenColourScheme, 9, 1, 1, 1)
        
        
        self.ScreenSimButton = pyqt.QPushButton('Simulate Screen Image')
        self.ScreenSimButton.clicked.connect(self.ScreenSim)
        self.ScreenSimTab.layout.addWidget(self.ScreenSimButton, 12, 0, 1, 2)
        
        self.Font = QFont("Times", 10)
        self.FontSelect = pyqt.QComboBox()
        self.FontOptions = ["Times","Arial", "Helvetica", "Serif"]
        self.FontSelect.addItems(self.FontOptions)
        self.FontSize = pyqt.QSpinBox()
        self.FontSize.setRange(1, 30)
        self.FontSize.setValue(10)
        self.FontSelect.currentIndexChanged.connect(self.ScreenFontChange)
        self.FontSize.valueChanged.connect(self.ScreenFontChange)
        self.ScreenSimTab.layout.addWidget(pyqt.QLabel('Font:'), 10, 0, 1, 1)
        self.ScreenSimTab.layout.addWidget(self.FontSelect, 10, 1, 1, 1)
        self.ScreenSimTab.layout.addWidget(pyqt.QLabel('Font Size:'), 11, 0, 1, 1)
        self.ScreenSimTab.layout.addWidget(self.FontSize, 11, 1, 1, 1)
        
        
        self.ScreenImage = pg.PlotWidget()
        self.ScreenImage.setMinimumSize(300,300)
        self.ScreenImage.setBackground('w')
        self.ScreenImage.showAxis('top')
        self.ScreenImage.showAxis('right')
        self.ScreenImage.getAxis('bottom').setPen('k')
        self.ScreenImage.getAxis('left').setPen('k')
        self.ScreenImage.getAxis('top').setPen('k')
        self.ScreenImage.getAxis('right').setPen('k')
        self.ScreenImage.getAxis('bottom').setTextPen('k')
        self.ScreenImage.getAxis('left').setTextPen('k')
        self.ScreenImage.getAxis('top').setTextPen('w')
        self.ScreenImage.getAxis('right').setTextPen('w')
        self.ScreenImage.getAxis("bottom").setStyle(tickFont = self.Font)
        self.ScreenImage.getAxis("left").setStyle(tickFont = self.Font)
        self.ScreenImage.getAxis("bottom").label.setFont(self.Font)
        self.ScreenImage.getAxis("left").label.setFont(self.Font)
        self.ScreenImage.getViewBox().setBackgroundColor(pg.colormap.getFromMatplotlib(self.ScreenColourScheme.currentText()).map(0.0))
        self.ScreenSimTab.layout.addWidget(self.ScreenImage, 0, 2, 13, 5)

    def ScreenSim(self):
        print('Simulating Screen')
        
        Scaling = UnitConversion(self.ScreenScale.currentText())
        screenDict = {}
        screenDict["n_pix_x"] = self.nPixX.value()
        screenDict["n_pix_y"] = self.nPixY.value()
        screenDict["pix_size_x"] = self.PixSizeX.value() * Scaling
        screenDict["pix_size_y"] = self.PixSizeY.value() * Scaling
        screenDict["central_pix_x"] = self.CentralPixX.value()
        screenDict["central_pix_y"] = self.CentralPixY.value()
        screenDict["base_constant"] = 0
        screenDict["beam_pix_pad"] = 5
        screenDict["physical_res"] = self.PhysicalResolution.value() * Scaling
        screenDict["KDE_method"] = "fastKDE"
        
        self.beamIn.driftBeam(self.ScreenDriftDistance.value() * UnitConversion(self.ScreenDriftScale.currentText()))
        
        screen = bts.BeamToScreen(screenDict, self.beamIn)
        
        self.beamIn.driftBeam(-1*self.ScreenDriftDistance.value() * UnitConversion(self.ScreenDriftScale.currentText()))
        
        #Do the digitisation
        bit_depth = 12
        peak_pix = 2**bit_depth
        screen.set_digitise_pix_by_max_value(bit_depth=bit_depth, max_pix_val=peak_pix)
        
        #Screen arrays are the pixel centres - we need the pixel edges
        xaxis = screen.beam_x_vals - 0.5*self.PixSizeX.value() * Scaling
        xaxis=np.append(xaxis, xaxis[-1]+self.PixSizeX.value() * Scaling)
        yaxis = screen.beam_y_vals - 0.5*self.PixSizeY.value() * Scaling
        yaxis=np.append(yaxis, yaxis[-1]+self.PixSizeY.value() * Scaling)
        X,Y = np.meshgrid(xaxis,yaxis)
        
        colourmap = pg.colormap.getFromMatplotlib(self.ScreenColourScheme.currentText())
        zerocolour = colourmap.map(0.0)
        
        self.ScreenImage.getViewBox().setBackgroundColor(zerocolour)
        
        self.ScreenImage.setLabel('left','y ['+ self.scaleunit[self.ScreenScale.currentIndex()] +self.scaleunit[self.SliceScale.currentIndex()] + 'm]')
        self.ScreenImage.setLabel('bottom','x ['+ self.scaleunit[self.ScreenScale.currentIndex()] +self.scaleunit[self.SliceScale.currentIndex()] + 'm]')
        self.Image = pg.PColorMeshItem(X * 1/Scaling,Y * 1/Scaling, screen.beam_processed_pix, colorMap = colourmap)
        self.ScreenImage.clear()
        self.ScreenImage.addItem(self.Image)
    
        
    def ScreenFontChange(self):
        self.Font = QFont(self.FontSelect.currentText(), self.FontSize.value())
        self.ScreenImage.getAxis("bottom").setStyle(tickFont = self.Font)
        self.ScreenImage.getAxis("left").setStyle(tickFont = self.Font)
        self.ScreenImage.getAxis("bottom").label.setFont(self.Font)
        self.ScreenImage.getAxis("left").label.setFont(self.Font)
        
    def OpenBeamFile(self):
        filename, ok = pyqt.QFileDialog.getOpenFileName(self,"Select a File", "", "HDF5 Files (*.h5)")
        if filename:
            self.filename = filename
            self.BeamPath.setText(self.filename)
            self.beamIn.read_WakeCode_beam_file(self.filename)        
            self.beamList = [self.beamIn];        
            self.beamListLabels = ["BeamIn"];
            self.CurrentBeam.setText("BeamIn")
            self.BeamListWidget.clear()
            self.BeamListWidget.addItem("BeamIn")
        else:
            print('No Beam Loaded')