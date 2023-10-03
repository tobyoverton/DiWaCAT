# -*- coding: utf-8 -*-
"""
Field Plotter GUI Window
T Overton
3/10/2023

Read in field file from a DiWaCAT simulation and plot
"""

from PyQt5.QtWidgets import QWidget
import PyQt5.QtWidgets as pyqt
import pyqtgraph as pg
import numpy as np
import scipy.constants as const

from Python_Tools.Modules import diwacat_tools as DWA

import matplotlib.pyplot as plt
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']


def UnitConversion(prefix):
    if prefix == 'unit' or prefix == 'Unit':
        return 1
    else:
        return 1/getattr(const, prefix)

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
