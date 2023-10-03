# -*- coding: utf-8 -*-
"""
Beam Plotter
T Overton
3/10/23

Read in a beam file (either before or after simulation)
Plot variables including histogram of each variable
Print beam properties
"""

from PyQt5.QtWidgets import QWidget
import PyQt5.QtWidgets as pyqt
import pyqtgraph as pg
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
    
class BeamPlotWindow(QWidget):
    """
    Load in a beam - plot or print variables, drift beam
    """
    
    def __init__(self):
        super().__init__()
        
        self.setWindowTitle("Plot Beam Properties")
        windowWidth = 1000
        windowHeight = 600
        self.setMinimumSize(windowWidth, windowHeight)
        
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
        self.LoadBeam.setFixedWidth(int(0.25*windowWidth))
        self.BeamsLayout.addWidget(self.LoadBeam,0,0,1,2)
        
        self.BeamList = pyqt.QListWidget()
        self.BeamList.addItems(self.beams)
        self.BeamList.setFixedWidth(int(0.25*windowWidth))
        self.BeamList.setMinimumHeight(int(0.2*windowHeight))
        self.BeamsLayout.addWidget(self.BeamList,1,0,1,2)
        
        self.BeamLabelText = pyqt.QLineEdit()
        self.BeamLabelText.setPlaceholderText('Beam Label')
        self.BeamLabelText.setFixedWidth(int(0.12*windowWidth))
        self.BeamsLayout.addWidget(self.BeamLabelText,2,0)
        
        self.ChangeLabelButton = pyqt.QPushButton('Change Label')
        self.ChangeLabelButton.clicked.connect(self.ChangeLabel)
        self.ChangeLabelButton.setFixedWidth(int(0.12*windowWidth))
        self.BeamsLayout.addWidget(self.ChangeLabelButton,2,1)
        
        self.RemoveBeamButton = pyqt.QPushButton('Remove Beam')
        self.RemoveBeamButton.clicked.connect(self.RemoveBeam)
        self.RemoveBeamButton.setFixedWidth(int(0.25*windowWidth))
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
        self.BeamPropertySelect.setFixedSize(int(0.15*windowWidth), int(0.2*windowHeight))
        self.PropertyPrintLayout.addWidget(self.BeamPropertySelect, 0, 0, 2, 1)
        self.PrintPropertyScale = pyqt.QComboBox()
        self.PrintPropertyScale.addItems(['unit','kilo', 'mega', 'giga','milli', 'micro', 'nano', 'pico', 'femto'])
        self.PrintPropertyScale.setFixedWidth(int(0.09*windowWidth))
        self.PropertyPrintLayout.addWidget(self.PrintPropertyScale,0,1,1,1)
        
        self.PrintPropertyButton = pyqt.QPushButton('Print Property')
        self.PrintPropertyButton.setFixedWidth(int(0.09*windowWidth))
        self.PrintPropertyButton.clicked.connect(self.PrintProperty)
        self.PropertyPrintLayout.addWidget(self.PrintPropertyButton,1,1,1,1)
        self.PropertyPrinted = pyqt.QListView()
        
        self.PropertyPrintBox = pyqt.QListWidget()
        self.PropertyPrintBox.setFixedWidth(int(0.25*windowWidth))
        self.PropertyPrintBox.setMinimumHeight(int(0.2*windowHeight))
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
        self.XAxisHistogram.setFixedHeight(int(windowHeight*0.3))
        self.XAxisHistogram.setMinimumWidth(int(windowWidth*0.48))
        self.XHistValues = [[[],[0]]]
        #Remember this one is rotated!
        self.XAxisHist = pg.PlotCurveItem( self.XHistValues[0][1],self.XHistValues[0][0], stepMode=True, fillLevel=0)
        self.XAxisHistogram.addItem(self.XAxisHist)
        
        self.YAxisHistogram = pg.PlotWidget()
        self.YAxisHistogram.setBackground('w')
        self.YAxisHistogram.setLabel('bottom','Intensity [arb.]')
        self.YAxisHistogram.setMinimumWidth(int(windowWidth*0.24))
        self.YAxisHistogram.setMinimumHeight(int(windowHeight*0.65))
        self.YHistValues = [[[],[0]]]
        #Remember this one is rotated!
        self.YAxisHist = pg.PlotCurveItem( self.YHistValues[0][1],self.YHistValues[0][0], stepMode=True, fillLevel=0)
        self.YAxisHistogram.addItem(self.YAxisHist)
        
        self.ScatterPlot = pg.PlotWidget()
        self.ScatterPlot.setBackground('w')
        self.ScatterPlot.setMinimumSize(int(windowWidth*0.48), int(windowHeight*0.65))
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