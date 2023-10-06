# -*- coding: utf-8 -*-
"""
DiWaCAT GUI using PyQT
"""

# Only needed for access to command line arguments

from PyQt5.QtWidgets import QApplication, QMainWindow
import PyQt5.QtWidgets as pyqt
import sys, os, subprocess

#Force the working directory to be the root directory of DiWaCAT files
os.chdir(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(os.path.dirname(__file__))

#Import the files with each GUI window and functions
from GUI_Source import BeamMakeWindow, BeamMeshWindow, FieldPlotter, BeamPlotter, DiWaCAT_Tracker, FieldCalcWindow


   
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
        self.b = BeamMakeWindow.BeamMakeWindow()
        self.b.show()
        
    def BeamMeshOpen(self,checked):
        self.c = BeamMeshWindow.BeamMeshWindow()
        self.c.show()

    def FieldPlotOpen(self, checked):
        self.d = FieldPlotter.FieldPlotWindow()
        self.d.show()
        
    def BeamPlotOpen(self,checked):
        self.e = BeamPlotter.BeamPlotWindow()
        self.e.show()
        
    def SimulateBeamWindow(self,checked):
        self.f = DiWaCAT_Tracker.BeamTrackWindow()
        self.f.show()
        
    def DiWaCATOpen(self,checked):
        self.g = FieldCalcWindow.CalculationWindow()
        self.g.show()
        #popen = subprocess.Popen(["./FieldSolver_Executable/DiWaCAT_FieldSolverUI.exe"])
        
if __name__ == '__main__':
    app = QApplication(sys.argv)
    w = MainWindow()
    w.show()
    app.exec()