# -*- coding: utf-8 -*-
"""
DiWaCAT GUI using PyQT
"""

# Only needed for access to command line arguments

from PyQt5.QtWidgets import QWidget
from PyQt5.QtWidgets import QApplication, QMainWindow
import PyQt5.QtWidgets as pyqt
from PyQt5.QtCore import QThread,pyqtSignal
import sys, os, subprocess
from pathos.helpers import mp
import scipy.constants as const
import numpy as np
import pyqtgraph as pg
from pathlib import Path
import copy
from Python_Tools.Modules import diwacat_tools as DWA
from Python_Tools.Modules import beam_tools as dbt

#Force the working directory to be the root directory of DiWaCAT files
os.chdir(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(os.path.dirname(__file__))

#Import the files with each GUI window and functions
from GUI_Source import BeamMakeWindow, BeamMeshWindow, FieldPlotter, BeamPlotter, FieldCalcWindow#, DiWaCAT_Tracker

#Issue with multithreading when using imported file. Temporary fix to include the functions in the GUI file
def UnitConversion(prefix):
    if prefix == 'unit' or prefix == 'Unit':
        return 1
    else:
        return 1/getattr(const, prefix)
   
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
            if __name__ == '__main__':
                pool = mp.Pool(10)
                args = [(CurrentBeam, pos) for pos in MacroPositions]
                KicksApplied = pool.starmap(WorkerFunction_Interpolation, args)       
                pool.close()
                pool.join()
            print(KicksApplied[100])
            KicksApplied = np.array(KicksApplied)
            xKick = np.array(KicksApplied[:,0,0])
            yKick = np.array(KicksApplied[:,1,0])
            zKick = np.array(KicksApplied[:,2,0])
            #Drift the beam to the end of the step and then apply the 
            
            CurrentBeam.driftBeam(0.5 * StepLength)
            #CurrentBeam.beam._beam['px'] += StepLength * KicksApplied[:,0] * CurrentBeam.beam.q_over_c
            #CurrentBeam.beam._beam['py'] += StepLength * KicksApplied[:,1] * CurrentBeam.beam.q_over_c
            #CurrentBeam.beam._beam['pz'] += StepLength * KicksApplied[:,2] * CurrentBeam.beam.q_over_c
            CurrentBeam.beam._beam['px'] += StepLength * xKick * CurrentBeam.beam.q_over_c
            CurrentBeam.beam._beam['py'] += StepLength * yKick * CurrentBeam.beam.q_over_c
            CurrentBeam.beam._beam['pz'] += StepLength * zKick * CurrentBeam.beam.q_over_c
            
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
    def stop(self):
        self._isRunning = False
            
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
        self.stop_thread()
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
    def stop_thread(self):
        #Make sure our thread is closed - attempt to stop memory leaks
        self.function.stop()
        self.function.quit()
        self.function.wait()
    def DriftBeamEnd(self):
        beamToDrift=DWA.DiWaCATOutput()
        beamToDrift.beam = copy.deepcopy(self.TrackedBeams[-1])
        driftDistance = self.DriftDistance.value() * 1/UnitConversion(self.DriftDistanceScale.currentText())
        beamToDrift.driftBeam(driftDistance)
        self.TrackedBeams.append(beamToDrift.beam)
        self.BeamList.addItem('Post Drift')
        
        self.LengthArray.append(self.DLWLength.value() * 1/UnitConversion(self.DLWLengthScale.currentText()) + driftDistance)
        self.SigmaX.append(beamToDrift.beam.Sx)
        self.SigmaY.append(beamToDrift.beam.Sy)
        self.Mx.append(beamToDrift.beam.Mx)
        self.My.append(beamToDrift.beam.My)
        self.epsx.append(beamToDrift.beam.normalized_horizontal_emittance)
        self.epsy.append(beamToDrift.beam.normalized_vertical_emittance)
        self.Charge.append(100 * beamToDrift.beam.total_charge * 1/self.TrackedBeams[0].total_charge)
        self.SigmaT.append(beamToDrift.beam.Sz * 1/const.c)
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

#-------- Definition of the main window -----------#
   
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
       #Temporary fix. Same reason as the definition of the tracing window in the main file
        #self.f = DiWaCAT_Tracker.BeamTrackWindow()
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
if __name__ == '__main__':
    app = QApplication(sys.argv)
    w = MainWindow()
    w.show()
    app.exec()
