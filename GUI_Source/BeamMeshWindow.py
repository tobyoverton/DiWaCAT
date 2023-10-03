# -*- coding: utf-8 -*-
"""
Beam Meshing GUI window
T Overton
3/10/23

Read in an existing beam, perform meshing and save the file
"""

from PyQt5.QtWidgets import QWidget
import PyQt5.QtWidgets as pyqt
from pathlib import Path
import scipy.constants as const

from Python_Tools.Modules import beam_tools as dbt



def UnitConversion(prefix):
    if prefix == 'unit' or prefix == 'Unit':
        return 1
    else:
        return 1/getattr(const, prefix)

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
            #Change z values so mean is at zero
            BeamIn.z += (-1) * BeamIn.Mz
            
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
        
