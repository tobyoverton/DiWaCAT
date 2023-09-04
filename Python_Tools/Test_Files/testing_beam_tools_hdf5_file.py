"""
Another great thingy for testing beam tools with hdf5 beam file
"""


import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const

from Python_Tools.Modules import beam_tools as dbt

import matplotlib.axes._axes as axes # this for pyCharm completion for axes objects
from matplotlib.lines import Line2D # this for custom legends


#from Python_Tools.Modules import read_beam_file as rbf
##

test_file = "import_files/EBT-BA1-COFFIN-FOC.hdf5"
# open a beam file and get it
myBeam = dbt.BeamFromFile()

myBeam.read_HDF5_beam_file(filename=test_file)
print('read a beam file from ', myBeam._beam['code'])

# Same test code as for Gaussian beam proceedure
myPlotter = dbt.BeamPlotter()

myPlotter.addBeamToPlotter('ImportedBeam',myBeam)

myPlotter.plotLPS(alpha = 0.1)
myPlotter.plotLPS_t(alpha = 0.1)


myPlotter.plot_transverse_properties(alpha = 0.1)

myPlotter.plot_physical_projections(alpha = 0.1)



plt.show()
#plt.close()

