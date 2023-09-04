"""
Another great thingy for testing beam tools with wakecode  beam file
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

test_file = "import_files/200fs_CLARA_PartialLoss.h5"


# open a beam file and get it
myBeam = dbt.BeamFromFile()

myBeam.read_WakeCode_beam_file(filename=test_file,)
print('read a beam file from ', myBeam._beam['code'])

print('total charge is ', myBeam.total_charge / const.pico)
print('num macros', myBeam.nMacros)

myPlotter = dbt.BeamPlotter()

myPlotter.addBeamToPlotter('ImportedWakeCode',myBeam)

myPlotter.plotLPS(alpha = 0.1)

myPlotter.plot_transverse_properties(alpha = 0.1)

myPlotter.plot_physical_projections(alpha = 0.1)

myPlotter.general_scatter_plot(x_val='x',x_scale='micro',y_val='cpz',y_scale='mega',alpha=0.1)

myPlotter.general_multi_plot(x_vals=['x','y'],x_scales=['micro','micro'],y_vals=['cpz','xp'],y_scales=['mega','milli'],alpha=0.1)


plt.show()
#plt.close()

