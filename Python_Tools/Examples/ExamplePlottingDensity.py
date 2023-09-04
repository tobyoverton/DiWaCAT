#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 15:02:30 2022

@author: mbcx4to2
"""

import os
import os, sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const

import matplotlib.axes._axes as axes # this for PyCharm completion for axes objects
from matplotlib.lines import Line2D # this for custom legends
from matplotlib.font_manager import FontProperties

from Python_Tools.Modules import beam_tools as dbt, read_beam_file as rbf, beam_to_screen as b2s

BeamFile_Input = "import_files/300fs_sigmar150_OnAxis.h5"

Beam_Input = dbt.BeamFromFile()
Beam_Input.read_WakeCode_beam_file(filename=BeamFile_Input)



myPlotter = dbt.BeamPlotter()
myPlotter.addBeamToPlotter('Input Beam',Beam_Input)
plot, axeslist = myPlotter.general_multi_plot(x_vals=('x', 'y'), y_vals=('xp', 'yp'), x_scales=('milli', 'milli'),
                                            y_scales=('milli', 'milli'), alpha=0.4,plotcolors=['red'])

# lines, labels = plot.axes[-1].get_legend_handles_labels()
# leg = plot.legend(lines,labels,loc=(0.2,0.78))
# for lh in leg.legendHandles:
#     lh._legmarker.set_markersize(6)
#     lh._legmarker.set_alpha(1)
# for ax in axeslist:
#     ax.get_legend().remove()

myPlotter2 = dbt.BeamPlotter()
myPlotter2.addBeamToPlotter('Input Beam',Beam_Input)
plot2,axeslist2 = myPlotter2.general_density_plot(x_val='z', y_val='y', x_scale='micro',
                                            y_scale='micro')
# for ax in axeslist2:
#     ax.set_xlim([-2800,700])
#     ax.set_ylim([-300,1000])
plot4, axeslist4 = myPlotter2.general_density_plot(x_val='x', y_val='y', x_scale='micro',
                                            y_scale='micro')
# for ax in axeslist4:
#     ax.set_xlim([-150,150])
#     ax.set_ylim([-250,1000])

plt.show()