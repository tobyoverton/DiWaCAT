"""
@author T Pacey

1) read beam dist from hdf5 file
2) Generate image from beam dist sand save to file
3) Import Image file
4) Analyse image file and produce a profile graph
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const

from Python_Tools.Modules import beam_tools as btools
from Python_Tools.Modules import beam_to_screen as b2s
from Python_Tools.Modules import image_tools as im_tools


test_file = "import_files/EBT-BA1-COFFIN-FOC.hdf5"
output_image_file = "export_files/Simulated_BA1_IP.hdf5"
PIX_SIZE = 32 * const.micro
CEN_PIX_Y = 600

#region Read beam dist from file

#---------------------------------
#READ BEAM FROM FILE
#-------------------------------

print("----- IMPORTING BEAM FROM FILE---------")

# open a beam file and get it
myBeam = btools.BeamFromFile()

myBeam.read_HDF5_beam_file(filename=test_file)
print('read a beam file from ', myBeam._beam['code'])
print(f"Number of macros is {myBeam.nMacros}")
myPlotter = btools.BeamPlotter()

myPlotter.addBeamToPlotter('ImportedBeam',myBeam)
myPlotter.general_density_plot(x_val='x',x_scale='milli',y_val='y',y_scale='milli')
myPlotter.general_scatter_plot(x_val='x',x_scale='milli',y_val='y',y_scale='milli')
plt.show()
#endregion



#region beam to image
print("--------SIMULATING A SCREEN ----------")

#Establish a screen dict
# Some example reasonable parameters

screenDict = {}
screenDict["n_pix_x"] = 1936
screenDict["n_pix_y"] = 1216
screenDict["pix_size_x"] = PIX_SIZE
screenDict["pix_size_y"] = PIX_SIZE
screenDict["central_pix_x"] = 800
screenDict["central_pix_y"] = CEN_PIX_Y
screenDict["base_constant"] = 0
screenDict["beam_pix_pad"] = 5

screenDict["physical_res"] = 10 * const.micro

#Optional
screenDict["KDE_method"] = "fastKDE"

#For digitisation
bit_depth = 12
peak_pix = 2**bit_depth # this sets us to "perfectly saturated"

testScreen = b2s.BeamToScreen(screenDict,myBeam)
testScreen.set_digitise_pix_by_max_value(bit_depth=bit_depth, max_pix_val=peak_pix)

plt.figure(1)
plt.pcolormesh(testScreen.beam_x_vals/const.milli, testScreen.beam_y_vals/const.milli, testScreen.beam_processed_pix, shading='auto',cmap = 'turbo')
plt.title("The pixels that are generated from the beam in physical space")
plt.xlabel("x (mm)")
plt.ylabel("y (mm)")

plt.figure(2)
plt.imshow(testScreen.pix_vals_digitised,aspect=testScreen.x_vals.ptp()/testScreen.y_vals.ptp(),cmap = 'jet',clim=(0,peak_pix))
plt.title("The simulated image in \"pixel space\"")


print("--------Saving image to file----------")
#First generate an image object
image_to_output = im_tools.ImageFromArray(array=testScreen.pix_vals_digitised,bit_depth=bit_depth,origin='Lower')

plt.figure(3)
plt.imshow(image_to_output.raw_image)
plt.title("This image will be saved")
image_to_output.save_to_hdf5(output_image_file)


#plt.show()



#endregion

#------------------IMPORTING IMAGE----------------------

#region Importing image

imported_image = im_tools.ImageFromFile(file_directory=output_image_file,bit_depth=12)
imported_image.read_image_from_HDF5()

plt.figure(4)
plt.imshow(imported_image.processed_image)
plt.title("this image has been imported")

plt.figure(5)
plt.plot(imported_image.y_pix_array,imported_image.y_projection)
plt.title("Projection in pixel space")

#normalise to central pixel. convert to real space. convert to m. convert to millimeters
physical_coords = (imported_image.y_pix_array - CEN_PIX_Y) * PIX_SIZE * (1/const.micro) * (const.milli)

plt.figure(6)
plt.plot(physical_coords,imported_image.y_projection)
plt.xlabel(" y (mm)")
plt.ylabel("Sum Counts")
plt.title("Projection in real space")
#endregion

plt.show()