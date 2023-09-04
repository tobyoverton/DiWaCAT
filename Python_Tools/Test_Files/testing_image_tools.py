"""
@author T Pacey

For testing functionality of image_tools module
"""


import os
import sys
import numpy as np
import matplotlib.pyplot as plt

from scipy.ndimage import gaussian_filter,median_filter,uniform_filter
from ruamel.yaml import YAML, representer
import h5py

from Python_Tools.Modules import image_tools as im_tools

#Little section for getting an image to open from my test YAML file -

#TODO change this reading section
#region opening from my YAML test file - not needed really
#TODO tidy up, different file types, reference test HDF5 needed

TEST_KEY = "iteration_5"
TEST_IMPORT_FILE = "../Test_files/import_files/BA1_slit_scan_metadata.yaml"

IMAGE_SERVER = "//claraserv3"

#specify the file and check it exists
assert os.path.isfile(TEST_IMPORT_FILE), "File does not exist"

#specify the yaml interpreter mode
yaml = YAML(typ='safe')

#open the YAML file and print everything
# with open(TEST_IMPORT_FILE, "r") as stream:
#     try:
#         print(yaml.load(stream))
#     except yaml.yamlError as exc:
#         print(exc)

with open(TEST_IMPORT_FILE, "r") as stream:
    meta_data = yaml.load(stream)
test_image_location_raw = meta_data[TEST_KEY]["image_location"]


test_image_location = test_image_location_raw
test_image_location = test_image_location.replace('//','/')



test_image_file = IMAGE_SERVER+test_image_location+".hdf5"

#endregion

test_image_file = "import_files/test_BA1_IP_jet.png"

test_image = im_tools.ImageFromFile(file_directory=test_image_file,bit_depth=12,origin='upper')

#print(test_image.get_all_HDF5_keys())

test_image.read_image_from_PNG(pilmode='P')

print(test_image.meta_data)

print(test_image.dimensions_pix)

#print(test_image.meta_data['hdf5_attributes'])

plt.figure(1)
plt.imshow(test_image.raw_image,cmap = 'jet',origin=test_image.origin)

y_mask = (450,650)
x_mask = (800,1000)

masked_test_image = im_tools.MaskedImage(test_image)

masked_test_image.mask_by_pixels(x_mask,y_mask)

plt.figure(2)
plt.imshow(masked_test_image.raw_image,cmap = 'jet',origin=masked_test_image.origin)

masked_test_image.perform_standard_cleanup(threshold_value=0)

plt.figure(3)
plt.imshow(masked_test_image.processed_image,cmap = 'jet',origin=masked_test_image.origin)

print("Mean pixel in x {:.1f}, standard deviation in x {:.1f}".format(masked_test_image.mean_pix_x,masked_test_image.std_dev_pix_x))
print("Mean pixel in y {:.1f}, standard deviation in y {:.1f}".format(masked_test_image.mean_pix_y,masked_test_image.std_dev_pix_y))

#plot the projections on mean centred axis
plt.figure(4)
plt.plot(masked_test_image.x_pix_array - masked_test_image.mean_pix_x,masked_test_image.x_projection)
plt.plot(masked_test_image.y_pix_array - masked_test_image.mean_pix_y,masked_test_image.y_projection)

plt.show()

