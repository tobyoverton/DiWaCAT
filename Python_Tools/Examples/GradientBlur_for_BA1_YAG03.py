"""
For testing the gradient blur fucntion on a YAg03 image to show how it blurs the edges of teh reticles
uses a custom reference image
@author Tom Pacey
"""
import numpy as np


from PIL import Image

import matplotlib.pyplot as plt
import Python_Tools.Modules.general_tools as gent

#TODO generalise to an read in image
# Apply to BA1 reference YAG image to tune spec for gradients
# Control for left/right blurring, centre out blurring etc.



ref_image = Image.open('../Test_files/import_files/YAG03 Reference.png')


yag_image = Image.open('../Test_files/import_files/YAG03_Reference_Modded_Reticules.png')
test_image = np.array(yag_image)

test_dict = {
    "num_slices":50,
    "split_axis":1,
    "max_sig":15,
    "base_blur":1,
    "process_blur":10
}

test_output_image = gent.gradient_gauss_blur(test_image,test_dict)

fig,ax1 = plt.subplots(1,2)
plt.gray()

ax1[0].imshow(test_image)
ax1[1].imshow(test_output_image)


fig,ax3 = plt.subplots(1,3)

ax3[0].imshow(ref_image)
ax3[1].imshow(test_image)
ax3[2].imshow(test_output_image)

for axe in ax3:
    axe.set_ylim([300, 750])
    axe.set_xlim([150, 300])

plt.gray()
plt.show()