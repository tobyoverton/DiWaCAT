"""
confirming functionality of gradient blur function
uses a uniform grid target form Thorlabs website
does a horizontal and vertical gradient blur from centre out
"""


from Python_Tools.Modules import general_tools as gent
import matplotlib.pyplot as plt

from PIL import Image
import numpy as np

test_file = "../Test_files/import_files/Distortion_Concentric_Circle_A1-780.jpg"
set_image = Image.open(test_file)
set_image.convert('L')
image_array = np.array(set_image)[:, :, 0]


test_dict_h_blur = {
    "num_slices":20,
    "split_axis":1,
    "max_sig":20,
    "base-blur":0,
    "process_blur":5
}

test_dict_v_blur = test_dict_h_blur.copy()
test_dict_v_blur["split_axis"] = 0

processed_image_h_blur = gent.gradient_gauss_blur(image_array,test_dict_h_blur)
processed_image_v_blur = gent.gradient_gauss_blur(image_array,test_dict_v_blur)


fig,ax = plt.subplots(1,3)
plt.gray()

ax[0].imshow(set_image)
ax[1].imshow(processed_image_h_blur)
ax[2].imshow(processed_image_v_blur)




plt.show()