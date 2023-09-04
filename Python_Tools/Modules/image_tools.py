"""
@author T Pacey

Set of tools for importing and analysing images
"""

import os
import sys
import numpy as np
from Python_Tools.Modules import general_tools
from scipy.ndimage import gaussian_filter, median_filter, uniform_filter
from scipy.ndimage import rotate as rotate_image
from scipy.signal import peak_widths
import h5py

# TODO can this be lived without??
from imageio import imread


class GeneralImage(object):
    """
    Parent Class that sets out the rules and properties of an image so it can be controlled and manipulated
    Child classes - image from array, image from hdf5, image from png (colour maps/LUT!)
    Image subset from array (masked image) - subselecting an image down and it having its properties recast
    Holds all the image processing functions, and  processing macros/method lists for the array
    Properties - metadata dict, image size, size in x, size in y, first moments, second moments.
    Projection in x, y , FWHM, FWXM, FWQM, FW-ANY-M,pixel histogram
    Gaussian fits
    """

    def __init__(self, bit_depth, origin):

        self._img_array = None
        self._proc_img_array = None

        self.processing_applied = False

        # comment strings for labelling, meta_data is a dict to hold extra parameters
        self._meta_data_dict = {}
        self._comment_string = ""

        self._origin = origin  # as in matplotlib pixel orientation
        self._bit_depth = bit_depth

    def check_properties(self):
        print("Checking image properties generated correctly from parameters")

        allowed_origins = ['lower', 'upper']
        if self.origin not in allowed_origins:
            print("specified origin not in allowed type of: ", allowed_origins)
            raise ValueError

        if self.bit_depth % 2 != 0:
            print("bit depth value of ", self.bit_depth, " does not divide by 2")
            raise ValueError

        print("Checks passed")

    def reset_processed_array(self):
        """
        Resets the processed array to the raw array and resets the processing flag
        :return: None
        """
        self._proc_img_array = self._img_array
        self.processing_applied = False

    def perform_median_filter(self, size=2):
        """
        performs median on the processed image array, and sets the  processing flag
        :param size: median filter size required by scipy
        :return:
        """
        # If filter applied flag False then the processed image is the raw image, else its the processed values
        # This is needed for chain processing operations
        self._proc_img_array = median_filter(self.processed_image, size)
        self.processing_applied = True

    def perform_thresholding(self, threshold_value):

        self._proc_img_array = (self.processed_image > threshold_value) * self.processed_image
        self.processing_applied = True

    def perform_uniform_filter(self,size):
        """

        :param size: size of uniform filter kernel
        :return:
        """
        # If filter applied flag False then the processed image is the raw image, else its the processed values
        # This is needed for chain processing operations
        self._proc_img_array = uniform_filter(self.processed_image, size,mode='constant')
        self.processing_applied = True

    def perform_gaussian_smoothing(self, sigma):
        """

        :param sigma:
        :return:
        """
        self._proc_img_array = gaussian_filter(self.processed_image, sigma)
        self.processing_applied = True

    def perform_standard_cleanup(self, threshold_value,reset_processing=True, median_size=2, gaussian_sigma=3):

        #first check if reset processing for this "macro"
        if reset_processing:
            self.reset_processed_array()

        # first x-rays
        self.perform_median_filter(size=median_size)

        # second threshold
        self.perform_thresholding(threshold_value)

        # finally clean up noise
        self.perform_gaussian_smoothing(sigma=gaussian_sigma)

        # belts and braces
        self.processing_applied = True

    def perform_rotation(self,angle):

        # using scipy ndimage rotate as rotate_image
        self._proc_img_array = rotate_image(self.processed_image,angle=angle,reshape=False)
        self.processing_applied = True


    #TODO some work is needed here
    def save_to_hdf5(self,filename,dataset_name='Capture000001',meta_data='',overwrite_existing=False):

        #Make sure extension is correct if not specified
        filename_info = os.path.splitext(filename)
        filename_output = filename_info[0] + '.hdf5'


        if general_tools.control_file_overwrite_logic(filename,overwrite_existing):

            with h5py.File(filename_output, "w") as f:
                dataset = f.create_dataset(dataset_name, data=self.processed_image)
                dataset.attrs["meta_data"] = meta_data


    # Properties

    @property
    def raw_image(self):
        return self._img_array

    @property
    def processed_image(self):
        # vital check for chain processing, if any processing has been applied then the processing_image must be used
        if self.processing_applied:
            return self._proc_img_array
        else:
            return self._img_array

    @property
    def meta_data(self):
        return self._meta_data_dict

    @property
    def comment(self):
        return self._comment_string

    @property
    def bit_depth(self):
        return self._bit_depth

    @property
    def origin(self):
        return self._origin

    @property
    def dimensions_pix(self):
        return self.raw_image.shape

    # TODO Calibrations!

    @property
    def x_dim_pix(self):
        return self.dimensions_pix[1]

    @property
    def y_dim_pix(self):
        return self.dimensions_pix[0]

    # Dervied properties

    @property
    def average_pixel_intensity(self):
        return np.average(self.processed_image)

    @property
    def total_count(self):
        return np.sum(self.processed_image)

    @property
    def x_projection(self):
        return np.sum(self.processed_image, axis=0)

    @property
    def y_projection(self):
        return np.sum(self.processed_image, axis=1)

    @property
    def pos_peaK_val_x_projection(self):
        return np.argmax(self.x_projection)

    @property
    def pos_peaK_val_y_projection(self):
        return np.argmax(self.y_projection)

    @property
    def x_pix_array(self):
        return np.arange(0, self.x_dim_pix)

    @property
    def y_pix_array(self):
        return np.arange(0, self.y_dim_pix)

    @property
    def mean_pix_x(self):
        return np.average(self.x_pix_array, weights=self.x_projection)

    @property
    def mean_pix_y(self):
        return np.average(self.y_pix_array, weights=self.y_projection)

    @property
    def std_dev_pix_x(self):
        return general_tools.weighted_std_dev(self.x_pix_array, self.x_projection, self.mean_pix_x)

    @property
    def std_dev_pix_y(self):
        return general_tools.weighted_std_dev(self.y_pix_array, self.y_projection, self.mean_pix_y)


class ImageFromFile(GeneralImage):
    def __init__(self, file_directory, bit_depth, origin="lower"):
        super().__init__(bit_depth=bit_depth, origin=origin)

        # apply general attributes

        self.file_directory = file_directory
        assert os.path.isfile(file_directory), "Cannot read image from file as File does not exist"

        self.extension = os.path.splitext(self.file_directory)[1]
        assert self.file_type != 'unknown', "file type is not yet handled by this object"

        self.check_properties()

    def read_image_from_PNG(self,pilmode):
        """

        :param pilmode: imread pilmode - for grayscale use 'L', if jet us 'P', not sure what others could be used
        :return:
        """
        #TODO what is pilmode for?

        assert self.file_type == "png","File type is not png, use other methods"

        read_image = imread(self.file_directory,pilmode=pilmode)

        self._meta_data_dict = read_image.meta
        self._img_array = np.array(read_image)




    def read_image_from_HDF5(self, key_to_open=None):

        assert self.file_type == "hdf5", "File type is not hdf5, use other methods"

        with h5py.File(self.file_directory, "r") as h5file:

            # Find what key to open
            if key_to_open == None:
                print("Key not specified, opening first dataset key")
                key_to_open = list(h5file.keys())[0]
                print("Opening: ", key_to_open)

            elif key_to_open in self.get_all_HDF5_keys():
                print("Opening: ", key_to_open)

            else:
                print("key: ", key_to_open, "not found in file")

            self._img_array = np.array(h5file.get(key_to_open))

            # open and format all attrbutes of image into dict
            # init the key for the metadata dict
            self._meta_data_dict['hdf5_attributes'] = {}
            for key in h5file[key_to_open].attrs.keys():
                self._meta_data_dict['hdf5_attributes'][key] = h5file[key_to_open].attrs[key]


    def get_all_HDF5_keys(self):
        assert self.file_type == "hdf5", "File type is not hdf5, use other methods"

        keys_in_file = []
        with h5py.File(self.file_directory, "r") as h5file:
            keys_in_file = list(h5file.keys())
        return keys_in_file

    # Properties

    @property
    def file_type(self):
        if self.extension in [".hdf5", ".h5"]:
            return "hdf5"

        elif self.extension in [".png"]:
            return "png"

        else:
            return "unknown"

        # TIFFs needed to be added


class ImageFromArray(GeneralImage):
    def __init__(self,array,bit_depth,origin):
        super().__init__(bit_depth= bit_depth,origin=origin)
        self._img_array = array


class MaskedImage(GeneralImage):
    def __init__(self, base_image: GeneralImage):
        super().__init__(bit_depth=base_image.bit_depth, origin=base_image.origin)
        # inital flag for applying the mask
        self.mask_applied = False
        self.base_image = base_image

        self.check_properties()

        self.x_min_crop = None
        self.x_max_crop = None
        self.y_min_crop = None
        self.y_max_crop = None

    def mask_by_pixels(self, x_pix_range, y_pix_range):
        """

        :param x_pix_range: tuple(x_min,x_max)
        :param y_pix_range: tuple(y_min,y_max)
        :return: masked applied as True if all code run

        """

        if self.mask_applied:
            print("Mask has already been applied, try clearing mask to reapply")
            return False

        assert len(x_pix_range) == 2, "x range must be length 2 tuple"
        assert len(y_pix_range) == 2, "y range must be length 2 tuple"

        #check if specifying min max with 0 -> -1
        if not (x_pix_range[0] == 0 and x_pix_range[1] ==-1):
            assert x_pix_range[0] < x_pix_range[1], "tuple must be of format min,max"

        if not (y_pix_range[0] == 0 and y_pix_range[1] ==-1):
            assert y_pix_range[0] < y_pix_range[1], "tuple must be of format min,max"

        x_min = x_pix_range[0]
        x_max = x_pix_range[1]

        y_min = y_pix_range[0]
        y_max = y_pix_range[1]

        # TODO factor in dimensions as well to check limits, can't mask on something larger than dimensions

        self._img_array = self.base_image._img_array[y_min:y_max, x_min:x_max]

        # Carry over any processing from the base image
        # carry the array of processed image
        if self.base_image._proc_img_array is not None:
            self._proc_img_array = self.base_image._proc_img_array[y_min:y_max, x_min:x_max]

        self.processing_applied = self.base_image.processing_applied  # carry the flag!

        self.mask_applied = True

        self.x_min_crop = x_min
        self.x_max_crop = x_max

        self.y_min_crop = y_min
        self.y_max_crop = y_max

        return self.mask_applied

    def clear_mask(self):
        self._img_array = self.base_image._img_array
        self._proc_img_array = self.base_image._proc_img_array

        self.mask_applied = False

    def mask_by_size(self):
        #No idea what i meant this to be now...
        pass

    def auto_mask_on_peak(self, min_mask_x=None,min_mask_y=None, pixel_padding = (0,0), forced_x_mask = None, forced_y_mask = None):
        """
        :param min_mask_x
        :param min_mask_y
        :param pixel_padding: additional pixels to pad around range found by autocropper, if int applies to both dims, if tuple applies to (x,y)
        :param forced_x_mask: if forcing a mask in x dimension use this with a tuple min,max
        :param forced_y_mask: if forcing a mask in y dimension use this with a tuple min,max
        passes to self.mask_by_pixels the pixels to mask over
        :return:
        """
        # # Storms method from here
        # https://gitlab.stfc.ac.uk/adig/python/clara-ba1-beamsize/-/blob/master/BeamSizes.py

        #set padding values for each dimension

        if type(pixel_padding) is int:
            pixel_pad_x = pixel_padding
            pixel_pad_y = pixel_padding
        elif type(pixel_padding) is tuple:
            pixel_pad_x = pixel_padding[1]
            pixel_pad_y = pixel_padding[0]


        if self.mask_applied:
            print("Mask has already been applied, try clearing mask to reapply")
            return False

        x_min = None
        x_max = None
        y_min = None
        y_max = None

        #TODO logic here is a bit garbled and is copypasty for both dimensions, needs rewriting as a method and run for each dim
        #deal with forced cases
        if forced_x_mask is not None:
            assert len(forced_x_mask) == 2, "x range must be length 2 tuple"
            assert forced_x_mask[0] < forced_x_mask[1], "tuple must be of format min,max"
            x_min = forced_x_mask[0]
            x_max = forced_x_mask[1]
            x_peak_width = 0 # Set to 0 for logic below

        else:
            print("Finding mask width in x")
            x_peak_width = int(
                peak_widths(self.base_image.x_projection, [self.base_image.pos_peaK_val_x_projection])[0][0])

        if forced_y_mask is not None:
            assert len(forced_y_mask) == 2, "y range must be length 2 tuple"
            assert forced_y_mask[0] < forced_y_mask[1], "tuple must be of format min,max"
            y_min = forced_y_mask[0]
            y_max = forced_y_mask[1]
            y_peak_width = 0 # Set to 0 for logic below


        else:
            print("Finding mask width in y")
            y_peak_width = int(peak_widths(self.base_image.y_projection, [self.base_image.pos_peaK_val_y_projection])[0][0])


        # SPLIT INTO DIMENSIONS X & Y and then check each one
        #X mask setting
        if x_peak_width !=0:
            if x_min is None and x_max is None:
                print( "Applying automatically detected mask to beam")
                x_min = self.base_image.pos_peaK_val_x_projection - x_peak_width - pixel_pad_x
                x_max = self.base_image.pos_peaK_val_x_projection + x_peak_width + pixel_pad_x

        elif min_mask_x is not None:
            print("----------Could not accurately locate peak width in x, setting mask specified minimum -------------")
            x_min = min_mask_x[0]
            x_max = min_mask_x[1]


        else:
            print("----------Could not accurately locate peak in x, no minimum specificed, setting mask to full range -------------")
            x_min = 0
            x_max = self.base_image.x_dim_pix

        #y dimension mask setting
        if y_peak_width != 0:
            if y_min is None and y_max is None:
                y_min = self.base_image.pos_peaK_val_y_projection - y_peak_width - pixel_pad_y
                y_max = self.base_image.pos_peaK_val_y_projection + y_peak_width + pixel_pad_y

        elif min_mask_y is not None:
            print("----------Could not accurately locate peak width in y, setting mask specified minimum -------------")
            y_min = min_mask_y[0]
            y_max = min_mask_y[1]

        else:
            print("----------Could not accurately locate peak in y, no minimum specificed, setting mask to full range -------------")
            y_min = 0
            y_max = self.base_image.y_dim_pix

        #Finally ready to reapply the mask

        self.mask_by_pixels((x_min,x_max),(y_min,y_max))

        return self.mask_applied


    def mask_on_radius(self):
        # Something smarter? mask on ellipse?
        pass

    # region Properties

    @property
    def cropped_x_at(self):
        if self.mask_applied and self.x_min_crop is not None:
            return self.x_min_crop
        else:
            return 0

    @property
    def cropped_x_to(self):
        if self.mask_applied and self.x_max_crop is not None:
            return self.x_max_crop
        else:
            return 0

    @property
    def cropped_y_at(self):
        if self.mask_applied and self.y_min_crop is not None:
            return self.y_min_crop
        else:
            return 0

    @property
    def cropped_y_to(self):
        if self.mask_applied and self.y_max_crop is not None:
            return self.y_max_crop
        else:
            return 0

    @property
    def x_global_shift(self):
        return self.cropped_x_at

    @property
    def y_global_shift(self):
        #TODO this could need to interact with ORIGIN in some way, need to watch how y axis can be flipped
        return self.cropped_y_at


    # endregion

