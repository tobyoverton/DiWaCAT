"""

Some general tools that don't belong in classes for specific uses

@author T Pacey
"""

import os
import math
import numpy as np

from PIL import Image

from scipy.ndimage import uniform_filter
from scipy.ndimage import gaussian_filter

def control_file_overwrite_logic(filename,overwrite_existing):
    """
    check if a file exists, and return logic on whether it should be overwritten, for various full_kick_output functions
    Should be used in some full_kick_output file function
    Saves copypasta
    :param filename: file to check
    :param overwrite_existing:
    :return: True if to continue with full_kick_output process, False if full_kick_output function should be stopped
    """

    continue_flag = None

    print("Controlling file logic for file at: ",filename)
    print("Overwriting file has been set to: ",overwrite_existing)

    if not overwrite_existing:
        if os.path.isfile(filename):
            print("File found, overwriting not selected, sending flag to stop file write process")
            continue_flag = False
        else:
            print("File not found, sending flag to continue file write process")
            continue_flag = True

    else:
        if os.path.isfile(filename):
            print("File found, overwriting is selected, sending flag to continue file write process")
            continue_flag = True

        else:
            print("File not found, sending flag to continue file write process")
            continue_flag = True

    return continue_flag

def check_file_string_has_type_extension(file):
    """

    :param file: a file name or file complete file directory
    :return: if a file string is identified via os.path.splittext returns true, additional full_kick_output of that extension if True
    """
    name, extension = os.path.splitext(file)

    if extension == '':
        return False, None
    else:
        return True,extension


def check_dict_params(test_dict,required_keys,dict_name = "<unknown>",output_type= "<undefined>"):
    """
    Check if a dictionary hold all the keys you want it to
    To help with print statements and/or logging - will print strings of the dict_name and the full_kick_output type
    The full_kick_output type cna be the class or the usage case for this run of check_dict_params
    :param test_dict: dict to be tested
    :param required_keys: list of keys to be tested
    :return: True if checks passed
    """

    missedKeys = []
    for key in required_keys:
        if key not in test_dict:
            missedKeys.append(key)

    if len(missedKeys) != 0:
        print("Missing the following keys from dict for " + dict_name)
        print(missedKeys)
        print("Output type " + output_type + " cannot be generated")
        raise KeyError

    print("Keys validated for dict: "+dict_name + " for use with "+output_type)
    return True


def gradient_gauss_blur(image,setup_dict):
    """

    :param image: image array (format?!) to apply bluur to
    :param setup_dict: dict that holds the vlaues for performing teh custom bluring
    :return: gray scale (!) gaussian blurr image
    :type image: numpy array
    :type setup_dict dict
    """
    required_vals = ["num_slices", "split_axis", "max_sig"]

    check_dict_params(setup_dict,required_vals,dict_name="image blur dict",output_type="gradient blur")

    # Optional tweaks
    if "base_blur" in setup_dict:
        base_blur = setup_dict["base_blur"]
    else:
        base_blur = 1

    if "process_blur" in setup_dict:
        process_blur = setup_dict["process_blur"]
    else:
        process_blur = 10

    #Check dimensions of input image

    if len(np.shape(image)) == 2:
        image_array = np.array(image)
    else:
        print("Converting image colour space")
        set_image = Image.fromarray(image)
        set_image.convert('L')
        image_array = np.array(set_image)[:, :, 0]

    num_slices = setup_dict["num_slices"]
    max_sig = setup_dict["max_sig"]
    split_axis = setup_dict["split_axis"]

    SIG_RATE = max_sig / num_slices

    imageslices = np.array_split(image_array, num_slices, axis=split_axis)

    sig_mask = np.arange(0, len(imageslices))
    sig_mask_split = np.array_split(sig_mask, 2)

    sig_vals_right_half = sig_mask_split[0] * SIG_RATE
    sig_vals_left_half = np.flip(sig_vals_right_half)

    SIG_VAL = np.concatenate([sig_vals_left_half, sig_vals_right_half])

    for i, slice in enumerate(imageslices):

        if split_axis == 0:
            imageslices[i] = gaussian_filter(slice, sigma=[SIG_VAL[i], base_blur], mode='nearest')
            reformedimage = np.concatenate(imageslices, axis=split_axis)
            processedimage = uniform_filter(reformedimage, size=[process_blur,0], mode='constant')

        elif split_axis == 1:
            imageslices[i] = gaussian_filter(slice, sigma=[base_blur, SIG_VAL[i]], mode='nearest')
            reformedimage = np.concatenate(imageslices, axis=split_axis)
            processedimage = uniform_filter(reformedimage, size=[0, process_blur], mode='constant')

        else:
            print("Could not define splitting axis, should be either 0 or 1")
            raise ValueError



    return processedimage



def weighted_avg_and_std_dev(values, weights):
    """
    Return the weighted average and standard deviation.

    from https://stackoverflow.com/questions/2413522/weighted-standard-deviation-in-numpy

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return average, np.sqrt(variance)

def weighted_std_dev(values,weights,average):
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return np.sqrt(variance)

def get_closest_larger_square(value):
    """
    based on this https://stackoverflow.com/questions/49875299/find-nearest-square-number-of-a-given-number/49876112
    handles edge case for when a square number is already provided
    :param value: number to gte clsoest largest square for
    :return: integer of closest larest square
    """

    assert type(value) == int, "Value must be an integer"

    root_of_n = math.sqrt(value)

    floor_integer = int(root_of_n)
    ceil_integer = floor_integer + 1


    floor_integer_square = floor_integer * floor_integer

    ceil_integer_square = ceil_integer * ceil_integer
    if floor_integer_square == value:
        return floor_integer_square
    else:
        return ceil_integer_square

def find_nearest(array,value):
    """
    from https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
    :param array:
    :param value:
    :return:
    """
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return array[idx-1]
    else:
        return array[idx]


def get_increment_list_to_value(value,increment):
    """
    This allows there to be a floating point error in the division of steps. this can compound if many steps are added up
    :param value: the value that the sum should go to
    :param increment: increment for the list
    :return: list of increments, where sum is equal to value, final increment is remainder
    """

    assert value > increment,ValueError

    remainder = value % increment

    quotient = value // increment

    #This is the safest way to check for a floating point error afaik
    clean_divisor = (quotient * increment) == value

    quotient_int = int(quotient)

    increment_array = np.full(quotient_int,increment)

    #this stops there being a small flaoting point error on the final step but ensuring we only add the remaidner when it is "large"
    if not clean_divisor:
        #add that final step
        #print("not a clean divisor, adding a final step to get to end of element")
        increment_array = np.append(increment_array,remainder)
        #assert np.sum(increment_array)==value

    return increment_array


def get_percentile_values(values,percentile = 95):
    """
    Centres beam on mean of values, reflects about 0, finds percentile Y that covers all macros within that percentile
    :param values: beam dimension to find RMS value of for given percentile
    :param percentile: percentile to find around cenrted to mean. Default is 95
    :return:
    """

    mean = np.mean(values)
    centerted_vals = values - mean
    reflected_vals = np.abs(centerted_vals)
    percentile_point = np.percentile(reflected_vals,percentile)
    left_bound = -1*percentile_point
    right_bound = percentile_point

    subset_array =  values[(values>left_bound) & (values<right_bound)]

    return subset_array