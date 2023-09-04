"""
@author T Pacey
"""

import os
import sys
import math
import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import interp1d

from ruamel.yaml import YAML
yaml = YAML(typ='safe')

from Python_Tools.Modules import image_tools as im_tools
from Python_Tools.Modules import general_tools as gent


class SlitScanAnalyser(object):
    """
    For analysis of slit scan data files
    takes single slit scan meta file and imports images for generating phase space plots and projections
    """
    def __init__(self, import_file,projection_plane):

        #create attributes
        self.meta_data_file = None
        self.meta_data = None
        self.all_slit_positions = []
        self.selected_slit_positions = None

        self.raw_slit_images = [] #clear-able
        self.slit_images = [] #clear-able

        #for phase space
        self.projection = None
        self.projected_counts = []
        self.u_prime_vals = []
        self.common_u_prime = None
        self.counts_on_common_u_prime = []
        self.u_coords = None
        self._zero_slit_position = None

        self._u_vals = None
        self._u_prime_vals = None

        #check file exists, then assign to attribute
        assert os.path.isfile(import_file), "Import file does not exist"
        self.meta_data_file = import_file

        #check we can have a projection plane
        if projection_plane not in ['x','y']:
            raise ValueError("Projection plane must be str of x or y")

        #set that plane
        self.projection = projection_plane

        #open the file and load the meta data
        print("Opening meta file at", self.meta_data_file)
        with open(self.meta_data_file, "r") as stream:
           self.meta_data = yaml.load(stream)

        print("Completed loading meta file data")

        for key in self.meta_data:
            self.all_slit_positions.append(self.meta_data[key]["slit_position"])

        #TODO - the units here or make a property of this

        self.all_slit_positions = sorted(self.all_slit_positions)

    def select_slit_positions(self,start_iter=0,end_iter=-1,step=1):
        self.selected_slit_positions = self.all_slit_positions[start_iter:end_iter:step]

    def clear_slit_images(self):
        self.slit_images = []
        #Destrcutor needed for the image objects?

    def clear_raw_slit_images(self):
        self.raw_slit_images = []
        # Destrcutor needed for the image objects?


    def read_in_raw_slit_images(self,image_server,bit_depth=16, origin='lower',key_to_open="Capture000001"):
        for pos in self.selected_slit_positions:
            for key in self.meta_data:
                pos_dict = self.meta_data[key]["slit_position"]
                if pos == pos_dict:
                    print("Accessing iteration : {}".format(key))
                    print(f"This should be slit position : {pos:.3f}")
                    image_loc = self.meta_data[key]['image_location']
                    image_file = image_server + image_loc + ".hdf5"
                    #print(image_file)
                    image = im_tools.ImageFromFile(file_directory=image_file,bit_depth=bit_depth,origin=origin)
                    image.read_image_from_HDF5(key_to_open=key_to_open)

                    self.raw_slit_images.append(image)

    def clean_up_images_and_threshold_images(self,rough_x_mask,rough_y_mask,percentile_threshold = 99):
        for image in self.raw_slit_images:
            masked_array = np.ma.array(image.raw_image,mask=False)
            masked_array.mask[rough_y_mask[0]:rough_y_mask[1],
                                rough_x_mask[0]:rough_x_mask[1]] = True

            threshold = np.percentile(masked_array,percentile_threshold)
            image.perform_standard_cleanup(threshold_value=threshold)



    def auto_mask_slit_images(self, pixel_padding = 0, x_mask=None,y_mask=None,rough_mask_x=None,rough_mask_y = None):
        for image in self.raw_slit_images:
            self.clean_up_images_and_threshold_images(rough_x_mask =rough_mask_x,rough_y_mask=rough_mask_y)
            masked_image = im_tools.MaskedImage(image)
            masked_image.auto_mask_on_peak(pixel_padding=pixel_padding,
                                           forced_x_mask=x_mask, forced_y_mask=y_mask,
                                           min_mask_x=rough_mask_x,min_mask_y=rough_mask_y)

            self.slit_images.append(masked_image)

    def mask_slit_images(self,x_mask,y_mask):
        for image in self.raw_slit_images:
            masked_image = im_tools.MaskedImage(image)
            masked_image.mask_by_pixels(x_mask, y_mask)

            self.slit_images.append(masked_image)

    def process_slit_images(self,threshold,**kwargs):
        for image in self.slit_images:
            image.perform_standard_cleanup(threshold_value=threshold,**kwargs)

    def rotate_slit_images(self,angle):
        for image in self.slit_images:
            image.perform_rotation(angle=angle)

# HERE we go to calcualting for phase space already

    def calc_divergences(self,screen_calib_factor,distance_slits2screen):

        if self.projection == 'x':
            scr_uvals = np.arange(0, len(self.slit_images[0].x_projection)) * screen_calib_factor

        elif self.projection == 'y':
            scr_uvals = np.arange(0, len(self.slit_images[0].y_projection)) * screen_calib_factor

        else:
            raise RuntimeError("Projection not recognised")

        #loop over each image and assign divergences
        for i, image in enumerate(self.slit_images):

            # CRUCIAL to select projection here
            if self.projection == 'x':
                projected_count = image.x_projection

            elif self.projection == 'y':
                projected_count = image.y_prection

            else:
                raise RuntimeError("Projection not recognised")

            u_prime = (scr_uvals - self.selected_slit_positions[i]) / distance_slits2screen

            #The positiosn of stage are in mm, the distance to screen must be mm, the calibration must be in mm
            # Therefore the divergence is rads, so we will convert back to mrads
            u_prime = u_prime * 1e3  # mrad

            self.u_prime_vals.append(u_prime)
            self.projected_counts.append(projected_count)


    def set_zero_slit_position(self,zero_position_index=None):

        if zero_position_index is None:
            zero_position_index = int(len(self.selected_slit_positions) / 2)

        self._zero_slit_position = self.selected_slit_positions[zero_position_index]

    def convert_to_common_u_prime_coords(self):
        """
        Each slit image has its own u' coord system for each slit position from the image
        These need to interpolated to a common coord system for each image
        :return:
        """

        self.common_u_prime = self.u_prime_vals[self.zero_slit_index]

        #for each set of projections, find an interpoaltion back to the common coord system
        # this is just resampling
        for i, u_prime in enumerate(self.u_prime_vals):
            x_vals = u_prime
            y_vals = self.projected_counts[i]
            f = interp1d(x_vals, y_vals,kind = "cubic", bounds_error=False, fill_value=(0.0,0.0))
            self.counts_on_common_u_prime.append(f(self.common_u_prime))

    def set_u_coords(self,flip_coords):

        if flip_coords:
            self.u_coords = np.flip(self.selected_slit_positions)
        else:
            self.u_coords = self.selected_slit_positions

        self.u_coords = self.u_coords - self.zero_slit_position

    def set_coord_mesh(self,flip_u_coords=True):

        self.set_u_coords(flip_coords=flip_u_coords)
        self.convert_to_common_u_prime_coords()

        self._u_vals, self._u_prime_vals = np.meshgrid(self.u_coords, self.common_u_prime)

    def shift_coord_mesh(self,shift_u,shift_u_prime):
        self._u_vals = self._u_vals - shift_u
        self._u_prime_vals = self._u_prime_vals - shift_u_prime

    def save_phase_space_plot_data(self,output_file,overwrite_existing_file = False,output_format=".npy"):

        has_ext,extension = gent.check_file_string_has_type_extension(output_file)

        if has_ext:
            complete_file_output = output_file
        else:
            complete_file_output = output_file + output_format
            extension = output_format

        print(f"Outputting phase space plot data to format {extension}")
        print(f"Outputting phase space plot data to {output_file}")

        if gent.control_file_overwrite_logic(complete_file_output, overwrite_existing=overwrite_existing_file):

            if extension == '.npy':
                with open(complete_file_output, 'wb') as f:
                    np.save(f, self.u)
                    np.save(f, self.u_prime)
                    np.save(f, self.intensities)

            elif extension == '.csv':
                with open(complete_file_output,'wb') as f:
                    output_array = np.column_stack((np.ravel(self.u),np.ravel(self.u_prime),np.ravel(self.intensities.T)))
                    print(output_array.shape)
                    print(output_array[0])
                    np.savetxt(f,output_array,delimiter=",",header=self.projection+' (mm)'+','+self.projection+"' (mrad)"+","+"Intensity (au)")

            print("Output completed")

#---------------------------------------------------
# Properties
#----------------------------------------------------

    @property
    def zero_slit_position(self):
        return self._zero_slit_position

    @property
    def zero_slit_index(self):
        return self.selected_slit_positions.index(self.zero_slit_position)

    @property
    def intensities(self):
        return np.asarray(self.counts_on_common_u_prime)

    @property
    def u(self):
        return self._u_vals

    @property
    def u_prime(self):
        return self._u_prime_vals

# METHODS in module for analysis
# Not the bets way to do this but can't see wodd for trees without this
#NOT THE ANALYSER CLASS

def plot_all_images_in_file(input_file,image_server,projection,x_mask=(0,-1), y_mask = (0,-1),cmap='turbo'):

    assert os.path.isfile(input_file), "File does not exist"
    all_images_in_scan = SlitScanAnalyser(input_file,projection)
    all_images_in_scan.select_slit_positions()#Selects all

    all_images_in_scan.read_in_raw_slit_images(image_server=image_server)

    all_images_in_scan.clean_up_images_and_threshold_images(rough_x_mask=x_mask,rough_y_mask=y_mask,percentile_threshold=99)

    all_images_in_scan.mask_slit_images(x_mask=x_mask,y_mask=y_mask)


    total_slit_images = len(all_images_in_scan.raw_slit_images)
    print(f"Total number of images = {total_slit_images}")

    # #TODO this is a very useful tool!
    # # making a plot of all the slit images
    plot_nums = np.arange(0, total_slit_images)
    square = gent.get_closest_larger_square(total_slit_images)
    side_len = int(np.sqrt(square))
    plot_square = np.append(plot_nums, np.full(square - len(plot_nums), np.nan))
    plot_square = np.reshape(plot_square, (side_len, side_len))

    fig = plt.figure()

    gs = fig.add_gridspec(side_len, side_len, hspace=0, wspace=0)
    axs = gs.subplots(sharex='col', sharey='row')

    for i,image in enumerate(all_images_in_scan.slit_images):
        squareloc = np.where(plot_square == i)
        axs[int(squareloc[0]), int(squareloc[1])].imshow(image.processed_image,cmap=cmap)
        #axs[int(squareloc[0]), int(squareloc[1])].set_title(f'i = {i:.0f}')

    return fig,axs


def process_images_to_output_file(input_file,image_server,projection,
                                  first_image=0,final_image=-1,skip_images=1,
                                  base_x_mask=None,base_y_mask=None,pixel_pad=0,rotation_angle = 0,
                                  forced_x_mask=None,forced_y_mask=None,
                                  output_file = None,overwrite_output=False,):

    scan = SlitScanAnalyser(input_file,projection)
    scan.select_slit_positions(first_image,final_image,skip_images)

    scan.read_in_raw_slit_images(image_server=image_server)

    print("image dimensions",scan.raw_slit_images[0].dimensions_pix)

    scan.auto_mask_slit_images(pixel_padding=pixel_pad,
                               x_mask=forced_x_mask,y_mask=forced_y_mask,
                               rough_mask_x=base_x_mask,rough_mask_y=base_y_mask)
    scan.rotate_slit_images(rotation_angle)
    total_slit_images = len(scan.slit_images)

    print(f"Total number of images = {total_slit_images}")

    fig1 = plt.figure()
    if projection == 'x':
        for image in scan.slit_images:
            plt.plot(image.x_pix_array + image.x_global_shift, image.x_projection)

    elif projection == 'y':
        for image in scan.slit_images:
            plt.plot(image.y_pix_array + image.y_global_shift, image.y_projection)

    # #TODO this is a very useful tool!
    # # making a plot of all the slit images
    plot_nums = np.arange(0, total_slit_images)
    square = gent.get_closest_larger_square(total_slit_images)
    side_len = int(np.sqrt(square))
    plot_square = np.append(plot_nums, np.full(square - len(plot_nums), np.nan))
    plot_square = np.reshape(plot_square, (side_len, side_len))


    fig2 = plt.figure()
    gs = fig2.add_gridspec(side_len, side_len, hspace=0, wspace=0)
    axs = gs.subplots(sharex='col', sharey='row')

    for i,image in enumerate(scan.slit_images):
        squareloc = np.where(plot_square == i)
        axs[int(squareloc[0]), int(squareloc[1])].imshow(image.processed_image,cmap='turbo')


    slit_pos_output = scan.selected_slit_positions

    total_counts_output = []
    centroids_output = []
    std_dev_output = []

    for image in scan.slit_images:
        total_counts_output.append(image.total_count)
        centroids_output.append(image.mean_pix_x + image.x_global_shift)
        std_dev_output.append(image.std_dev_pix_x)

    print(centroids_output)

    np.asarray(total_counts_output)

    #save those arrays
    if gent.control_file_overwrite_logic(output_file,overwrite_existing=overwrite_output):
        with open(output_file,'wb') as f:
            print("starting output of data to file:",output_file)
            np.save(f,np.asarray(slit_pos_output))
            np.save(f,np.asarray(total_counts_output))
            np.save(f,np.asarray(centroids_output))
            np.save(f,np.asarray(std_dev_output))

            print("Output to file  complete")


    return fig1,fig2


def calculate_emittance_from_processed_file(input_file,distance_to_screen,screen_calib_factor,proj,beam_gamma = 70,output_file = None,overwrite_output = False):
    # specify the file and check it exists
    assert os.path.isfile(input_file), "File does not exist"

    image_dict = {}

    with open(input_file, 'rb') as f:
        slit_positions = np.load(f)
        # slit_pos.sort()
        input_counts = np.load(f)
        input_centroids = np.load(f)
        input_std_devs = np.load(f)

    for i, pos in enumerate(slit_positions):
        image_dict[pos] = {}
        image_dict[pos]["total_count"] = input_counts[i]
        image_dict[pos]["centroid"] = input_centroids[i] * screen_calib_factor
        image_dict[pos]["std_dev"] = input_std_devs[i] * screen_calib_factor

    print(image_dict)


    Ntot = 0
    av_u_list = []

    # Establish Ntot, av_x,av_xp

    for pos in image_dict:
        nj = image_dict[pos]['total_count']
        Ntot += nj

        xsj = pos
        av_u_list.append(nj * xsj)

    av_x_raw = sum(av_u_list) / Ntot

    print(f"Raw <x> = {av_x_raw:.5f}")

    mid_slit_pos = gent.find_nearest(list(image_dict.keys()), av_x_raw)

    print(f"Nearest measured slit pos to centroid = {mid_slit_pos:.5f}")
    print(f"Difference = {mid_slit_pos - av_x_raw:.5f}")

    mid_centroid_pos = image_dict[mid_slit_pos]['centroid']


    correction_factor = av_x_raw

    Ntot = 0
    av_x_list = []
    av_xp_list = []

    # Establish Ntot, av_x,av_xp

    for pos in image_dict:
        nj = image_dict[pos]['total_count']
        Ntot += nj

        xsj = pos - correction_factor

        # Screen centroid for beamlet - eqn 9
        Xbar_j = image_dict[pos]['centroid'] - mid_centroid_pos

        # Beamlet divergence - eq 10
        xp_bar_j = (Xbar_j - xsj) / distance_to_screen

        av_xp_list.append(nj * xp_bar_j)
        av_x_list.append(nj * xsj)

    av_x = sum(av_x_list) / Ntot
    av_xp = sum(av_xp_list) / Ntot

    # av_x = 0
    # av_xp=0

    av_x2_list = []
    av_xp2_list = []
    av_x_xp_list = []
    for pos in image_dict:
        nj = image_dict[pos]['total_count']
        xsj = pos - correction_factor

        Xbar_j = image_dict[pos]['centroid'] - mid_centroid_pos

        # Beamlet divergence - eq 10
        xp_bar_j = (Xbar_j - xsj) / distance_to_screen

        sig_j = image_dict[pos]['std_dev']

        print(f"at pos {pos:.2f}, xsj =  {xsj:.2f},Xbar = {Xbar_j:.1e}mm & sig_j = {sig_j:.1e} mm")

        sig_xp_j = sig_j / distance_to_screen

        av_x2_list.append(nj * ((xsj - av_x) ** 2))
        av_xp2_list.append((nj * (sig_xp_j ** 2)) + (nj * ((xp_bar_j - av_xp) ** 2)))
        av_x_xp_list.append((nj * xsj * xp_bar_j) - (Ntot * av_x * av_xp))

    av_x2 = sum(av_x2_list) / Ntot
    av_xp2 = sum(av_xp2_list) / Ntot
    av_x_xp = sum(av_x_xp_list) / Ntot

    print(f"<"+proj+"> = {av_x:.5f}")

    print(f"<"+proj+"^2> = {av_x2:.5f}")

    print(f"<"+proj+"'> = {av_xp:.5f}")

    print(f"<"+proj+"'^2> = {av_xp2:.3e}")

    print(f"<"+proj+""+proj+"'> = {av_x_xp:.3e}")

    print(f"emit_g ^2 = {(av_x2 * av_xp2) - (av_x_xp ** 2)}")

    try:
        emit_x = math.sqrt((av_x2 * av_xp2) - (av_x_xp ** 2))

        print(f"Geometric emittance is = {emit_x:.3e} mm rad")

        #print(f"Normalised emittance is = {emit_x * 70 * 1e3:.4f} mm mrad")
    except:
        print("Emittance could not be calucalted")

    # endregion

    # Quick and dirty emittance calc
    print("-------------------")

    print("Q&D sanity check:")

    beam_width_by_slits = max(image_dict.keys()) - min(image_dict.keys())

    rms_size_at_slits = beam_width_by_slits / 4
    calc_rms_size_at_slits = math.sqrt(av_x2)

    print(f"Beam width at slits = {beam_width_by_slits:.1f}")
    print(f"RMS size from width estimate = {rms_size_at_slits:.1f}")
    print(f"RMS size from calclation = {calc_rms_size_at_slits}")

    rms_size_central_slit = image_dict[mid_slit_pos]['std_dev']

    q_and_d_emit = (rms_size_at_slits * rms_size_central_slit) / distance_to_screen

    print(f"Geometric Q&D emittance is = {q_and_d_emit:.3e} mm rad")

    print("------------")

    print(f"Normalised emittance is = {emit_x * beam_gamma * 1e3:.1f} mm mrad")
    print(f"Normalised Q&D emittance is = {q_and_d_emit * beam_gamma * 1e3:.1f} mm mrad")

    print("------------")

    if output_file is not None:
        print("Writing to output file")
        if gent.control_file_overwrite_logic(output_file,overwrite_existing=overwrite_output):
            with open(output_file, 'w') as outfile:
                outfile.write("--------------------\n")
                outfile.write(f"Screen calibration factor = {screen_calib_factor} \n")
                outfile.write(f"Distance from slit to screen = {distance_to_screen} \n")
                outfile.write(f"Beam gamma value = {beam_gamma} \n")
                outfile.write("-------------------- \n")
                outfile.write(f"Raw geometric RMS emittance x calculated {emit_x} \n")
                outfile.write("-------------------- \n")
                outfile.write(f"Normalised emittance  is = {emit_x * beam_gamma * 1e3:.1f} mm mrad\n")
                outfile.write(
                    f"Normalised Q&D  emittance is = {q_and_d_emit * beam_gamma * 1e3:.1f} mm mrad \n")
                outfile.write("--------------------")


    return emit_x,q_and_d_emit