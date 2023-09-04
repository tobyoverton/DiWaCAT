"""
@author Tom Pacey
"""
import numpy as np
from scipy.stats import gaussian_kde
from scipy.ndimage import gaussian_filter
from scipy.interpolate import interp2d

from fastkde import fastKDE
from Python_Tools.Modules import general_tools
#from Modules import general_tools

#TODO - add up items form DW discussion
# gradient blur
# PSF from CLARA YAG report
# Noise Equivalent charge - i.e. base charge on pix that is noise floor
# Noise floor modelling - random noise? Gaussian Noise?
# Digitisation method - need to think about how to modify background here - flow could be wrong
# The minus sign issue for vertical - some kind of dict control needed?
# coord should be based on how images are read in

#TODO
# Pixel intesnity histograms
# Is saturation too aggressive - i.e. if 1 hot pixel with beam, then drops intensity too much...
# Some sort of saturation factor control?


#TODO
# Export methods - controlled PNG full_kick_output (scaling and colour!)
# HDF5 in same structure as a CLARA camera image

# Saturation process:
# https://stackoverflow.com/questions/19666626/replace-all-elements-of-python-numpy-array-that-are-greater-than-some-value

#TODO method that adds beam to screen, and clears beam from screen
class BeamToScreen(object):

    def __init__(self,screenDict,beam):

        """
        :param screenDict: dict that decribes the screen properties. require parameters...
        :param beam: a beam object, must have properties beam.x beam.y
        :type beam GeneralBeam
        :type screenDict: dict
        """

        self.beam = beam

        self.required_keys = ["n_pix_x","n_pix_y","pix_size_x","pix_size_y","central_pix_x","central_pix_y",
                              "beam_pix_pad","base_constant","beam_pix_pad","physical_res","KDE_method"]

        general_tools.check_dict_params(screenDict,self.required_keys,"screen data","beamToScreen")

        #TODO this can come from a BG file as well
        #Bg file can be PNg or other
        self.Nx = screenDict["n_pix_x"]
        self.Ny = screenDict["n_pix_y"]


        self.dx = screenDict["pix_size_x"]
        self.dy = screenDict["pix_size_y"]
        self.cen_x = screenDict["central_pix_x"]
        self.cen_y = screenDict["central_pix_y"]

        self.build_base_arrays()

        # TODO implement random noise, array of other vals (images!) etc
        self.base_constant = screenDict["base_constant"]
        self.set_background()

        self.beam_pixel_pad = screenDict["beam_pix_pad"]

        self.set_beam_grid()

        #TODO this process could be other KDE types with switch control

        if "KDE_method" in screenDict:
            self.KDE_method = screenDict["KDE_method"]
        else:
            self.KDE_method = "auto"

        print("Evaluating KDE")
        self.evaluate_beam_KDE()

        #TODO implement different smoothing filters
        #TODO implement different res in x and y
        self.phys_res_x = screenDict["physical_res"]
        self.phys_res_y = screenDict["physical_res"]

        if "post_process_dict" not in screenDict:
            self.post_process_dict = None
        else:
            self.post_process_dict = screenDict["post_process_dict"]

        #used for blurring, doesn't need to be defined?
        self.pix_res_x = self.phys_res_x * self. dx
        self.pix_res_y = self.phys_res_y * self. dy

        #TODO add blurr for image response DOF
        #This is like a physical YAG limit
        self.apply_physical_res_limit()


        #Add the beam pix to the screen pix
        #find overlaps
        self.xi_overlaps,self.yi_overlaps = self.get_pix_overlap_inds(xx_big=self.xx_vals,yy_big=self.yy_vals,xx_small=self.beam_xx,yy_small=self.beam_yy)



        #TODO this will break if background isn't zeroes!
        # Lots of redundant scaling here, but not sure what full_kick_output is best...

        #Just raw full_kick_output
        self.pix_vals_raw = self.add_beam_to_background(sorted_index_x=self.xi_overlaps, sorted_index_y=self.yi_overlaps, z_small=self.beam_processed_pix, z_big=self.bckgrnd)

        # TODO add blurr for image response DOF

        self.pix_vals_processed = self.post_process_image_pix()


        #TODO add a method here for post processing of image whole image


        #Scaled 0 - 1
        self.pix_vals_normalised = self.normalise_pix(self.pix_vals_raw)

        # Scaled so sum = 1
        self.pix_vals_pdf = self.pix_vals_normalised / np.sum(self.pix_vals_normalised)

        # Scaled by charge of input beam
        self.pix_vals_q_dens = self.pix_vals_pdf * self.beam.total_charge

        self.pix_vals_digitised = None


    def build_base_arrays(self):
        self.x_vals = np.arange(0, self.Nx) - self.cen_x
        self.x_vals = self.x_vals * self.dx

        self.y_vals = np.arange(0, self.Ny) - self.cen_y
        self.y_vals = self.y_vals * self.dy

        self.xx_vals,self.yy_vals = np.meshgrid(self.x_vals,self.y_vals)

    #TODO method for different background, from image file, add salt and pepper noise etc.
    def set_background(self):
        self.bckgrnd = np.full(self.xx_vals.shape,self.base_constant)

    def set_beam_grid(self):
        """

        :return:
        """
        beam_x_min = np.min(self.beam.x)
        beam_y_min = np.min(self.beam.y)

        beam_x_max = np.max(self.beam.x)
        beam_y_max = np.max(self.beam.y)

        ind_x_min = self.find_nearest_index(self.x_vals,beam_x_min)
        ind_x_min -= self.beam_pixel_pad

        ind_y_min = self.find_nearest_index(self.y_vals,beam_y_min)
        ind_y_min -= self.beam_pixel_pad

        ind_x_max = self.find_nearest_index(self.x_vals,beam_x_max)
        ind_x_max += self.beam_pixel_pad

        ind_y_max = self.find_nearest_index(self.y_vals,beam_y_max)
        ind_y_max += self.beam_pixel_pad

        ind_x_min,ind_x_max,ind_y_min,ind_y_max = self.correct_indices(ind_x_min,ind_x_max,ind_y_min,ind_y_max)

        self.beam_x_vals = self.x_vals[ind_x_min:ind_x_max]
        self.beam_y_vals = self.y_vals[ind_y_min:ind_y_max]


        self.beam_xx,self.beam_yy = np.meshgrid(self.beam_x_vals,self.beam_y_vals)

        self.beam_grid_coords = np.append(self.beam_xx.reshape(-1,1),self.beam_yy.reshape(-1,1),axis=1)



    def evaluate_beam_KDE(self):

        if self.KDE_method == "scipy_gaussian":
            xy_data = np.vstack((self.beam.x,self.beam.y))
            self.beam_KDE = gaussian_kde(xy_data)
            beam_pix_vals = self.beam_KDE(self.beam_grid_coords.T)
            self.beam_pix_vals_pdf = beam_pix_vals.reshape(len(self.beam_y_vals), len(self.beam_x_vals))

        elif self.KDE_method == "fastKDE" or self.KDE_method == "auto":
            fastPDF, axes = fastKDE.pdf(self.beam.x, self.beam.y)
            x_ax,y_ax = axes
            interped_PDF = interp2d(x_ax,y_ax,fastPDF)
            self.beam_pix_vals_pdf = interped_PDF(self.beam_x_vals, self.beam_y_vals)

        else: #default method - Do i want to have as fastKDe or SciPy...fastKDe needs some C
            fastPDF, axes = fastKDE.pdf(self.beam.x, self.beam.y)
            x_ax,y_ax = axes
            interped_PDF = interp2d(x_ax,y_ax,fastPDF)
            self.beam_pix_vals_pdf = interped_PDF(self.beam_x_vals, self.beam_y_vals)





    def find_nearest_index(self,array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return idx


    #TODO mildly hacky approach but does the job
    def correct_indices(self,xmin,xmax,ymin,ymax):
        """
        correcting for when macros are our of screen range and self.find_nearest_index screws the pooch
        We are hacking a problem here in np.argmin i think
        :param xmin:
        :param xmax:
        :param ymin:
        :param ymax:
        :return:
        """

        #hacky corrector code

        if xmin <0:
            xmin=0

        if ymin <0:
            ymin = 0

        if xmax > len(self.x_vals):
            xmax = -1

        if ymax > len(self.y_vals):
            ymax = -1


        return xmin,xmax,ymin,ymax

    def get_pix_overlap_inds(self,xx_big,yy_big,xx_small,yy_small):
        #For x
        index = np.argsort(xx_big[0])
        sorted_xx_big = xx_big[0][index]
        sorted_index_x = np.searchsorted(sorted_xx_big, xx_small[0])

        #repeat for y
        index = np.argsort(yy_big[:, 0])
        sorted_yy_big = yy_big[:, 0][index]
        sorted_index_y = np.searchsorted(sorted_yy_big, yy_small[:, 0])

        return sorted_index_x,sorted_index_y

    def add_beam_to_background(self,sorted_index_x,sorted_index_y,z_big,z_small):

        for i, xi in enumerate(sorted_index_x):  # we know where to find the overlapping columns in x
            for j, yj in enumerate(sorted_index_y):  # same for y

                z_big[yj,xi] += z_small[j,i]  # the index in big corresponds to the value at the index in small

        return z_big

    def normalise_pix(self,pixel_vals):

        return pixel_vals / np.max(pixel_vals)

    #TODO - rename and rework - optical PSF and then YAG PSF!!
    # property that is the type of smoothing
    # control structure for post process
    # test and validate postprocessing
    def apply_physical_res_limit(self):
        self.beam_processed_pix = gaussian_filter(self.beam_pix_vals_pdf, sigma=[self.pix_res_x, self.pix_res_y])

    def post_process_image_pix(self):
        """

        :return: post processed pixels based on other post processing methods. Has a custom gradient blurr
        """

        accepted_methods = ["GaussianGradient"]

        if self.post_process_dict == None: # Just do nothing
            print("No post processing method defined")
            return self.pix_vals_raw

        if "method" not in self.post_process_dict:
            print("method key is not defined in post processing dict")
            raise KeyError
        elif self.post_process_dict["method"] not in accepted_methods:
            print("method key is not one of the accepted methods")
            raise ValueError

        if self.post_process_dict["method"] == "GaussianGradient":
            return gradient_gauss_blur()






    def set_digitise_pix_by_max_value(self,bit_depth,max_pix_val):
        '''
        Digitise the normalsied pixels to a set peak pixel value
        If set to maximum of bit depth, this will allow for exact saturation
        :param bit_depth:
        :type bit_depth: int
        :param max_pix_val: the pixel value for maximum. Manual set
        :type max_pix_val: int
        :return: none - setter function
        '''

        max_pix = 2**bit_depth

        if max_pix_val > max_pix:
            raise Exception("Maximum pixel value cannot be higher than allowed by bit depth")

        pix_bins = np.arange(0,max_pix)

        scaled_pix = self.pix_vals_normalised * max_pix

        print("scaled max pixel is: ",np.max(scaled_pix))

        self.pix_vals_digitised = np.digitize(scaled_pix, pix_bins,right=True)

    def set_digitise_pix_by_charge_dens(self,bit_depth,saturation_charge_dens):
        '''
        Digitise the normalsied pixels based on a saturation charge density (best defined and extracted from another test screen)

        :param bit_depth:
        :param saturation_charge_dens: the charge density that will hvae the "saturated value" at the given bit depth
        :return: none - setter
        '''

        max_pix = 2 ** bit_depth

        pix_bins = np.arange(0, max_pix)

        #Saturated the charge density
        q_dens = np.minimum(self.pix_vals_q_dens,saturation_charge_dens)

        #q_dens[q_dens > saturation_charge_dens] = saturation_charge_dens

        #Scale the charge density
        print("saturation charge set =",saturation_charge_dens)
        print("saturation charge in array = ",np.max(q_dens))


        #Vital to do this correctly, takes into accoutn if under saturated!
        q_dens_normalised = q_dens / saturation_charge_dens

        q_dens_scaled = q_dens_normalised * max_pix

        print("scaled max pix is ", np.max(q_dens_scaled))

        self.pix_vals_digitised = np.digitize(q_dens_scaled, pix_bins,right=True)



def gradient_gauss_blur():
    pass




if __name__ == "__main__":
    print("Hello")
