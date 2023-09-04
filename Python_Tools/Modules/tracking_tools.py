"""
Module for tracking of particle beams from beam tools standard file
@author T Pacey
"""

import os
import numpy as np
import scipy.constants as constants
import copy

from Python_Tools.Modules import beam_tools
from Python_Tools.Modules import general_tools as gent




class LatticeTracker(object):
    def __init__(self, beam: beam_tools.GeneralBeam, lattice: list):
        """

        :param beam: a beam_tools General beam type
        :param lattice: a list of lattice elements (dicts)
        """

        self._beam = beam
        self._initial_beam = self._beam.returnBeam() # copies a beam over for comparative purposes

        self._lattice = lattice

        self._beam.check_beam_keys()

        self.parameters_to_save = {}

        self.defined_elements = [] #TODO write - need a dict with the required args for each element

        self.saved_beams = []
        self.save_beam_element_end = False

        self._global_step_value = None

        self.zero_length_elements = ['save_beam','thin_quad']


    #region General boilerplate housekeeping
    def check_lattice(self):
        print("Checking lattice element list")
        for element in self._lattice:
            for element_name in element:
                assert hasattr(self,element_name),"undefined element for lattice"


    def clear_lattice(self):
        self._lattice = None

    def set_lattice(self,lattice):
        self._lattice = lattice


    def clear_saved_parameters(self):
        self._parameters = None

    def clear_inital_beam(self):
        """
        If you want to save memory
        :return: None
        """
        self._initial_beam = None

    def set_inital_beam(self):
        self._initial_beam = self._beam.returnBeam()

    #endregion

    #region tracking functions

    def set_saved_parameters(self,parameters: list):
        """
        :param parameters: list of beam parameters that can be measured and then saved
        :return:
        """
        for parameter in parameters:
            assert hasattr(self._beam,parameter),"parameter must be defined in the beam"
            self.parameters_to_save[parameter] = []


    def get_steps_in_element(self,element_params):

        if "step_size" in element_params:
            return gent.get_increment_list_to_value(value=element_params['L'], increment=element_params["step_size"])

        elif self.global_step_value is not None:
            return gent.get_increment_list_to_value(value=element_params['L'], increment=self.global_step_value)

        else:
            return [element_params['L']]

    def perform_track(self):

        self.check_lattice()

        for element in self._lattice:
            print(f"Tracking in element: {element}")
            for element_name in element:
                element_parameters = element[element_name]

                if element_name not in self.zero_length_elements:
                    #We have to do trackingover its length in steps defined by element, globally, or just to the end of the element
                    dz_steps = self.get_steps_in_element(element_parameters)

                    for dz in dz_steps:
                        temp_parameter_dict = element_parameters # a dict to hold the small values of L for the steps
                        temp_parameter_dict['L'] = dz # set L to what we want
                        getattr(self, element_name)(**temp_parameter_dict) # magical kwargs

                        for parameter in self.parameters_to_save:
                            self.parameters_to_save[parameter].append(getattr(self._beam,parameter))

                else: #element must be of zero length
                    getattr(self, element_name)(**element_parameters)  # magical kwargs
                    for parameter in self.parameters_to_save:
                        self.parameters_to_save[parameter].append(getattr(self._beam, parameter))

                if self.save_beam_element_end:
                    self.save_beam()
    #endregion

    #region The elements

    def save_beam(self):
        self.saved_beams.append(self._beam.returnBeam())

    def drift(self,L,step_size = None):
        """
        Simple drift, no ballistic compression
        :param L: drift length in m
        :param step_size: needed only for step size setting per element, does not act in function
        :return:
        """
        self._beam._beam['x'] = self._beam.x + L * self._beam.xp
        self._beam._beam['y'] = self._beam.y + L * self._beam.yp
        self._beam._beam['z'] = self._beam.z + L
        # xp,yp,pz doesn't change

    def thin_quad(self, f):
        """
        Thin les approximation of quad, F>0 focusing horizontally
        :param f: focal length
        :param step_size:
        :return:
        """
        temp_xp = (-1 * self._beam.x / f) + self._beam.xp
        temp_yp = (self._beam.y / f + self._beam.yp)

        # convert back to our set variables of px, py
        self._beam._beam['px'] = self._beam._angle_to_mom(temp_xp, self._beam.pz)
        self._beam._beam['py'] = self._beam._angle_to_mom(temp_yp, self._beam.pz)

    def thick_quad(self,K,L,step_size=None):
        """
        Thick quad, K>0 means focusing applied in X, defocusing applied in y
        :param K: Quad strength
        :param L: quad length
        :return: none, sets x,y,px,py,z based on quad
        """
        if K >= 0:
            new_x, new_xp = self.__transportMatrixThickFQuad(self._beam.x, self._beam.xp, K, L)
            new_y, new_yp = self.__transportMatrixThickDQuad(self._beam.y, self._beam.yp, K, L)
        else:
            K = -1 * K
            new_x, new_xp = self.__transportMatrixThickDQuad(self._beam.x, self._beam.xp, K, L)
            new_y, new_yp = self.__transportMatrixThickFQuad(self._beam.y, self._beam.yp, K, L)

        # Apply to beam
        self._beam._beam['x'] = new_x
        self._beam._beam['px'] = self._beam._angle_to_mom(new_xp, self._beam.pz)

        self._beam._beam['y'] = new_y
        self._beam._beam['py'] = self._beam._angle_to_mom(new_yp, self._beam.pz)

        self._beam._beam['z'] = self._beam.z + L

    def __transportMatrixThickFQuad(self, u, uprime, K, L):
        """
        Take canonical variables u, u' and transform through thick focusing quad of strength K, length L
        :param u: x or y
        :param uprime: x' or y'
        :param K: quad strength
        :param L: length
        :return: transformed u,u'
        """

        phi = np.sqrt(K) * L
        mat = np.array([[np.cos(phi), (1 / np.sqrt(K) * np.sin(phi))],
                        [(-1 * np.sqrt(K) * np.sin(phi)), np.cos(phi)]])

        output = np.matmul(mat, np.array([u, uprime]))

        return output[0], output[1]

    def __transportMatrixThickDQuad(self, u, uprime, K, L):
        """
        Take canonical variables u, u' and transform through thick DE-focusing quad of strength K, length L
        :param u: x or y
        :param uprime: x' or y'
        :param K: quad strength
        :param L: length
        :return: transformed u,u'
        """

        K = np.abs(K)  # first clean up K if used as negative
        phi = np.sqrt(K) * L
        mat = np.array([[np.cosh(phi), (1 / np.sqrt(K) * np.sinh(phi))],
                        [(1 * np.sqrt(K) * np.sinh(phi)), np.cosh(phi)]])

        output = np.matmul(mat, np.array([u, uprime]))

        return output[0], output[1]

    #endregion


    #properties

    @property
    def global_step_value(self):
        return self._global_step_value

    @global_step_value.setter
    def global_step_value(self,step):
        assert step>0,"Steps must be larger than 0"
        self._global_step_value = step

    @property
    def saved_parameters(self):
        return self.parameters_to_save