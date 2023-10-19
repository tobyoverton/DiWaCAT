"""
Based on J Jones read beam file from SimFrame

class - general beam
class - gaussianBunch - creates your Gaussian bunch
class - beamBinner - bins beams, smooths,ouputs, will work on James file or gaussianBunch
class - beam kicker

Tom P March 2021
"""
import os, time, csv, sys, subprocess

import numpy as np
import scipy.constants as constants
from scipy.ndimage import gaussian_filter, median_filter, uniform_filter
from scipy.stats import skewnorm
import h5py
import matplotlib.pyplot as plt
import copy
from fastkde import fastKDE

from Python_Tools.Modules import Arb_Dist_Sampler as ADS
from Python_Tools.Modules import general_tools

from Python_Tools.Modules import minimumVolumeEllipse as mve

MVE = mve.EllipsoidTool() # Should this be here?


# from scipy.spatial.distance import cdist
# from scipy.spatial import ConvexHull
# from scipy.stats import gaussian_kde



class GeneralBeam(object):
    def __init__(self, particle_mass=constants.m_e):

        self.E0 = particle_mass * constants.speed_of_light ** 2
        self.E0_eV = self.E0 / constants.elementary_charge
        self.q_over_c = (constants.elementary_charge / constants.speed_of_light)
        self.speed_of_light = constants.speed_of_light
        self.E0_MeV = self.E0_eV / constants.mega

        self.required_beam_keys = ['x', 'y', 'z', 'px', 'py', 'pz', 't', 'charge', 'macro_type', 'code']

    def resetDicts(self):
        self._beam = {}
        self.twiss = {}  # TODO not sure what Twiss is needed for as a dict, can't remember

    def check_beam_keys(self, usage="GeneralBeam class"):

        return general_tools.check_dict_params(self._beam, self.required_beam_keys, "beam", usage)

        # TODO add the length checks here

    def export_beam(self):
        pass

    # so all beam files will be written the same way
    # must match read from HDF5 setup so that method will import directly back to the class

    def write_HDF5_beam_file(self, filename, overwrite_existing_file=False, centered=False, mass=constants.m_e,
                             sourcefilename=None, pos=None, rotation=None, longitudinal_reference='t', xoffset=0,
                             yoffset=0, zoffset=0, toffset=0, ):
        """
        Code for outputing to hdf5 file compatible with simframe
        :param filename:
        :param centered:
        :param mass:
        :param sourcefilename:
        :param pos:
        :param rotation:
        :param longitudinal_reference:
        :param xoffset:
        :param yoffset:
        :param zoffset:
        :param toffset:
        :return:
        """

        complete_keys = self.check_beam_keys(usage="outputting to hdf5 file")

        if not complete_keys:  # don't even write an full_kick_output as the keys aren't complete (for some reason)
            raise RuntimeError("The keys could not be validated for the beam dict to full_kick_output to file")

        continue_flag = general_tools.control_file_overwrite_logic(filename, overwrite_existing_file)

        print("I have this continue flag", continue_flag)

        if not continue_flag:  # if this is false, our control logic doesn't want us to continue!
            return False

        if isinstance(zoffset, (list, np.ndarray)) and len(zoffset) == 3:
            xoffset = zoffset[0]
            yoffset = zoffset[1]
            zoffset = zoffset[2]

        with h5py.File(filename, "w") as f:
            inputgrp = f.create_group("Parameters")

            # shouldn't be needed because of my control at the generation steps
            # if not 'total_charge' in self._beam or self._beam['total_charge'] == 0:
            #    self._beam['total_charge'] = np.sum(self._beam['charge'])

            if sourcefilename is not None:
                inputgrp['Source'] = sourcefilename
            if pos is not None:
                inputgrp['Starting_Position'] = pos
            else:
                inputgrp['Starting_Position'] = [0, 0, 0]
            if rotation is not None:
                inputgrp['Rotation'] = rotation
            else:
                inputgrp['Rotation'] = 0

            # inputgrp['total_charge'] = self._beam['total_charge']
            inputgrp['total_charge'] = self.total_charge
            inputgrp['npart'] = len(self.x)
            inputgrp['centered'] = centered
            inputgrp['code'] = self.code
            inputgrp['particle_mass'] = mass
            beamgrp = f.create_group("beam")

            if 'reference_particle' in self._beam:
                beamgrp['reference_particle'] = self._beam['reference_particle']
            if 'status' in self._beam:
                beamgrp['status'] = self._beam['status']
            beamgrp['longitudinal_reference'] = longitudinal_reference

            # a default arg removed from args - could be critical for simframe (?) so leaving ina s full_kick_output
            beamgrp['cathode'] = False
            # print('hdf5 write cathode', cathode, np.array(beamgrp['cathode']))
            if len(self._beam['charge']) == len(self.x):
                chargevector = self._beam['charge']
            else:
                chargevector = np.full(len(self.x), self.charge / len(self.x))
            array = np.array(
                [self.x + xoffset, self.y + yoffset, self.z + zoffset, self.cpx, self.cpy, self.cpz, self.t + toffset,
                 chargevector]).transpose()
            beamgrp['columns'] = np.array(['x', 'y', 'z', 'cpx', 'cpy', 'cpz', 't', 'q'], dtype='S')
            beamgrp['units'] = np.array(['m', 'm', 'm', 'eV', 'eV', 'eV', 's', 'e'], dtype='S')
            beamgrp.create_dataset("beam", data=array)

        return os.path.isfile(filename)  # because we wrote a file!

    # region General methods
    # -----------------GENERAL METHODs------------------------

    def emittance(self, x, xp, p=None):
        '''

        :param x:
        :param xp:
        :param p:
        :return:
        '''
        emittance = np.sqrt(self.covariance(x, x) * self.covariance(xp, xp) - self.covariance(x, xp) ** 2)
        if p is None:
            return emittance
        else:
            gamma = np.mean(p) / self.E0_eV
            return gamma * emittance

    def covariance(self, u, up):
        u2 = u - np.mean(u)
        up2 = up - np.mean(up)
        return np.mean(u2 * up2) - np.mean(u2) * np.mean(up2)

    def _mom_to_angle(self, pu, pz):
        return np.arctan(pu / pz)

    def _angle_to_mom(self, u_prime, pz):
        return pz * np.tan(u_prime)

    def get_time_coord(self):
        # James Jones list comprehension
        assert hasattr(self, 'z')
        assert hasattr(self, 'Bz')
        t_temp = [(z / (-1 * Bz * constants.speed_of_light)) for z, Bz in zip(self.z, self.Bz)]
        return t_temp

    def setBeamCharge(self, Q):
        self._beam["charge"] = Q

    def returnBeam(self):
        # Like dumping a beam out to a fixed reference name

        return copy.deepcopy(self)

    def mve_emittance(self, x, xp, p=None):
        (center, radii, rotation, hullP) = MVE.getMinVolEllipse(list(zip(x, xp)), .01)
        emittance = radii[0] * radii[1]
        if p is None:
            return emittance
        else:
            gamma = np.mean(p) / self.E0_eV
            return gamma * emittance

    def eta_correlation(self, u):
        return self.covariance(u, self.p) / self.covariance(self.p, self.p)

    def eta_corrected(self, u):
        return u - self.eta_correlation(u) * self.p

    def get_2nd_moment_percentile(self,parameter,percentile = 95):
        """
        For more information see general_tools get_percentile_values
        :param parameter: string to get variable from, eg 'x'
        :param percentile: percentile to centre on
        :return: 2nd moment of that selected percentile
        """
        percentile_array = general_tools.get_percentile_values(values=getattr(self,parameter),percentile=percentile)

        return np.sqrt(self.covariance(percentile_array, percentile_array))



    # -----Optics functions---------
    def shiftBeamVertically(self, dy):
        self._beam['y'] = self.y + dy

    def shiftBeam(self, parameter: str, delta):
        """
        Shift a beam by a given amount in any of 6D space
        :param parameter: str for beam key to be shifted
        :param delta: number for amount to shift by
        :return: nothing
        """
        assert parameter in self.required_beam_keys, "must shift by a corrected labelled parameter of the beam"

        self._beam[parameter] += delta

    def driftBeam(self, L):
        self._beam['x'] = self.x + L * self.xp
        self._beam['y'] = self.y + L * self.yp
        self._beam['z'] = self.z + L
        # xp,yp,pz doesn't change

    def thinLensFQuadBeam(self, f):
        temp_xp = (-1 * self.x / f) + self.xp
        temp_yp = (self.y / f + self.yp)

        # convert back to our set variables of px, py
        self._beam['px'] = self._angle_to_mom(temp_xp, self.pz)
        self._beam['py'] = self._angle_to_mom(temp_yp, self.pz)

    def thickQuadBeam(self, K, L):
        """
        Thick quad, K>0 means focusing applied in X, defocusing applied in y
        :param K: Quad strength
        :param L: quad length
        :return: none, sets x,y,px,py,z based on quad
        """
        if K >= 0:
            new_x,new_xp = self.__transportMatrixThickFQuad(self.x,self.xp,K,L)
            new_y,new_yp = self.__transportMatrixThickDQuad(self.y,self.yp,K,L)
        else:
            K = -1*K
            new_x,new_xp = self.__transportMatrixThickDQuad(self.x,self.xp,K,L)
            new_y,new_yp = self.__transportMatrixThickFQuad(self.y,self.yp,K,L)

        #Apply to beam
        self._beam['x'] = new_x
        self._beam['px'] = self._angle_to_mom(new_xp, self.pz)

        self._beam['y'] = new_y
        self._beam['py'] = self._angle_to_mom(new_yp, self.pz)

        self._beam['z'] = self.z + L

    def thickWakeQuadBeam(self,dk_dz,k0,L):
        """
        Special sauce thick quad that varies linearly along bunch - using normalsied coord z0
        Only acts horizontally for now
        :param dk_dz: gradient of variation
        :param k0: intetercept of variation at z0
        :param L: Length of interation for special quad
        :return:
        """
        new_x, new_xp = self.__transportMatrixThickWakeFQuad(self.x, self.xp, self.zn, dk_dz,k0, L)
        #new_y, new_yp = self.__transportMatrixThickDQuad(self.y, self.yp, K, L)

        #apply back to beam
        self._beam['x'] = new_x
        self._beam['px'] = self._angle_to_mom(new_xp, self.pz)
        self._beam['z'] = self.z + L

    def __transportMatrixThickWakeFQuad(self,u,uprime,dz,dk_dz,k0,L):

        final_output = np.zeros((2,self.nMacros))
        i=0
        for u_part,u_prime_part,dz_part in zip(u,uprime,dz):
            K = (dz_part * dk_dz) + k0
            if K <0:
                K=0

            phi = np.sqrt(K) * L
            mat = np.array([[np.cos(phi), (1 / np.sqrt(K) * np.sin(phi))],
                            [(-1 * np.sqrt(K) * np.sin(phi)), np.cos(phi)]])

            output = np.matmul(mat, np.array([u_part, u_prime_part]))

            final_output[:,i] += output

            i +=1

        return final_output[0], final_output[1]

    def __transportMatrixThickFQuad(self,u,uprime,K,L):
        """
        Take canonical variables u, u' and transform through thick focusing quad of strength K, length L
        :param u: x or y
        :param uprime: x' or y'
        :param K: quad strength
        :param L: length
        :return: transformed u,u'
        """

        phi = np.sqrt(K) * L
        mat = np.array([[np.cos(phi),(1/np.sqrt(K)*np.sin(phi))],
                        [(-1*np.sqrt(K)*np.sin(phi)),np.cos(phi)]])

        output = np.matmul(mat,np.array([u,uprime]))

        return output[0],output[1]


    def __transportMatrixThickDQuad(self,u,uprime,K,L):
        """
        Take canonical variables u, u' and transform through thick DE-focusing quad of strength K, length L
        :param u: x or y
        :param uprime: x' or y'
        :param K: quad strength
        :param L: length
        :return: transformed u,u'
        """

        K = np.abs(K) # first clean up K if used as negative
        phi = np.sqrt(K) * L
        mat = np.array([[np.cosh(phi),(1/np.sqrt(K)*np.sinh(phi))],
                        [(1*np.sqrt(K)*np.sinh(phi)),np.cosh(phi)]])

        output = np.matmul(mat,np.array([u,uprime]))

        return output[0], output[1]

    # endregion

    # = -------------GENERAL Properties-----------------------

    @property
    def parameter_units(self):
        # This is explanatory, and can be called by other functions for full_kick_output, plotting etc
        dict = {
            'x': 'm',
            'y': 'm',
            'z': 'm',
            'zn': 'm',
            't': 's',
            'charge': 'C',
            'total_charge': 'C',
            'cpx': 'eV/c',
            'cpy': 'eV/c',
            'cpz': 'eV/c',
            'xp': 'rad',
            'yp': 'rad',
            'px': 'kgm/s',
            'py': 'kgm/s',
            'pz': 'kgm/s',

        }
        return dict

    # region required variables
    # --------------------Variables that must be defined----------------------
    @property
    def x(self):
        return self._beam['x']

    @property
    def y(self):
        return self._beam['y']
    @property
    def r(self):
        return np.sqrt(self._beam['x']**2 + self._beam['y']**2)
    @property
    def z(self):
        return self._beam['z']

    @property
    def px(self):
        return self._beam['px']

    @property
    def py(self):
        return self._beam['py']

    @property
    def pz(self):
        return self._beam['pz']

    @property
    def t(self):
        return self._beam['t']

    @property
    def charge(self):
        return self._beam['charge']  # array of charges

    @property
    def macro_type(self):
        if self._beam['macro_type'] == 'fixed_weight' \
                or self._beam['macro_type'] == 'variable_weight':
            return self._beam['macro_type']
        else:
            return 'undefined'

    @property
    def code(self):
        if 'code' in self._beam:
            return self._beam['code']
        else:
            return 'unknown'

    # endregion

    # region derived parameters

    # ---------------------------------------------------------------------------
    # = properties that are defined from others

    @property
    def charge_per_macro(self):
        if self.macro_type == 'fixed_weight':
            return self.charge[0]  # constant if fixed weight
        else:
            return self.charge  # array if variable weight

    @property
    def total_charge(self):
        return np.sum(self.charge)

    @property
    def cpx(self):
        return (self.px) / self.q_over_c

    @property
    def cpy(self):
        return (self.py) / self.q_over_c

    @property
    def cpz(self):
        return (self.pz) / self.q_over_c

    @property
    def xp(self):
        return self._mom_to_angle(self.px, self.pz)

    @property
    def yp(self):
        return self._mom_to_angle(self.py, self.pz)

    @property
    def zn(self):
        return self.z - self.z0

    @property
    def z0(self):
        return np.mean(self.z)

    @property
    def Sx(self):
        return np.sqrt(self.covariance(self.x, self.x))

    @property
    def Sy(self):
        return np.sqrt(self.covariance(self.y, self.y))

    @property
    def Sz(self):
        return np.sqrt(self.covariance(self.z, self.z))

    @property
    def Sx_95(self):
        return self.get_2nd_moment_percentile('x',95)

    @property
    def Sy_95(self):
        return self.get_2nd_moment_percentile('y',95)

    @property
    def Sz_95(self):
        return self.get_2nd_moment_percentile('z',95)

    @property
    def Mx(self):
        return np.mean(self._beam['x'])

    @property
    def My(self):
        return np.mean(self._beam['y'])

    @property
    def Mz(self):
        return np.mean(self._beam['z'])

    @property
    def Mzn(self):
        return 0.0

    @property
    def cp(self):
        return np.sqrt(self.cpx ** 2 + self.cpy ** 2 + self.cpz ** 2)

    @property
    def p(self):
        return self.cp * self.q_over_c

    @property
    def Brho(self):
        return np.mean(self.p) / constants.elementary_charge

    @property
    def gamma(self):
        return np.sqrt(1 + (self.cp / self.E0_eV) ** 2)

    @property
    def BetaGamma(self):
        return self.cp / self.E0_eV

    @property
    def vx(self):
        velocity_conversion = 1 / (constants.m_e * self.gamma)
        return velocity_conversion * self.px

    @property
    def vy(self):
        velocity_conversion = 1 / (constants.m_e * self.gamma)
        return velocity_conversion * self.py

    @property
    def vz(self):
        velocity_conversion = 1 / (constants.m_e * self.gamma)
        return velocity_conversion * self.pz

    @property
    def Bx(self):
        return self.vx / constants.speed_of_light

    @property
    def By(self):
        return self.vy / constants.speed_of_light

    @property
    def Bz(self):
        return self.vz / constants.speed_of_light

    @property
    def nMacros(self):
        return len(self.x)

    # endregion

    # region emittances

    @property
    def normalized_horizontal_emittance(self):
        return self.emittance(self.x, self.xp, self.cp)

    @property
    def normalized_vertical_emittance(self):
        return self.emittance(self.y, self.yp, self.cp)

    @property
    def horizontal_emittance(self):
        return self.emittance(self.x, self.xp)

    @property
    def vertical_emittance(self):
        return self.emittance(self.y, self.yp)

    @property
    def normalized_mve_horizontal_emittance(self):
        return self.mve_emittance(self.x, self.xp, self.cp)

    @property
    def normalized_mve_vertical_emittance(self):
        return self.mve_emittance(self.y, self.yp, self.cp)

    @property
    def horizontal_mve_emittance(self):
        return self.mve_emittance(self.x, self.xp)

    @property
    def vertical_mve_emittance(self):
        return self.mve_emittance(self.y, self.yp)

    @property
    def horizontal_emittance_90(self):
        emit = self.horizontal_emittance
        alpha = self.alpha_x
        beta = self.beta_x
        gamma = self.gamma_x
        emiti = gamma * self.x ** 2 + 2 * alpha * self.x * self.xp + beta * self.xp * self.xp
        # for each particle, sort them in amplitude, and then take the one at 90% in the sorted list
        return sorted(emiti)[int(0.9 * len(emiti) - 0.5)]

    @property
    def normalized_horizontal_emittance_90(self):
        emit = self.horizontal_emittance_90
        return np.mean(self.cp) / self.E0_eV * emit

    @property
    def vertical_emittance_90(self):
        emit = self.vertical_emittance
        alpha = self.alpha_y
        beta = self.beta_y
        gamma = self.gamma_y
        emiti = gamma * self.y ** 2 + 2 * alpha * self.y * self.yp + beta * self.yp * self.yp
        # for each particle, sort them in amplitude, and then take the one at 90 % in the sorted list
        return sorted(emiti)[int(0.9 * len(emiti) - 0.5)]

    @property
    def normalized_vertical_emittance_90(self):
        emit = self.vertical_emittance_90
        return np.mean(self.cp) / self.E0_eV * emit

    # endregion

    # region twiss parameters

    # TODO add some tasty print function with lines and hashes

    @property
    def beta_x(self):
        return self.covariance(self.x, self.x) / self.horizontal_emittance

    @property
    def alpha_x(self):
        return -1 * self.covariance(self.x, self.xp) / self.horizontal_emittance

    @property
    def gamma_x(self):
        return self.covariance(self.xp, self.xp) / self.horizontal_emittance

    @property
    def beta_y(self):
        return self.covariance(self.y, self.y) / self.vertical_emittance

    @property
    def alpha_y(self):
        return -1 * self.covariance(self.y, self.yp) / self.vertical_emittance

    @property
    def gamma_y(self):
        return self.covariance(self.yp, self.yp) / self.vertical_emittance

    @property
    def twiss_analysis(self):
        return self.horizontal_emittance, self.alpha_x, self.beta_x, self.gamma_x, self.vertical_emittance, self.alpha_y, self.beta_y, self.gamma_y

    # endregion

    # TODO print functions here need formatting
    def print_all_beam_properties_table(self):
        print("-----Beam properties----- ")
        print("Beam generated from code: " + self.code)
        print(f"Number of macro particles = {self.nMacros} ")

        self.print_charge_properties()

        print("-------------------------")
        print("Longitudinal properties")

    def print_charge_properties(self):
        print("Macro particle type: " + self.macro_type)
        print(f"Total charge in beam: {self.total_charge} C")
        print(f"Charge per macro : {self.charge_per_macro} C")


class GaussianBeam(GeneralBeam):
    """
    A class to generate, drift,focus a beam of fixed weight macros
    should write a beam dict to the GeneralBeam class spec...TBC
    """

    def __init__(self):
        super().__init__()
        self.resetDicts()

        self.acceptedLongitudinalProfiles = ["Gaussian", "Uniform", "Plateau", "DoubleGauss", "SkewGaussian"]
        self.requiredParameters_waistGenerator = ["pz_MeV", "eps_x_N", "eps_y_N", "sig_x_0", "sig_y_0", "sig_z_0",
                                                  "sig_pz_0", "lCorr_Fac",
                                                  "x_0", "y_0", "xp_0", "yp_0", "z_0", "charge_per_macro",
                                                  "LongitudinalProfile"]

    def setBeamParameters(self, bp):

        print("Checking parameters for beam")
        self.checkBeamParameters(bp, self.requiredParameters_waistGenerator)
        self.testLongitudinalProfile(bp)
        print("Setting parameters for beam")
        self._beam = bp

    def checkBeamParameters(self, bp, requiredKeys):

        general_tools.check_dict_params(bp, requiredKeys, "beam parameters", "beam generator")

    def generateTransverseMacros(self, numMacros):
        x, xp, y, yp = np.random.default_rng().multivariate_normal(self.meanArray_waist, self.covarianceMatrix_waist,
                                                                   numMacros).T
        return (x, xp, y, yp)

    def testLongitudinalProfile(self, bp):

        print("Testing longitudinal profile")

        try:
            bp["LongitudinalProfile"]
        except KeyError:
            print("Longitudinal profile is not defined")
            raise

        if bp["LongitudinalProfile"] not in self.acceptedLongitudinalProfiles:
            print("The set longitudinal profile of ", bp["LongitudinalProfile"], " is not accepted")
            print("Accepted profiles are", self.acceptedLongitudinalProfiles)
            raise ValueError

        # Now to test profiles individually to make sure they have required values

        if bp["LongitudinalProfile"] == "Gaussian":
            print("Checking parameters for Gaussian profile")
            extra_keys = []
            self.checkBeamParameters(bp, extra_keys)


        elif bp["LongitudinalProfile"] == "Plateau":
            print("Checking parameters for Plateau profile")
            extra_keys = ["plat_rise"]
            self.checkBeamParameters(bp, extra_keys)

        elif bp["LongitudinalProfile"] == "DoubleGauss":
            print("Checking parameters for DoubleGauss profile")
            extra_keys = ["sig_z_2", "offset", "rel_amp"]
            self.checkBeamParameters(bp, extra_keys)
            # extra check is needed
            if bp["offset"] < 0:
                print("Offset must be a positive, relative amplitude can be altered to get desired shape")
                raise ValueError

        elif bp["LongitudinalProfile"] == "SkewGaussian":
            print("Checking parameters for Skew-Gaussian profile")
            extra_keys = []
            self.checkBeamParameters(bp, extra_keys)
        elif bp["LongitudinalProfile"] == "Uniform":
            print("Checking parameters for Uniform profile")
            extra_keys = []
            self.checkBeamParameters(bp, extra_keys)
        else:
            print("Unexpected error when setting longitudinal profile, accepted profiles are")
            print(self.acceptedLongitudinalProfiles)
            raise RuntimeError

    def generateLongitudinalMacros(self, numMacros):

        if self._beam['LongitudinalProfile'] == "Gaussian":
            z = np.random.default_rng().normal(self._beam['z_0'], self._beam['sig_z_0'], numMacros)
            return z

        elif self._beam['LongitudinalProfile'] == 'SkewGaussian':
            z = skewnorm.rvs(loc=self._beam['z_0'], scale=self._beam['sig_z_0'], a=self._beam['skew'], size=numMacros)
            return z

        elif self._beam['LongitudinalProfile'] == "Uniform":
            z = np.random.default_rng().uniform(-2 * self._beam["sig_z_0"], 2 * self._beam["sig_z_0"], numMacros)
            return z

        elif self._beam['LongitudinalProfile'] == "Plateau":
            z = self.getPlateauValues(self._beam["sig_z_0"], self._beam["plat_rise"], numMacros)
            return z

        elif self._beam['LongitudinalProfile'] == "DoubleGauss":
            z = self.getDoubleGaussValues(self._beam["sig_z_0"], self._beam["sig_z_2"], self._beam["offset"],
                                          self._beam["rel_amp"], numMacros)
            return z

        else:
            print("Something has gone wrong and longitudinal macros not generated")
            raise RuntimeError

    def generateMomentumMacros(self):

        try:
            self._beam['z']
        except KeyError:
            print("Longitudinal macros not yet defined!")
            raise

        # TODO tidy this syntax up
        # Empty array to start
        print("Generating momentum macros")
        # cpz = []

        uncorr = np.random.default_rng().normal(loc=0, scale=self._beam["sig_pz_0"], size=self.zn.shape)
        corr = (self._beam["lCorr_Fac"] * self.zn) / (4 * self.Sz)

        cpz = (corr + uncorr) + self._beam['pz_MeV']
        # cpz.append(cpzMacro)

        # Convert to numpy array
        # cpz = np.asarray(cpz)
        return (cpz)

    def generate6DMacros(self, numMacros):

        print("Generating macros")

        self._beam['z'] = self.generateLongitudinalMacros(numMacros)

        # important conversion
        self._beam['pz'] = self.generateMomentumMacros() * self.q_over_c * constants.mega

        x, xp, y, yp = self.generateTransverseMacros(numMacros)

        self._beam['x'] = x
        self._beam['y'] = y
        self._beam['px'] = self._angle_to_mom(xp, self.pz)
        self._beam['py'] = self._angle_to_mom(yp, self.pz)

        # TODO replace with general beam method
        # James Jones list comprehension
        # t_temp = [(z / (-1 * Bz * constants.speed_of_light)) for z, Bz in zip(self.z, self.Bz)]

        t_temp = self.get_time_coord()

        # convert it to our lovely jubbly np arrays
        self._beam['t'] = np.asarray(t_temp)

        self._beam['code'] = "GaussianGenerator"

        self._beam['charge'] = np.full_like(self._beam['x'], self._beam['charge_per_macro'])
        self._beam['macro_type'] = 'fixed_weight'

        self.check_beam_keys(usage="generating Gaussian Beam")

    def calcWaistOpticFunctions(self):

        self._beam["alpha_x"] = 0
        self._beam["alpha_y"] = 0

        self._beam["gamma"] = np.sqrt(1 + (self._beam["pz_MeV"] / self.E0_MeV) ** 2)

        self._beam["eps_x_g"] = self._beam["eps_x_N"] / self._beam["gamma"]
        self._beam["eps_y_g"] = self._beam["eps_y_N"] / self._beam["gamma"]

        self._beam["beta_x_0"] = (self._beam["sig_x_0"]) ** 2 / self._beam["eps_x_g"]
        self._beam["beta_y_0"] = (self._beam["sig_y_0"]) ** 2 / self._beam["eps_y_g"]

        self._beam["gamma_x_0"] = 1 / self._beam["beta_x_0"]
        self._beam["gamma_y_0"] = 1 / self._beam["beta_y_0"]

    def getSampleFromArbDist(self, pdf, num_samples, max_z, nMacros):
        dist = ADS.Distribution(pdf, transform=lambda x: ((x / num_samples) - 0.5) * (max_z / 0.5))
        return dist(nMacros)

    def getPlateauValues(self, sig_z, plat_rise, nMacros, num_samples=5000):
        min_z = -2 * sig_z - 3 * plat_rise
        max_z = 2 * sig_z + 3 * plat_rise

        x = np.linspace(min_z, max_z, num_samples)

        pdf = np.piecewise(x,
                           [x < -2 * sig_z, x > 2 * sig_z],
                           [lambda x: np.exp(-1 / 2 * ((x + 2 * sig_z) ** 2) / plat_rise ** 2),
                            lambda x: np.exp(-1 / 2 * ((x - 2 * sig_z) ** 2) / plat_rise ** 2),
                            1
                            ]
                           )

        data = self.getSampleFromArbDist(pdf, num_samples, max_z, nMacros)

        return np.asarray(*data)

    def getDoubleGaussValues(self, sig_z, sig_z_2, offset, rel_amp, nMacros, num_samples=10000):

        rel_amp = abs(rel_amp)  # just in case a negative is passed

        min_z = min(-3 * sig_z, offset - (3 * sig_z_2))
        max_z = max(3 * sig_z, offset + (3 * sig_z_2))

        print(min_z, max_z)

        x = np.linspace(min_z, max_z, num_samples)

        pdf = np.exp((-1 / 2) * ((x ** 2) / (sig_z ** 2))) + (
                rel_amp * np.exp((-1 / 2) * (((x - offset) ** 2) / (sig_z_2 ** 2))))

        # pdf = np.exp( (-1/2) * ((x**2)/(sig_z**2)) ) * rel_amp

        data = self.getSampleFromArbDist(pdf, num_samples, max_z, nMacros)

        return np.asarray(*data)

    # -------Properties----------

    @property
    def meanArray_waist(self):
        return [self._beam["x_0"], self._beam["xp_0"],
                self._beam["y_0"], self._beam["yp_0"]]

    @property
    def covarianceMatrix_waist(self):  # TODO change this to not take dicts and use beta_0 but just use betas
        return [
            [self._beam["eps_x_g"] * self._beam["beta_x_0"], -1 * self._beam["eps_x_g"] * self._beam["alpha_x"], 0, 0],
            [-1 * self._beam["eps_x_g"] * self._beam["alpha_x"], self._beam["eps_x_g"] * self._beam["gamma_x_0"], 0, 0],
            [0, 0, self._beam["eps_y_g"] * self._beam["beta_y_0"], -1 * self._beam["eps_y_g"] * self._beam["alpha_y"]],
            [0, 0, -1 * self._beam["eps_y_g"] * self._beam["alpha_y"], self._beam["eps_y_g"] * self._beam["gamma_y_0"]]]


# -----------------------------------------------------------------


class BeamFromFile(GeneralBeam):
    def __init__(self):
        super().__init__()
        self.resetDicts()

    def read_HDF5_beam_file(self, filename, local=False):
        self.resetDicts()
        with h5py.File(filename, "r") as h5file:
            if h5file.get('beam/reference_particle') is not None:
                self._beam['reference_particle'] = np.array(h5file.get('beam/reference_particle'))
            if h5file.get('beam/longitudinal_reference') is not None:
                self._beam['longitudinal_reference'] = np.array(h5file.get('beam/longitudinal_reference'))
            else:
                self._beam['longitudinal_reference'] = 't'
            if h5file.get('beam/status') is not None:
                self._beam['status'] = np.array(h5file.get('beam/status'))
            x, y, z, cpx, cpy, cpz, t, charge = np.array(h5file.get('beam/beam')).transpose()
            cp = np.sqrt(cpx ** 2 + cpy ** 2 + cpz ** 2)
            self._beam['x'] = x
            self._beam['y'] = y
            self._beam['z'] = z
            # self.beam['cpx'] = cpx
            # self.beam['cpy'] = cpy
            # self.beam['cpz'] = cpz
            self._beam['px'] = cpx * self.q_over_c
            self._beam['py'] = cpy * self.q_over_c
            self._beam['pz'] = cpz * self.q_over_c
            # self.beam['cp'] = cp
            # self.beam['p'] = cp * self.q_over_c
            # self.beam['xp'] = np.arctan(self.px/self.pz)
            # self.beam['yp'] = np.arctan(self.py/self.pz)
            self._beam['clock'] = np.full(len(self.x), 0)
            # self.beam['gamma'] = np.sqrt(1+(self.cp/self.E0_eV)**2)
            # velocity_conversion = 1 / (constants.m_e * self.gamma)
            # self.beam['vx'] = velocity_conversion * self.px
            # self.beam['vy'] = velocity_conversion * self.py
            # self.beam['vz'] = velocity_conversion * self.pz
            # self.beam['Bx'] = self.vx / constants.speed_of_light
            # self.beam['By'] = self.vy / constants.speed_of_light
            # self.beam['Bz'] = self.vz / constants.speed_of_light
            self._beam['t'] = t
            self._beam['charge'] = charge

            self._beam['macro_type'] = 'fixed_weight'

            startposition = np.array(h5file.get('/Parameters/Starting_Position'))
            startposition = startposition if startposition is not None else [0, 0, 0]
            self._beam['starting_position'] = startposition
            theta = np.array(h5file.get('/Parameters/Rotation'))
            theta = theta if theta is not None else 0
            self._beam['rotation'] = theta
            code = np.array(h5file.get('/Parameters/code'))  # Added for completeness as its in the write method
            self._beam['code'] = code.item().decode("utf-8")  # wee bit hacky to get it to behave nicely
            if local == True:
                self.rotate_beamXZ(self._beam['rotation'], preOffset=self._beam['starting_position'])

            self.check_beam_keys(usage="importing hdf5 beam file")

    def read_WakeCode_beam_file(self, filename, remove_losses=True, force_charge=None):
        """

        :param filename:
        :param remove_losses:
        :param force_charge:
        :return:
        """
        # Each macro with 0 charge has been lost
        # Import all macros
        # assess all ones with charge = 0
        # delete from ALL dict arrays indexes with charge = 0
        self.resetDicts()
        with h5py.File(filename, "r") as h5file:

            x, y, z, cpx, cpy, cpz, t, charge = np.array(h5file.get('beam/beam')).transpose()

            # TODO this is a hack - ask toby to fix his files
            if remove_losses:
                print("Checking for lost particles")
                if np.any(charge == 0):
                    print("Beam losses detected")
                    macro_select = np.nonzero(charge)
                else:
                    print("No losses detected in file read in")
                    macro_select = ...
            else:
                print("Removing lost particles not specified, keeping all macros in file ")
                macro_select = ...

            self._beam['x'] = x[macro_select]
            self._beam['y'] = y[macro_select]
            self._beam['z'] = z[macro_select]
            # self.beam['cpx'] = cpx
            # self.beam['cpy'] = cpy
            # self.beam['cpz'] = cpz
            self._beam['px'] = cpx[macro_select] * self.q_over_c
            self._beam['py'] = cpy[macro_select] * self.q_over_c
            self._beam['pz'] = cpz[macro_select] * self.q_over_c
            # self.beam['cp'] = cp
            # self.beam['p'] = cp * self.q_over_c
            # self.beam['xp'] = np.arctan(self.px/self.pz)
            # self.beam['yp'] = np.arctan(self.py/self.pz)
            # self.beam['clock'] = np.full(len(self.x), 0)

            # convert it to our lovely jubbly np arrays
            self._beam['t'] = t[macro_select]

            self._beam['charge'] = charge[macro_select]
            # self.beam['total_charge'] = total_charge

            # self.beam['charge_per_macro'] = self.beam['total_charge'] / len(self.x)

            self._beam['code'] = "WakeCodeOutput"
            self._beam['macro_type'] = 'fixed_weight'

            self.check_beam_keys(usage="importing WakeCode beam file")

    def read_beam_from_fbPIC_file(self, filename, dump_number=0):
        self.resetDicts()

        with h5py.File(filename, "r") as h5file:
            data_subsets = []
            for name in h5file['data']:
                data_subsets.append(name)

            data_value = data_subsets[dump_number]

            # print('data/' + data_value + '/particles/beam/weighting' in h5file)
            weights = np.array(h5file.get('data/' + data_value + '/particles/beam/weighting'))

            fixed_weight = np.all(weights == weights[0])

            if fixed_weight:
                print("These are fixed weight macros")
                self._beam['macro_type'] = 'fixed_weight'
                self._beam['code'] = "fbPIC"
            else:
                raise ValueError("Not all macors have same weight, variable weight macros not supported")

            # print('data/' + data_value + '/particles/beam/position' in h5file)
            positions = {'x': np.empty(weights.shape),
                         'y': np.empty(weights.shape),
                         'z': np.empty(weights.shape)}

            # these are for the dimensions in the data set - labelled in x,y,z
            momenta = {'x': np.empty(weights.shape),
                       'y': np.empty(weights.shape),
                       'z': np.empty(weights.shape)}

            for key in positions:
                data_loc = 'data/' + data_value + '/particles/beam/position/' + key
                if data_loc in h5file:
                    SI_convert_factor = h5file.get(data_loc).attrs.get('unitSI')
                    positions[key] = np.array(h5file.get(data_loc)) * SI_convert_factor

            for key in momenta:
                data_loc = 'data/' + data_value + '/particles/beam/momentum/' + key
                if data_loc in h5file:
                    print(h5file.get(data_loc).attrs.keys())
                    SI_convert_factor = h5file.get(data_loc).attrs.get('unitSI')
                    momenta[key] = np.array(h5file.get(data_loc)) * SI_convert_factor

            # set the values from the file to beam dict
            self._beam['x'] = positions['x']
            self._beam['y'] = positions['y']
            self._beam['z'] = positions['z']

            self._beam['px'] = momenta['x']
            self._beam['py'] = momenta['y']
            self._beam['pz'] = momenta['z']

            self._beam['charge'] = weights * constants.elementary_charge

            # James Jones list comprehension
            # t_temp = [(z / (-1 * Bz * constants.speed_of_light)) for z, Bz in zip(self.z, self.Bz)]

            t_temp = self.get_time_coord()

            # convert it to our lovely jubbly np arrays
            self._beam['t'] = np.asarray(t_temp)

        self.check_beam_keys(usage="import fbPIC file")


class beamBinner(object):
    """
    Takes a beam object from either gaussian generator bunch or from Simframe read_beam_file (to be tested)
    """

    def __init__(self, beam):
        self.beam = beam
        self.mesh = {}
        self.smoothedMesh = {}

    def clearDicts(self):
        self.mesh = {}
        self.smoothedMesh = {}

    def binTheBeam(self, Lx, Dx, Ly, Dy, Lz, Dz):
        """
        Method for binning beam. Works on Li,Di only, if Li not set correctly could clip beam charge significantly!

        :param Lx: Length around mean for bins in x
        :param Dx: Separation between bins in x
        :param Ly: Length around mean for bins in y
        :param Dy: Separation between bins in y
        :param Lz: Length around mean for bins in z
        :param Dz: Separation between bins in z
        :return: Fraction of charge captured in mesh
        """
        # Setup a mesh grid in each dim

        xMin = self.beam.Mx - (Lx / 2)
        xMax = self.beam.Mx + (Lx / 2)
        nBinsx = int(Lx / Dx)

        yMin = self.beam.My - (Ly / 2)
        yMax = self.beam.My + (Ly / 2)
        nBinsy = int(Ly / Dy)

        zMin = self.beam.Mzn - (Lz / 2)
        zMax = self.beam.Mzn + (Lz / 2)
        nBinsz = int(Lz / Dz)

        # Sanity check mesh parameters
        self.mesh['cps_x'] = self.beam.Sx / Dx
        self.mesh['cps_y'] = self.beam.Sy / Dy
        self.mesh['cps_z'] = self.beam.Sz / Dz

        self.mesh['xRange'] = (xMin, xMax)
        self.mesh['yRange'] = (yMin, yMax)
        self.mesh['zRange'] = (zMin, zMax)

        # Make 3D hist to bin the beam
        # Must use beam.zn here!
        H, edges = np.histogramdd((self.beam.x, self.beam.y, self.beam.zn),
                                  bins=(nBinsx, nBinsy, nBinsz),
                                  range=((xMin, xMax), (yMin, yMax), (zMin, zMax)))

        self.mesh['beamHist'] = H
        self.mesh['histEdges'] = edges

        # Calculate the 'central' mesh points from the bin edges
        self.mesh['cenValsX'] = np.empty(len(edges[0][:-1]))  # empty array
        self.mesh['dx'] = edges[0][1] - edges[0][0]
        for i, val in enumerate(edges[0][:-1]):
            self.mesh['cenValsX'][i] = val + (self.mesh['dx'] / 2)

        self.mesh['cenValsY'] = np.empty(len(edges[1][:-1]))  # empty array
        self.mesh['dy'] = edges[1][1] - edges[1][0]
        for i, val in enumerate(edges[1][:-1]):
            self.mesh['cenValsY'][i] = val + (self.mesh['dy'] / 2)

        self.mesh['cenValsZ'] = np.empty(len(edges[2][:-1]))  # empty array
        self.mesh['dz'] = edges[2][1] - edges[2][0]
        for i, val in enumerate(edges[2][:-1]):
            self.mesh['cenValsZ'][i] = val + (self.mesh['dz'] / 2)

        self.mesh['nMacroHist'] = np.sum(self.mesh['beamHist'], axis=(0, 1, 2))

        return self.chargeFracOnMesh

    # -------Smoothing methods---------

    def smoothGaussian(self, **kwargs):

        self.smoothedMesh["smoothed"] = True

        self.smoothedMesh["smoothedHist"] = gaussian_filter(self.beamHist, **kwargs)
        self.smoothedMesh["smoothMethod"] = "Gaussian Smoothing"

        self.smoothedMesh["kwargs"] = str(kwargs)
        self.smoothedMesh["normFactor"] = np.sum(self.smoothedMesh["smoothedHist"], axis=(0, 1, 2))

    def smoothMedian(self, **kwargs):

        self.smoothedMesh["smoothed"] = True

        self.smoothedMesh["smoothedHist"] = median_filter(self.beamHist, **kwargs)
        self.smoothedMesh["smoothMethod"] = "Median Smoothing"

        self.smoothedMesh["kwargs"] = str(kwargs)
        self.smoothedMesh["normFactor"] = np.sum(self.smoothedMesh["smoothedHist"], axis=(0, 1, 2))

    def smoothUniform(self, **kwargs):

        self.smoothedMesh["smoothed"] = True

        self.smoothedMesh["smoothedHist"] = uniform_filter(self.beamHist, **kwargs)
        self.smoothedMesh["smoothMethod"] = "Uniform Smoothing"

        self.smoothedMesh["kwargs"] = str(kwargs)
        self.smoothedMesh["normFactor"] = np.sum(self.smoothedMesh["smoothedHist"], axis=(0, 1, 2))

    # --------Output-------------

    def write_hd5f_mesh_file(self, filename, outputMacros=False, overwrite_existing_file=False, includeSmoothed=True,
                             includeUnSmoothed=False):
        """
        Outputs the mesh to an HDF5 file

        :param filename: the filename to full_kick_output to
        :param outputMacros: control for saving addition dataset of the raw macros, default False
        :param includeSmoothed: include the smoothed dataset if created, default True
        :param includeUnSmoothed: include raw unsmoothed data in addition to smoothed, default False
        :return: Flag if file written
        """

        # First check if file exists

        # TODO replace with general tools version

        continue_flag = general_tools.control_file_overwrite_logic(filename, overwrite_existing_file)

        if not continue_flag:  # if this is false, our control logic doesn't want us to continue!
            return False

        with h5py.File(filename, "w") as f:

            inputgrp = f.create_group("MeshParameters")
            inputgrp['code'] = self.beam.code
            inputgrp['chargePerMacro'] = self.beam.charge_per_macro
            inputgrp['macrosInMesh'] = self.nMacroHist
            inputgrp['macroFractionOnMesh'] = self.nMacroHist / self.beam.nMacros

            inputgrp['xRange'] = self.xRange
            inputgrp['yRange'] = self.yRange
            inputgrp['zRange'] = self.zRange

            inputgrp['dx'] = self.dx
            inputgrp['dy'] = self.dy
            inputgrp['dz'] = self.dz

            meshgrp = f.create_group("MeshValues")
            meshgrp['xPoints'] = self.cenValsX
            meshgrp['yPoints'] = self.cenValsY
            meshgrp['zPoints'] = self.cenValsZ

            if includeSmoothed and self.smoothed:  # make sure you want to include the smoothed mesh (it will be default if available)

                print("Writing smoothed data to file")
                # Deal with smoothing if used
                inputgrp["Smoothed"] = self.smoothed
                meshgrp.create_dataset("macrosAtPoint",
                                       data=self.smoothedHist)  # if you've smoothed it and asked for it, then

                inputgrp['smoothingMethod'] = self.smoothMethod
                inputgrp['smooth_kwargs'] = self.smooth_kwargs

                if includeUnSmoothed:  # add this in if you want to as a backup
                    print("Writing unsmoothed data to file in addition to smoothed")
                    meshgrp.create_dataset("unSmoothed_macrosAtPoint", data=self.beamHist)


            else:  # just do the default mesh
                print("Writing unsmoothed data to file")
                meshgrp.create_dataset("macrosAtPoint", data=self.beamHist)

            # Optional pass out all the macros to also be implemented
            if outputMacros:
                beamgrp = f.create_group("beam")
                chargevector = np.full(self.beam.nMacros, self.beam.charge_per_macro)
                array = np.array(
                    [self.beam.x, self.beam.y, self.beam.z, self.beam.cpx, self.beam.cpy, self.beam.cpz, self.beam.t,
                     chargevector]).transpose()
                beamgrp['columns'] = np.array(['x', 'y', 'z', 'cpx', 'cpy', 'cpz', 't', 'q'], dtype='S')
                beamgrp['units'] = np.array(['m', 'm', 'm', 'eV', 'eV', 'eV', 's', 'e'], dtype='S')
                beamgrp.create_dataset("beam", data=array)

        return True

    # ------Properties---------

    @property
    def dx(self):
        return self.mesh['dx']

    @property
    def dy(self):
        return self.mesh['dy']

    @property
    def dz(self):
        return self.mesh['dz']

    @property
    def xRange(self):
        return self.mesh['xRange']

    @property
    def yRange(self):
        return self.mesh['yRange']

    @property
    def zRange(self):
        return self.mesh['zRange']

    @property
    def cenValsX(self):
        return self.mesh['cenValsX']

    @property
    def cenValsY(self):
        return self.mesh['cenValsY']

    @property
    def cenValsZ(self):
        return self.mesh['cenValsZ']

    @property
    def histEdges(self):
        return self.mesh['histEdges']

    @property
    def beamHist(self):
        return self.mesh['beamHist']

    @property
    def nMacroHist(self):
        return self.mesh['nMacroHist']

    @property
    def chargeFracOnMesh(self):
        return self.nMacroHist / self.beam.nMacros

    @property
    def smoothed(self):
        try:
            return self.smoothedMesh["smoothed"]

        except KeyError:
            return False

        finally:
            return self.smoothedMesh["smoothed"]

    @property
    def smoothedHist(self):
        return self.smoothedMesh["smoothedHist"] * (self.nMacroHist / self.smoothedMesh["normFactor"])

    @property
    def smoothMethod(self):
        return self.smoothedMesh["smoothMethod"]

    @property
    def smooth_kwargs(self):
        return self.smoothedMesh["kwargs"]


# --------------------------------------------------------


# ------------------------------------------------------------------------------------------------------------------


# -------------------------------------------------------

# Check this with a linear kick with one of my files thats used for beam binning

class BeamKicker(object):

    def __init__(self, beam):
        self.inputBeam = beam  # type:GeneralBeam  # a complete instance of the gaussianBunch object
        self.outputBeam = self.inputBeam.returnBeam()
        self.kicker = {}

        self.requiredParameters = ['freq', 'amplitude', 'length', 'phase', 'transverse_profile']
        self.acceptedTransverseProfiles = ['Uniform', 'Gaussian']

    def resetOutputBeam(self):
        self.outputBeam.resetDicts()

    def checkParameters(self, kp, requiredKeys):

        missedKeys = []
        for key in requiredKeys:
            if key not in kp:
                missedKeys.append(key)

        if len(missedKeys) != 0:
            print("Missing the following keys from kicker parameters dict")
            print(missedKeys)
            print("kicker cannot be generated")
            raise KeyError

        print("Check passed")

    def checkTransverseProfile(self, kp):

        print("Testing transverse kicker profile")

        try:
            kp["transverse_profile"]
        except KeyError:
            print("transverse profile is not defined")
            raise

        if kp["transverse_profile"] not in self.acceptedTransverseProfiles:
            print("The set transverse profile of ", kp["transverse_profile"], " is not accepted")
            print("Accepted profiles are", self.acceptedTransverseProfiles)
            raise ValueError

        # Now to test profiles individually to make sure they have required values

        if kp["transverse_profile"] == "Uniform":
            print("Uniform transverse kicker profile, checking parameters")
            extra_keys = []
            self.checkParameters(kp, extra_keys)

        if kp["transverse_profile"] == "Gaussian":
            print("Gaussian transverse kicker profile, checking parameters")
            extra_keys = ['sig_r']
            self.checkParameters(kp, extra_keys)

    def setKickParameters(self, kp):

        print("Checking parameters for kicker")
        self.checkParameters(kp, self.requiredParameters)
        self.checkTransverseProfile(kp)
        print("Setting kicker parameters")
        self.kicker = kp

    def kickFunction(self, z_pos, x_pos, y_pos):

        # TODO switch here on type of kicker longitinally
        # Can have cosinusoidal
        # Maybe zn convolved
        #

        return self.amp_z * np.cos((self.k_z * z_pos) + self.phi) * self.l_z * self.transverseFactor(x_pos, y_pos)

    def transverseFactor(self, x_pos, y_pos):

        if self.kicker['transverse_profile'] == 'Uniform':
            return 1

        if self.kicker['transverse_profile'] == 'Gaussian':
            return np.exp((-1 / 2) * ((x_pos) / (self.kicker['sig_r'])) ** 2)

        else:
            print('Something has gone wrong with transverse factor and kick cannot be applied')
            raise RuntimeError

    def kickOutput(self):

        print("Applying kick")

        new_cpz = []
        for i, pz in enumerate(self.outputBeam._beam['pz']):
            cpz = pz / self.outputBeam.q_over_c
            delta_cpz = self.kickFunction(z_pos=self.outputBeam._beam['z'][i],
                                          x_pos=self.outputBeam._beam['x'][i],
                                          y_pos=self.outputBeam._beam['y'][i])

            new_cpz.append(cpz + delta_cpz)

        new_cpz = np.asarray(new_cpz)
        new_pz = new_cpz * self.outputBeam.q_over_c
        self.outputBeam._beam['pz'] = None
        self.outputBeam._beam['pz'] = np.asarray(new_pz)

    # Should not be used
    def returnOutputBeam(self):
        # TODO Test this
        return self.outputBeam.returnBeam()

    # -------------Properties------------------------

    @property
    def amp_z(self):
        return self.kicker['amplitude']

    @property
    def l_z(self):
        return self.kicker['length']

    @property
    def k_z(self):
        return (2 * constants.pi * self.kicker['freq']) / (constants.speed_of_light)

    @property
    def phi(self):
        return self.kicker['phase']

    @property
    def peakVoltage(self):
        return self.amp_z * self.l_z


# --------------------------------------------------------------------

class BeamPlotter(object):
    '''
    Class to do plots of beam objects
    '''

    # TODO
    # Add generalised KDe method for single beams - but make it any val against any val with nice labels!

    def __init__(self):
        self.beams = {}  # dict of beams to plot

    @property
    def variable_formats(self):
        # using latex style as we are plotting in matplotlib
        formats = {
            'xp': "x\'",
            'yp': "y\'",
            'zn': "$\zeta$",
            'cpz': "p$_z$",
            'cpx': "p$_x$",
            'cpy': "p$_y$",
            'pz': "p$_z$",
            'px': "p$_x$",
            'py': "p$_y$",
        }
        return formats

    def addBeamToPlotter(self, name, beam):
        """
        :param name: Name of beam that is being added to plotter
        :type name: str
        :param beam: Beam object that is being added to plotter
        :type beam: GeneralBeam
        :return: None, appends to a dict of named beams
        """
        self.beams[name] = beam

    def clearBeamFromPlotter(self, name):
        del self.beams[name]

    def clearAllBeams(self):
        self.beams = {}

    def get_unit_scale_str(self, unit):
        assert type(unit) == str
        assert hasattr(constants, unit)

        val = getattr(constants, unit)

        unit_sub = unit[0]

        if val > 1:
            unit_sub = unit_sub.upper()

        # edge cases
        # it is stupid the system works this way
        if val == 1000:
            unit_sub = 'k'

        if val == 1e-6:
            unit_sub = '\u03BC'  # mu

        return unit_sub

    def variable_formatter(self, variable):
        if variable in self.variable_formats:
            return self.variable_formats[variable]
        else:
            return variable

    def plotLPS(self, z_scale='micro', p_scale='mega', **kwargs):
        ax = self.general_scatter_plot(x_val='zn', x_scale=z_scale, y_val='cpz', y_scale=p_scale, **kwargs)
        return ax

    def plotLPS_t(self, t_scale='pico', p_scale='mega', **kwargs):
        ax = self.general_scatter_plot(x_val='t', x_scale=t_scale, y_val='cpz', y_scale=p_scale, **kwargs)

        ax.set_xlabel(r't ' + '(' + self.get_unit_scale_str(t_scale) + 's)')
        ax.set_ylabel(r'p$_z$ ' + '(' + self.get_unit_scale_str(p_scale) + 'eV/c)')

        return ax

    def plotHPS(self, x_scale='milli', xp_scale='milli', **kwargs):
        ax = self.general_scatter_plot(x_val='x', x_scale=x_scale, y_val='xp', y_scale=xp_scale, **kwargs)

        ax.set_xlabel(r'x ' + '(' + self.get_unit_scale_str(x_scale) + 'm)')
        ax.set_ylabel('x\' ' + '(' + self.get_unit_scale_str(xp_scale) + 'rad)')

        return ax

    def plotVPS(self, y_scale='milli', yp_scale='milli', **kwargs):
        ax = self.general_scatter_plot(x_val='y', x_scale=y_scale, y_val='yp', y_scale=yp_scale, **kwargs)

        ax.set_xlabel(r'y ' + '(' + self.get_unit_scale_str(y_scale) + 'm)')
        ax.set_ylabel('y\' ' + '(' + self.get_unit_scale_str(yp_scale) + 'rad)')
        return ax

    # TODO get the general plots/combined plots with the nice labels as above

    def plot_transverse_phase_space(self, x_scales=(constants.milli, constants.milli),
                                    y_scales=(constants.milli, constants.milli), **kwargs):
        fig, axes = self.general_multi_plot(x_vals=('x', 'y'), y_vals=('xp', 'yp'), x_scales=x_scales,
                                            y_scales=y_scales, **kwargs)
        return fig, axes

    def plot_transverse_properties(self, x_scales=('milli', 'milli', 'milli'), y_scales=('milli', 'milli', 'milli'),
                                   **kwargs):
        fig, axes = self.general_multi_plot(x_vals=('x', 'x', 'y'), y_vals=('y', 'xp', 'yp'), x_scales=x_scales,
                                            y_scales=y_scales, **kwargs)
        return fig, axes

    def plot_physical_projections(self, x_scales=('micro', 'micro', 'micro'), y_scales=('micro', 'micro', 'micro'),
                                  **kwargs):
        fig, axes = self.general_multi_plot(x_vals=('x', 'zn', 'zn'), y_vals=('y', 'x', 'y'), x_scales=x_scales,
                                            y_scales=y_scales, **kwargs)
        return fig, axes

    def general_scatter_plot(self, x_val, x_scale, y_val, y_scale, alpha=0.3, beams_to_plot="all",plotcolors=plt.rcParams['axes.prop_cycle'].by_key()['color']):

        assert type(x_val == str)
        assert type(y_val == str)
        fig, ax = plt.subplots()
        legend = []

        if beams_to_plot == "all":
            names_to_plot = list(self.beams.keys())

        else:
            assert type(beams_to_plot) == list

            for name in beams_to_plot:
                assert name in self.beams

            names_to_plot = beams_to_plot

        set_x_scale = getattr(constants, x_scale)
        set_y_scale = getattr(constants, y_scale)
        j=0;
        for name, beam in self.beams.items():
            if name in names_to_plot:
                ax.plot(getattr(beam, x_val) / set_x_scale, getattr(beam, y_val) / set_y_scale, '.', alpha=alpha,
                        markersize=1, color=plotcolors[j])
                legend.append(name)
                j+=1

        lgnd = ax.legend(legend)

        self.format_legends_size_opacity(lgnd)

        xvar = self.variable_formatter(x_val)
        x_unit = self.beams[names_to_plot[0]].parameter_units[x_val]
        yvar = self.variable_formatter(y_val)
        y_unit = self.beams[names_to_plot[0]].parameter_units[y_val]

        ax.set_xlabel(xvar + ' (' + self.get_unit_scale_str(x_scale) + x_unit + ')')
        ax.set_ylabel(yvar + ' (' + self.get_unit_scale_str(y_scale) + y_unit + ')')
        return fig, ax

    def general_multi_plot(self, x_vals, x_scales, y_vals, y_scales, alpha=0.3, beams_to_plot="all", plotcolors=plt.rcParams['axes.prop_cycle'].by_key()['color']):

        assert type(x_vals == list)
        assert type(y_vals == list)

        assert len(x_vals) == len(y_vals)
        assert len(x_vals) == len(x_scales)
        assert len(x_vals) == len(y_scales)

        fig, ax_list = plt.subplots(1, len(x_vals))
        legend = []

        if beams_to_plot == "all":
            names_to_plot = list(self.beams.keys())
        else:
            assert type(beams_to_plot) == list

            for name in beams_to_plot:
                assert name in self.beams

            names_to_plot = beams_to_plot

        for i, x_val in enumerate(x_vals):

            y_val = y_vals[i]
            set_x_scale = getattr(constants, x_scales[i])
            set_y_scale = getattr(constants, y_scales[i])
            j=0;
            for name, beam in self.beams.items():
                if name in names_to_plot:
                    ax_list[i].plot(getattr(beam, x_val) / set_x_scale, getattr(beam, y_val) / set_y_scale, '.',
                                    alpha=alpha, markersize=1,label=name, color = plotcolors[j])
                    legend.append(name)
                    j+=1
            ax_list[i].legend(legend)

            xvar = self.variable_formatter(x_val)
            x_unit = self.beams[names_to_plot[0]].parameter_units[x_val]

            yvar = self.variable_formatter(y_val)
            y_unit = self.beams[names_to_plot[0]].parameter_units[y_val]

            ax_list[i].set_xlabel(xvar + ' (' + self.get_unit_scale_str(x_scales[i]) + x_unit + ')')
            ax_list[i].set_ylabel(yvar + ' (' + self.get_unit_scale_str(y_scales[i]) + y_unit + ')')
        #fig.tight_layout()

        for ax in ax_list:
            lgnd = ax.legend(legend)

            lgnd =ax.legend(legend)
            self.format_legends_size_opacity(lgnd)


        return fig, ax_list


    def general_density_plot(self, x_val, x_scale, y_val, y_scale, cmap='turbo', beams_to_plot="all",hspace=0, namelabel=True):
        assert type(x_val == str)
        assert type(y_val == str)


        if beams_to_plot == "all":
            names_to_plot = list(self.beams.keys())

        else:
            assert type(beams_to_plot) == list

            for name in beams_to_plot:
                assert name in self.beams

            names_to_plot = beams_to_plot

        fig, ax_list = plt.subplots(len(names_to_plot), sharex=True)
        set_x_scale = getattr(constants, x_scale)
        set_y_scale = getattr(constants, y_scale)

        i=0
        for name, beam in self.beams.items():
            if name in names_to_plot:
                myPDF, ax = fastKDE.pdf(getattr(beam, x_val) / set_x_scale,getattr(beam, y_val) / set_y_scale)
                v1, v2 = ax
                if len(names_to_plot) == 1:
                    ax_list.pcolormesh(v1,v2,myPDF,cmap=cmap)
                    if namelabel == True:
                        ax_list.text(0.05,.9,name,horizontalalignment='left',transform=ax_list.transAxes, color='white', fontweight='bold')
                elif len(names_to_plot) > 1:
                    ax_list[i].pcolormesh(v1,v2,myPDF,cmap=cmap)
                    if namelabel==True:
                        ax_list[i].text(0.01,.85,name,horizontalalignment='left',transform=ax_list[i].transAxes, color='white', fontweight='bold')
                    i+=1

        xvar = self.variable_formatter(x_val)
        x_unit = self.beams[names_to_plot[0]].parameter_units[x_val]
        yvar = self.variable_formatter(y_val)
        y_unit = self.beams[names_to_plot[0]].parameter_units[y_val]

        if len(names_to_plot) > 1:
            for ax in ax_list:
                ax.set_facecolor(plt.cm.jet(1))
                ax.label_outer()
                ax.set_ylabel(yvar + ' (' + self.get_unit_scale_str(y_scale) + y_unit + ')')
            ax_list[-1].set_xlabel(xvar + ' (' + self.get_unit_scale_str(x_scale) + x_unit + ')')
        elif len(names_to_plot) == 1:
            ax_list.set_facecolor(plt.cm.jet(1))
            ax_list.label_outer()
            ax_list.set_ylabel(yvar + ' (' + self.get_unit_scale_str(y_scale) + y_unit + ')')
            ax_list.set_xlabel(xvar + ' (' + self.get_unit_scale_str(x_scale) + x_unit + ')')

        fig.subplots_adjust(hspace=hspace)
        fig.tight_layout()
        return fig,ax_list

    def general_profile_plot(self, x_val, x_scale, beams_to_plot="all", namelabel=True, plotcolors=plt.rcParams['axes.prop_cycle'].by_key()['color']):
        assert type(x_val == str)
        fig, ax_list = plt.subplots()
        legend = []

        if beams_to_plot == "all":
            names_to_plot = list(self.beams.keys())

        else:
            assert type(beams_to_plot) == list

            for name in beams_to_plot:
                assert name in self.beams

            names_to_plot = beams_to_plot

        set_x_scale = getattr(constants, x_scale)

        j=0;
        for name, beam in self.beams.items():
            if name in names_to_plot:
                myPDF, ax = fastKDE.pdf(getattr(beam, x_val) / set_x_scale)
                NormalisationFactor = 1 / max(myPDF)
                ax_list.plot(ax, NormalisationFactor * myPDF, color=plotcolors[j])
                legend.append(name)
                j+=1

        lgnd = ax_list.legend(legend)

        self.format_legends_size_opacity(lgnd)


        xvar = self.variable_formatter(x_val)
        x_unit = self.beams[names_to_plot[0]].parameter_units[x_val]

        ax_list.set_xlabel(xvar + ' (' + self.get_unit_scale_str(x_scale) + x_unit + ')')
        ax_list.set_ylabel('Charge (arb.)')
        return fig,ax_list


    def format_legends_size_opacity(self,legend,marker_size=6,alpha=1):
        """
        For setting plot legends
        :param legend: matplotlib ax legend
        :return:
        """
        for legendmarkers in legend.legendHandles:

            try:
                legendmarkers._legmarker.set_markersize(6)
                legendmarkers._legmarker.set_alpha(1)

            except AttributeError:
                # Toms fix - runs on latest matplotlib and windows...
                legendmarkers.set_markersize(6)
                legendmarkers.set_alpha(1)

            # Toby's fix
            # legendmarkers._legmarker.set_markersize(6)
            # legendmarkers._legmarker.set_alpha(1)


# -----------------------------------------------------------------------------

class BeamSlicer(object):
    '''
    Class to slice up a beam object
    Should return beam objects, so they can be drifted etc and turned into an image

    '''

    def __init__(self, beam, slice_parameter: str):
        """

        :param beam: beam object that holds macros and other items
        :type beam: GeneralBeam
        """
        # TODO docstrings
        self.beam = beam

        self.parameter = slice_parameter
        assert self.parameter in self.beam._beam

        # store the slices in an array...
        self.sliced_beams = []

        self.min_p = None
        self.max_p = None
        self.N_slices = None
        self.slice_width = None
        self.slice_centres = None

    def set_slice_parameters(self, min_p, max_p, N_slices, slice_width):
        # TODO add min and max base don min and max of array and N slice specified,
        # Todo slice width based on non-overlapping if not specified
        assert N_slices > 1, "Number of slices must be > 1, undefined behaviour otherwise."
        self.min_p = min_p
        self.max_p = max_p
        self.N_slices = N_slices
        self.slice_width = slice_width

        self.slice_centres = np.linspace(self.min_p, self.max_p, self.N_slices, endpoint=True)

    def perform_slicing(self):
        half_width = self.slice_width / 2

        for slice_cen in self.slice_centres:
            self.add_beam_slice(p_min=slice_cen - half_width, p_max=slice_cen + half_width)

    def add_beam_slice(self, p_min, p_max):
        # check the parameter to slice is in the beam dict!
        assert p_min < p_max

        selected_macros_index = np.where(np.logical_and(self.beam._beam[self.parameter] > p_min,
                                                        self.beam._beam[self.parameter] < p_max))

        print("selected this many macros", np.shape(selected_macros_index))

        # TODO Find efficient way to do this - beam from array
        sliced_beam = self.beam.returnBeam()

        sliced_beam._beam['x'] = self.beam._beam['x'][selected_macros_index]
        sliced_beam._beam['px'] = self.beam._beam['px'][selected_macros_index]

        sliced_beam._beam['y'] = self.beam._beam['y'][selected_macros_index]
        sliced_beam._beam['py'] = self.beam._beam['py'][selected_macros_index]

        sliced_beam._beam['z'] = self.beam._beam['z'][selected_macros_index]
        sliced_beam._beam['pz'] = self.beam._beam['pz'][selected_macros_index]

        self.sliced_beams.append(sliced_beam)
