import os, time, csv, sys, subprocess
import copy
import h5py
import numpy as np
import scipy.constants as constants
from scipy.spatial.distance import cdist
from scipy.spatial import ConvexHull
from scipy.stats import gaussian_kde
from itertools import compress
try:
    import sdds
except:
    print('sdds failed to load')
    pass
#sys.path.append(os.path.dirname(__file__))
#import minimumVolumeEllipse as mve
##import read_gdf_file as rgf
#MVE = mve.EllipsoidTool()

class beam(object):

    particle_mass = constants.m_e
    E0 = particle_mass * constants.speed_of_light**2
    E0_eV = E0 / constants.elementary_charge
    q_over_c = (constants.elementary_charge / constants.speed_of_light)
    speed_of_light = constants.speed_of_light

    def __init__(self, sddsindex=0):
        self.beam = {}
        self.sddsindex = sddsindex

    def set_particle_mass(self, mass=constants.m_e):
        self.particle_mass = mass

    def normalise_to_ref_particle(self, array, index=0,subtractmean=False):
        array = copy.copy(array)
        array[1:] = array[0] + array[1:]
        if subtractmean:
            array = array - array[0]#np.mean(array)
        return array

    def reset_dicts(self):
        self.beam = {}
        self.twiss = {}
        self.slice = {}
        self._tbins = []
        self._pbins = []

    def read_SDDS_beam_file(self, fileName, charge=None, ascii=False):
        self.reset_dicts()
        self.sdds = sdds.SDDS(self.sddsindex)
        self.sdds.load(fileName)
        for col in range(len(self.sdds.columnName)):
            if len(self.sdds.columnData[col]) == 1:
                self.beam[self.sdds.columnName[col]] = np.array(self.sdds.columnData[col][0])
            else:
                self.beam[self.sdds.columnName[col]] = np.array(self.sdds.columnData[col])
        self.SDDSparameters = dict()
        for param in range(len(self.sdds.parameterName)):
            self.SDDSparameters[self.sdds.parameterName[param]] = self.sdds.parameterData[param]
        # print 'self.SDDSparameterNames = ', self.SDDSparameterNames
        self.beam['code'] = "SDDS"
        cp = self.beam['p'] * self.E0_eV
        cpz = cp / np.sqrt(self.beam['xp']**2 + self.beam['yp']**2 + 1)
        cpx = self.beam['xp'] * cpz
        cpy = self.beam['yp'] * cpz
        self.beam['px'] = cpx * self.q_over_c
        self.beam['py'] = cpy * self.q_over_c
        self.beam['pz'] = cpz * self.q_over_c
        self.beam['t'] = self.beam['t']
        self.beam['z'] = (-1*self.Bz * constants.speed_of_light) * (self.t-np.mean(self.t)) #np.full(len(self.t), 0)
        if 'Charge' in self.SDDSparameters and len(self.SDDSparameters['Charge']) > 0:
            self.beam['total_charge'] = self.SDDSparameters['Charge'][0]
        elif charge is None:
            self.beam['total_charge'] = 0
        else:
            self.beam['total_charge'] = charge
        self.beam['charge'] = []

    def write_SDDS_file(self, filename, ascii=False, xyzoffset=[0,0,0]):
        """Save an SDDS file using the SDDS class."""
        xoffset = xyzoffset[0]
        yoffset = xyzoffset[1]
        zoffset = xyzoffset[2] # Don't think I need this because we are using t anyway...
        x = sdds.SDDS(self.sddsindex)
        if ascii:
            x.mode = x.SDDS_ASCII
        else:
            x.mode = x.SDDS_BINARY
        # {x, xp, y, yp, t, p, particleID}
        Cnames = ["x", "xp", "y", "yp", "t","p"]
        Ccolumns = ['x', 'xp', 'y', 'yp', 't', 'BetaGamma']
        Ctypes = [x.SDDS_DOUBLE, x.SDDS_DOUBLE, x.SDDS_DOUBLE, x.SDDS_DOUBLE, x.SDDS_DOUBLE, x.SDDS_DOUBLE]
        Csymbols = ["", "x'","","y'","",""]
        Cunits = ["m","","m","","s","m$be$nc"]
        Ccolumns = [np.array(self.x) - float(xoffset), self.xp, np.array(self.y) - float(yoffset), self.yp, self.t , self.cp/self.E0_eV]
        # {Step, pCentral, Charge, Particles, IDSlotsPerBunch}
        Pnames = ["pCentral", "Charge", "Particles"]
        Ptypes = [x.SDDS_DOUBLE, x.SDDS_DOUBLE, x.SDDS_LONG]
        Psymbols = ["p$bcen$n", "", ""]
        Punits = ["m$be$nc", "C", ""]
        parameterData = [[np.mean(self.BetaGamma)], [abs(self.beam['total_charge'])], [len(self.x)]]
        for i in range(len(Ptypes)):
            x.defineParameter(Pnames[i], Psymbols[i], Punits[i],"","", Ptypes[i], "")
            x.setParameterValueList(Pnames[i], parameterData[i])
        for i in range(len(Ctypes)):
            # name, symbol, units, description, formatString, type, fieldLength
            x.defineColumn(Cnames[i], Csymbols[i], Cunits[i],"","", Ctypes[i], 0)
            x.setColumnValueLists(Cnames[i], [list(Ccolumns[i])])
        x.save(filename)

    def set_beam_charge(self, charge):
        self.beam['total_charge'] = charge

    def read_csv_file(self, file, delimiter=' '):
        with open(file, 'r') as f:
            data = np.array([l for l in csv.reader(f, delimiter=delimiter,  quoting=csv.QUOTE_NONNUMERIC, skipinitialspace=True)])
        return data

    def write_csv_file(self, file, data):
        if sys.version_info[0] > 2:
            with open(file, 'w', newline='') as f:
                writer = csv.writer(f, delimiter=' ',  quoting=csv.QUOTE_NONNUMERIC, skipinitialspace=True)
                [writer.writerow(l) for l in data]
        else:
            with open(file, 'wb') as f:
                writer = csv.writer(f, delimiter=' ',  quoting=csv.QUOTE_NONNUMERIC, skipinitialspace=True)
                [writer.writerow(l) for l in data]

    def read_astra_beam_file(self, file, normaliseZ=False):
        starttime = time.time()
        self.reset_dicts()
        data = self.read_csv_file(file)
        # datanp = np.loadtxt(file)
        self.interpret_astra_data(data, normaliseZ=normaliseZ)

    # def read_hdf5_beam(self, data):
    #     self.reset_dicts()
    #     self.interpret_astra_data(data)

    def interpret_astra_data(self, data, normaliseZ=False):
        x, y, z, cpx, cpy, cpz, clock, charge, index, status = np.transpose(data)
        self.beam['code'] = "ASTRA"
        self.beam['reference_particle'] = data[0]
        # if normaliseZ:
        #     self.beam['reference_particle'][2] = 0
        self.beam['longitudinal_reference'] = 'z'
        znorm = self.normalise_to_ref_particle(z, subtractmean=True)
        z = self.normalise_to_ref_particle(z, subtractmean=False)
        cpz = self.normalise_to_ref_particle(cpz, subtractmean=False)
        clock = self.normalise_to_ref_particle(clock, subtractmean=True)
        cp = np.sqrt(cpx**2 + cpy**2 + cpz**2)
        self.beam['x'] = x
        self.beam['y'] = y
        self.beam['z'] = z
        self.beam['px'] = cpx * self.q_over_c
        self.beam['py'] = cpy * self.q_over_c
        self.beam['pz'] = cpz * self.q_over_c
        self.beam['clock'] = 1.0e-9*clock
        self.beam['charge'] = 1.0e-9*charge
        self.beam['index'] = index
        self.beam['status'] = status
        # print self.Bz
        self.beam['t'] = [clock if status == -1 else (z / (-1 * Bz * constants.speed_of_light)) for status, z, Bz, clock in zip(self.beam['status'], znorm, self.Bz, self.beam['clock'])]
        # self.beam['t'] = self.z / (1 * self.Bz * constants.speed_of_light)#[time if status is -1 else 0 for time, status in zip(clock, status)]#
        self.beam['total_charge'] = np.sum(self.beam['charge'])

    def read_csrtrack_beam_file(self, file):
        self.reset_dicts()
        data = self.read_csv_file(file)
        self.beam['code'] = "CSRTrack"
        self.beam['reference_particle'] = data[0]
        self.beam['longitudinal_reference'] = 'z'
        z, x, y, cpz, cpx, cpy, charge = np.transpose(data[1:])
        z = self.normalise_to_ref_particle(z, subtractmean=False)
        cpz = self.normalise_to_ref_particle(cpz, subtractmean=False)
        cp = np.sqrt(cpx**2 + cpy**2 + cpz**2)
        self.beam['x'] = x
        self.beam['y'] = y
        self.beam['z'] = z
        self.beam['px'] = cpx * self.q_over_c
        self.beam['py'] = cpy * self.q_over_c
        self.beam['pz'] = cpz * self.q_over_c
        self.beam['clock'] = np.full(len(self.x), 0)
        self.beam['clock'][0] = data[0, 0] * 1e-9
        self.beam['index'] = np.full(len(self.x), 5)
        self.beam['status'] = np.full(len(self.x), 1)
        self.beam['t'] = self.z / (-1 * self.Bz * constants.speed_of_light)# [time if status is -1 else 0 for time, status in zip(clock, self.beam['status'])]
        self.beam['charge'] = charge
        self.beam['total_charge'] = np.sum(self.beam['charge'])

    def read_vsim_h5_beam_file(self, filename, charge=70e-12, interval=1):
        self.reset_dicts()
        with h5py.File(filename, "r") as h5file:
            data = np.array(h5file.get('/BeamElectrons'))[1:-1:interval]
            z, y, x, cpz, cpy, cpx = data.transpose()
        self.beam['code'] = "VSIM"
        self.beam['longitudinal_reference'] = 'z'
        cp = np.sqrt(cpx**2 + cpy**2 + cpz**2)
        self.beam['x'] = x
        self.beam['y'] = y
        self.beam['z'] = z
        self.beam['px'] = cpx * self.particle_mass
        self.beam['py'] = cpy * self.particle_mass
        self.beam['pz'] = cpz * self.particle_mass
        self.beam['t'] = [(z / (1 * Bz * constants.speed_of_light)) for z, Bz in zip(self.z, self.Bz)]
        # self.beam['t'] = self.z / (1 * self.Bz * constants.speed_of_light)#[time if status is -1 else 0 for time, status in zip(clock, status)]#
        self.beam['total_charge'] = charge
        self.beam['charge'] = []

    def read_pacey_beam_file(self, file, charge=250e-12):
        self.reset_dicts()
        data = self.read_csv_file(file, delimiter='\t')
        self.beam['code'] = "TPaceyASTRA"
        self.beam['longitudinal_reference'] = 'z'
        x, y, z, cpx, cpy, cpz = np.transpose(data)
        cp = np.sqrt(cpx**2 + cpy**2 + cpz**2)
        self.beam['x'] = x
        self.beam['y'] = y
        self.beam['z'] = z
        self.beam['px'] = cpx * self.q_over_c
        self.beam['py'] = cpy * self.q_over_c
        self.beam['pz'] = cpz * self.q_over_c
        self.beam['t'] = [(z / (1 * Bz * constants.speed_of_light)) for z, Bz in zip(self.z, self.Bz)]
        # self.beam['t'] = self.z / (1 * self.Bz * constants.speed_of_light)#[time if status is -1 else 0 for time, status in zip(clock, status)]#
        self.beam['total_charge'] = charge
        self.beam['charge'] = []

    def convert_csrtrackfile_to_astrafile(self, csrtrackfile, astrafile):
        data = self.read_csv_file(csrtrackfile)
        z, x, y, cpz, cpx, cpy, charge = np.transpose(data[1:])
        charge = -charge*1e9
        clock0 = (data[0, 0] / constants.speed_of_light) * 1e9
        clock = np.full(len(x), 0)
        clock[0] = clock0
        index = np.full(len(x), 1)
        status = np.full(len(x), 5)
        array = np.array([x, y, z, cpx, cpy, cpz, clock, charge, index, status]).transpose()
        self.write_csv_file(astrafile, array)

    def find_nearest_vector(self, nodes, node):
        return cdist([node], nodes).argmin()

    def rms(self, x, axis=None):
        return np.sqrt(np.mean(x**2, axis=axis))

    def create_ref_particle(self, array, index=0, subtractmean=False):
        array[1:] = array[0] + array[1:]
        if subtractmean:
            array = array - np.mean(array)
        return array

    def write_astra_beam_file(self, file, index=1, status=5, charge=None, normaliseZ=False):
        if not isinstance(index,(list, tuple, np.ndarray)):
            if len(self.beam['charge']) == len(self.x):
                chargevector = 1e9*self.beam['charge']
            else:
                chargevector = np.full(len(self.x), 1e9*self.charge/len(self.x))
        if not isinstance(index,(list, tuple, np.ndarray)):
            indexvector = np.full(len(self.x), index)
        statusvector = self.beam['status'] if 'status' in self.beam else status if isinstance(status,(list, tuple, np.ndarray)) else np.full(len(self.x), status)
        ''' if a particle is emitting from the cathode it's z value is 0 and it's clock value is finite, otherwise z is finite and clock is irrelevant (thus zero) '''
        if self.beam['longitudinal_reference'] == 't':
            zvector = [0 if status == -1 and t == 0 else z for status, z, t in zip(statusvector, self.z, self.t)]
        else:
            zvector = self.z
        ''' if the clock value is finite, we calculate it from the z value, using Betaz '''
        # clockvector = [1e9*z / (1 * Bz * constants.speed_of_light) if status == -1 and t == 0 else 1.0e9*t for status, z, t, Bz in zip(statusvector, self.z, self.t, self.Bz)]
        clockvector = [1.0e9*t for status, z, t, Bz in zip(statusvector, self.z, self.t, self.Bz)]
        ''' this is the ASTRA array in all it's glory '''
        array = np.array([self.x, self.y, zvector, self.cpx, self.cpy, self.cpz, clockvector, chargevector, indexvector, statusvector]).transpose()
        if 'reference_particle' in self.beam:
            ref_particle = self.beam['reference_particle']
            # print 'we have a reference particle! ', ref_particle
            # np.insert(array, 0, ref_particle, axis=0)
        else:
            ''' take the rms - if the rms is 0 set it to 1, so we don't get a divide by error '''
            rms_vector = [a if abs(a) > 0 else 1 for a in self.rms(array, axis=0)]
            ''' normalise the array '''
            norm_array = array / rms_vector
            ''' take the meen of the normalised array '''
            mean_vector = np.mean(norm_array, axis=0)
            ''' find the index of the vector that is closest to the mean - if you read in an ASTRA file, this should actually return the reference particle! '''
            nearest_idx = self.find_nearest_vector(norm_array, mean_vector);
            ref_particle = array[nearest_idx]
            ''' set the closest mean vector to be in position 0 in the array '''
            array = np.roll(array, -1*nearest_idx, axis=0)

        ''' normalise Z to the reference particle '''
        array[1:,2] = array[1:,2] - ref_particle[2]
        ''' should we leave Z as the reference value, set it to 0, or set it to be some offset? '''
        if not normaliseZ is False:
            array[0,2] = 0
        if not isinstance(normaliseZ,(bool)):
            array[0,2] += normaliseZ
        ''' normalise pz and the clock '''
        # print('Mean pz = ', np.mean(array[:,5]))
        array[1:,5] = array[1:,5] - ref_particle[5]
        array[0,6] = array[0,6] + ref_particle[6]
        np.savetxt(file, array, fmt=('%.12e','%.12e','%.12e','%.12e','%.12e','%.12e','%.12e','%.12e','%d','%d'))

    def write_vsim_beam_file(self, file, normaliseT=False):
        if len(self.beam['charge']) == len(self.x):
            chargevector = self.beam['charge']
        else:
            chargevector = np.full(len(self.x), self.beam['total_charge']/len(self.x))
        if normaliseT:
            tvector = self.t - np.mean(self.t)
            zvector = self.z - np.mean(self.z)
        else:
            tvector = self.t
            zvector = self.z
        zvector = [t * (1 * Bz * constants.speed_of_light) if z == 0 else z for z, t, Bz in zip(zvector, tvector, self.Bz)]
        ''' this is the VSIM array in all it's glory '''
        array = np.array([zvector, self.y, self.x, self.Bz*self.gamma*constants.speed_of_light, self.By*self.gamma*constants.speed_of_light, self.Bx*self.gamma*constants.speed_of_light]).transpose()
        ''' take the rms - if the rms is 0 set it to 1, so we don't get a divide by error '''
        np.savetxt(file, array, fmt=('%.12e','%.12e','%.12e','%.12e','%.12e','%.12e'))

    def write_gdf_beam_file(self, filename, normaliseZ=False):
        q = np.full(len(self.x), -1 * constants.elementary_charge)
        m = np.full(len(self.x), constants.electron_mass)
        nmacro = np.full(len(self.x), abs(self.beam['total_charge'] / constants.elementary_charge / len(self.x)))
        toffset = np.mean(self.z / (self.Bz * constants.speed_of_light))
        z = self.z if not normaliseZ else (self.z - np.mean(self.z))
        dataarray = np.array([self.x, self.y, z, q, m, nmacro, self.gamma*self.Bx, self.gamma*self.By, self.gamma*self.Bz]).transpose()
        namearray = 'x y z q m nmacro GBx GBy GBz'
        np.savetxt(filename, dataarray, fmt=('%.12e','%.12e','%.12e','%.12e','%.12e','%.12e','%.12e','%.12e','%.12e'), header=namearray, comments='')

    #def read_gdf_beam_file_object(self, file):
    #    if isinstance(file, (str)):
    #        gdfbeam = rgf.read_gdf_file(file)
    #    elif isinstance(file, (rgf.read_gdf_file)):
    #        gdfbeam = file
    #    else:
    #        raise Exception('file is not str or gdf object!')
    #    return gdfbeam

    def calculate_gdf_s(self, file):
        gdfbeam = self.read_gdf_beam_file_object(file)
        datagrab = gdfbeam.get_grab(0)
        avgt = [datagrab.avgt]
        position = [datagrab.position]
        sposition = list(zip(*list(sorted(zip(avgt[0], position[0])))))[1]
        ssposition = list(zip(sposition, list(sposition[1:])+[0]))
        offset = 0
        spos = []
        for p1,p2 in ssposition:
            spos += [p1 + offset]
            if p2 < p1:
                offset += p1
        return spos

    def calculate_gdf_eta(self, file):
        gdfbeam = self.read_gdf_beam_file_object(file)
        etax = []
        etaxp = []
        tp = []
        for p in gdfbeam.positions:
            self.read_gdf_beam_file(gdfbeam=gdfbeam, position=p)
            if len(self.x) > 0:
                e, ep, t = self.calculate_etax()
                etax += [e]
                etaxp += [ep]
                tp += [t]
        etax, etaxp = list(zip(*list(sorted(zip(tp, etax, etaxp)))))[1:]
        return etax, etaxp

    def read_gdf_beam_file_info(self, file):
        self.reset_dicts()
        gdfbeamdata = None
        gdfbeam = self.read_gdf_beam_file_object(file)
        print('grab_groups = ',  gdfbeam.grab_groups)
        print(( 'Positions = ', gdfbeam.positions))
        print(( 'Times = ', gdfbeam.times))

    def read_gdf_beam_file(self, file=None, position=None, time=None, block=None, charge=None, longitudinal_reference='t', gdfbeam=None):
        self.reset_dicts()
        if gdfbeam is None and not file is None:
            gdfbeam = self.read_gdf_beam_file_object(file)
        elif gdfbeam is None and file is None:
            return None

        if position is not None:# and (time is not None or block is not None):
            # print 'Assuming position over time!'
            self.beam['longitudinal_reference'] = 't'
            gdfbeamdata = gdfbeam.get_position(position)
            if gdfbeamdata is not None:
                # print('GDF found position ', position)
                time = None
                block = None
            else:
                print('GDF DID NOT find position ', position)
                position = None
        elif position is None and time is not None and block is not None:
            # print 'Assuming time over block!'
            self.beam['longitudinal_reference'] = 'p'
            gdfbeamdata = gdfbeam.get_time(time)
            if gdfbeamdata is not None:
                block = None
            else:
                 time = None
        elif position is None and time is None and block is not None:
            gdfbeamdata = gdfbeam.get_grab(block)
            if gdfbeamdata is None:
                block = None
        elif position is None and time is None and block is None:
            gdfbeamdata = gdfbeam.get_grab(0)
        self.beam['code'] = "GPT"
        self.beam['x'] = gdfbeamdata.x
        self.beam['y'] = gdfbeamdata.y
        if hasattr(gdfbeamdata,'z') and longitudinal_reference == 'z':
            # print( 'z!')
            # print(( gdfbeamdata.z))
            self.beam['z'] = gdfbeamdata.z
            self.beam['t'] = np.full(len(self.z), 0)# self.z / (-1 * self.Bz * constants.speed_of_light)
        elif hasattr(gdfbeamdata,'t') and longitudinal_reference == 't':
            # print( 't!')
            self.beam['t'] = gdfbeamdata.t
            self.beam['z'] = (-1 * gdfbeamdata.Bz * constants.speed_of_light) * (gdfbeamdata.t-np.mean(gdfbeamdata.t)) + gdfbeamdata.z
        self.beam['gamma'] = gdfbeamdata.G
        if hasattr(gdfbeamdata,'q') and  hasattr(gdfbeamdata,'nmacro'):
            self.beam['charge'] = gdfbeamdata.q * gdfbeamdata.nmacro
            self.beam['total_charge'] = np.sum(self.beam['charge'])
        else:
            if charge is None:
                self.beam['total_charge'] = 0
            else:
                self.beam['total_charge'] = charge
        # print(( self.beam['charge']))
        vx = gdfbeamdata.Bx * constants.speed_of_light
        vy = gdfbeamdata.By * constants.speed_of_light
        vz = gdfbeamdata.Bz * constants.speed_of_light
        velocity_conversion = 1 / (constants.m_e * self.beam['gamma'])
        self.beam['px'] = vx / velocity_conversion
        self.beam['py'] = vy / velocity_conversion
        self.beam['pz'] = vz / velocity_conversion
        return gdfbeam

    def rotate_beamXZ(self, theta, preOffset=[0,0,0], postOffset=[0,0,0]):
        preOffset=np.array(preOffset)
        postOffset=np.array(postOffset)

        rotation_matrix = np.array([[np.cos(theta), 0, np.sin(theta)], [0, 1, 0], [-1*np.sin(theta), 0, np.cos(theta)]])
        beam = np.array([self.x,self.y,self.z]).transpose()
        self.beam['x'],self.beam['y'],self.beam['z'] = (np.dot(beam-preOffset, rotation_matrix)-postOffset).transpose()

        beam = np.array([self.px, self.py, self.pz]).transpose()
        self.beam['px'], self.beam['py'], self.beam['pz'] = np.dot(beam, rotation_matrix).transpose()

        if 'reference_particle' in self.beam:
            beam = np.array([self.beam['reference_particle'][0], self.beam['reference_particle'][1], self.beam['reference_particle'][2]])
            self.beam['reference_particle'][0], self.beam['reference_particle'][1], self.beam['reference_particle'][2] = (np.dot([beam-preOffset], rotation_matrix)[0]-postOffset)
            # print 'rotated ref part = ', np.dot([beam-preOffset], rotation_matrix)[0]
            beam = np.array([self.beam['reference_particle'][3], self.beam['reference_particle'][4], self.beam['reference_particle'][5]])
            self.beam['reference_particle'][3], self.beam['reference_particle'][4], self.beam['reference_particle'][5] = np.dot([beam], rotation_matrix)[0]

        self.beam['rotation'] = theta
        self.beam['offset'] = preOffset

    def unrotate_beamXZ(self):
        offset = self.beam['offset'] if 'offset' in self.beam else np.array([0,0,0])
        if 'rotation' in self.beam or abs(self.beam['rotation']) > 0:
            self.rotate_beamXZ(-1*self.beam['rotation'], -1*offset)

    def write_HDF5_beam_file(self, filename, centered=False, mass=constants.m_e, sourcefilename=None, pos=None, rotation=None, longitudinal_reference='t', zoffset=0):
        # print('zoffset = ', zoffset, type(zoffset))
        if isinstance(zoffset,(list, np.ndarray)) and len(zoffset) == 3:
            xoffset = zoffset[0]
            yoffset = zoffset[1]
            zoffset = zoffset[2]
        else:
            xoffset = 0
            yoffset = 0
        # print('xoffset = ', xoffset)
        # print('yoffset = ', yoffset)
        # print('zoffset = ', zoffset)
        with h5py.File(filename, "w") as f:
            inputgrp = f.create_group("Parameters")
            if not 'total_charge' in self.beam or self.beam['total_charge'] == 0:
                self.beam['total_charge'] = np.sum(self.beam['charge'])
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
            inputgrp['total_charge'] = self.beam['total_charge']
            inputgrp['npart'] = len(self.x)
            inputgrp['centered'] = centered
            inputgrp['code'] = self.beam['code'] # This could be changed here to specify teh type maybe?
            inputgrp['particle_mass'] = mass
            beamgrp = f.create_group("beam")
            if 'reference_particle' in self.beam:
                beamgrp['reference_particle'] = self.beam['reference_particle']
            if 'status' in self.beam:
                beamgrp['status'] = self.beam['status']
            beamgrp['longitudinal_reference'] = longitudinal_reference
            if len(self.beam['charge']) == len(self.x):
                chargevector = self.beam['charge']
            else:
                chargevector = np.full(len(self.x), self.charge/len(self.x))
            array = np.array([self.x + xoffset, self.y + yoffset, self.z + zoffset, self.cpx, self.cpy, self.cpz, self.t, chargevector]).transpose()
            beamgrp['columns'] = np.array(['x','y','z', 'cpx', 'cpy', 'cpz', 't', 'q'], dtype='S')
            beamgrp['units'] = np.array(['m','m','m','eV','eV','eV','s','e'], dtype='S')
            beamgrp.create_dataset("beam", data=array)

    def read_HDF5_beam_file(self, filename, local=False):
        self.reset_dicts()
        with h5py.File(filename, "r") as h5file:
            if h5file.get('beam/reference_particle') is not None:
                self.beam['reference_particle'] = np.array(h5file.get('beam/reference_particle'))
            if h5file.get('beam/longitudinal_reference') is not None:
                self.beam['longitudinal_reference'] = np.array(h5file.get('beam/longitudinal_reference'))
            else:
                self.beam['longitudinal_reference'] = 't'
            if h5file.get('beam/status') is not None:
                self.beam['status'] = np.array(h5file.get('beam/status'))
            x, y, z, cpx, cpy, cpz, t, charge = np.array(h5file.get('beam/beam')).transpose()
            cp = np.sqrt(cpx**2 + cpy**2 + cpz**2)
            self.beam['x'] = x
            self.beam['y'] = y
            self.beam['z'] = z
            # self.beam['cpx'] = cpx
            # self.beam['cpy'] = cpy
            # self.beam['cpz'] = cpz
            self.beam['px'] = cpx * self.q_over_c
            self.beam['py'] = cpy * self.q_over_c
            self.beam['pz'] = cpz * self.q_over_c
            # self.beam['cp'] = cp
            # self.beam['p'] = cp * self.q_over_c
            # self.beam['xp'] = np.arctan(self.px/self.pz)
            # self.beam['yp'] = np.arctan(self.py/self.pz)
            self.beam['clock'] = np.full(len(self.x), 0)
            # self.beam['gamma'] = np.sqrt(1+(self.cp/self.E0_eV)**2)
            # velocity_conversion = 1 / (constants.m_e * self.gamma)
            # self.beam['vx'] = velocity_conversion * self.px
            # self.beam['vy'] = velocity_conversion * self.py
            # self.beam['vz'] = velocity_conversion * self.pz
            # self.beam['Bx'] = self.vx / constants.speed_of_light
            # self.beam['By'] = self.vy / constants.speed_of_light
            # self.beam['Bz'] = self.vz / constants.speed_of_light
            self.beam['t'] = t
            self.beam['charge'] = charge
            self.beam['total_charge'] = np.sum(self.beam['charge'])
            #startposition = np.array(h5file.get('/Parameters/Starting_Position'))
            #startposition = startposition if startposition is not None else [0,0,0]
            #self.beam['starting_position'] = startposition
            #theta = np.array(h5file.get('/Parameters/Rotation'))
            #theta = theta if theta is not None else 0
            #self.beam['rotation'] = theta
            #code = np.array(h5file.get('/Parameters/code')) # Added for completeness as its in the write method
            #self.beam['code'] = code.item().decode("utf-8") # wee bit hacky to get it to behave nicely
            if local == True:
                self.rotate_beamXZ(self.beam['rotation'], preOffset=self.beam['starting_position'])

    ''' ********************  Statistical Parameters  ************************* '''

    def kde_function(self, x, bandwidth=0.2, **kwargs):
        """Kernel Density Estimation with Scipy"""
        # Note that scipy weights its bandwidth by the covariance of the
        # input data.  To make the results comparable to the other methods,
        # we divide the bandwidth by the sample standard deviation here.
        # Taken from https://jakevdp.github.io/blog/2013/12/01/kernel-density-estimation/
        if not hasattr(self, '_kde_x') or not len(x) == len(self._kde_x) or not np.allclose(x, self._kde_x) or not bandwidth == self._kde_bandwidth:
            self._kde_x = x
            self._kde_bandwidth = bandwidth
            self._kde_function = gaussian_kde(x, bw_method=bandwidth / x.std(ddof=1), **kwargs)
        return self._kde_function

    def PDF(self, x, x_grid, bandwidth=0.2, **kwargs):
        kde = self.kde_function(x, bandwidth, **kwargs)
        return kde.evaluate(x_grid)

    def PDFI(self, x, x_grid, bandwidth=0.2, **kwargs):
        kde = self.kde_function(x, bandwidth, **kwargs)
        vals = kde.evaluate(x_grid)
        f = lambda bin, val: self.charge / len(self.t) * (val / bin)
        return vals#self.charge * vals / (2*abs(x_grid[1] - x_grid[0])) / len(self.t) #[f(x_grid[1] - x_grid[0], v) for v in vals]

    def CDF(self, x, x_grid, bandwidth=0.2, **kwargs):
        kde = self.kde_function(x, bandwidth, **kwargs)
        cdf = np.vectorize(lambda e: kde.integrate_box_1d(x_grid[0], e))
        return cdf(x_grid)

    def FWHM(self, X, Y, frac=0.5):
        frac = 1.0/frac if frac > 1 else frac
        d = Y - (max(Y) * frac)
        indexes = np.where(d > 0)[0]
        return abs(X[indexes][-1] - X[indexes][0]), indexes

    def covariance(self, u, up):
        u2 = u - np.mean(u)
        up2 = up - np.mean(up)
        return np.mean(u2*up2) - np.mean(u2)*np.mean(up2)

    def emittance(self, x, xp, p=None):
        emittance = np.sqrt(self.covariance(x, x)*self.covariance(xp, xp) - self.covariance(x, xp)**2)
        if p is None:
            return emittance
        else:
            gamma = np.mean(p)/self.E0_eV
            return gamma*emittance

    @property
    def volume(self):
        return self.volume6D(self.x, self.y, self.z-np.mean(self.z), self.cpx/self.cpz, self.cpy/self.cpz, ((self.cpz/np.mean(self.cp)) - 1))

    @property
    def density(self):
        return len(self.x) / self.volume

    def volume6D(self, x, y, t, xp, yp, cp):
        if len(x) < 10:
            return 1e6
        else:
            beam = list(zip(x, y, t, xp, yp, cp))
            return ConvexHull(beam, qhull_options='QJ').volume

    #def mve_emittance(self, x, xp, p=None):
    #    (center, radii, rotation, hullP) = MVE.getMinVolEllipse(list(zip(x,xp)), .01)
    #    emittance = radii[0] * radii[1]
    #    if p is None:
    #        return emittance
    #    else:
    #        gamma = np.mean(p)/self.E0_eV
    #        return gamma*emittance

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
        emiti = gamma * self.x**2 + 2 * alpha * self.x * self.xp + beta * self.xp * self.xp
        return sorted(emiti)[int(0.9*len(emiti)-0.5)]
    @property
    def normalized_horizontal_emittance_90(self):
        emit = self.horizontal_emittance_90
        return np.mean(self.cp)/self.E0_eV * emit
    @property
    def vertical_emittance_90(self):
        emit = self.vertical_emittance
        alpha = self.alpha_y
        beta = self.beta_y
        gamma = self.gamma_y
        emiti = gamma * self.y**2 + 2 * alpha * self.y * self.yp + beta * self.yp * self.yp
        return sorted(emiti)[int(0.9*len(emiti)-0.5)]
    @property
    def normalized_vertical_emittance_90(self):
        emit = self.vertical_emittance_90
        return np.mean(self.cp)/self.E0_eV * emit

    @property
    def beta_x(self):
        self.twiss['beta_x'] = self.covariance(self.x,self.x) / self.horizontal_emittance
        return self.twiss['beta_x']
    @property
    def alpha_x(self):
        self.twiss['alpha_x'] = -1*self.covariance(self.x,self.xp) / self.horizontal_emittance
        return self.twiss['alpha_x']
    @property
    def gamma_x(self):
        self.twiss['gamma_x'] = self.covariance(self.xp,self.xp) / self.horizontal_emittance
        return self.twiss['gamma_x']
    @property
    def beta_y(self):
        self.twiss['beta_y'] = self.covariance(self.y,self.y) / self.vertical_emittance
        return self.twiss['beta_y']
    @property
    def alpha_y(self):
        self.twiss['alpha_y'] = -1*self.covariance(self.y,self.yp) / self.vertical_emittance
        return self.twiss['alpha_y']
    @property
    def gamma_y(self):
        self.twiss['gamma_y'] = self.covariance(self.yp,self.yp) / self.vertical_emittance
        return self.twiss['gamma_y']
    @property
    def twiss_analysis(self):
        return self.horizontal_emittance, self.alpha_x, self.beta_x, self.gamma_x, self.vertical_emittance, self.alpha_y, self.beta_y, self.gamma_y

    def eta_correlation(self, u):
        return self.covariance(u,self.p) / self.covariance(self.p, self.p)
    def eta_corrected(self, u):
        return u - self.eta_correlation(u)*self.p
    @property
    def horizontal_emittance_corrected(self):
        xc = self.eta_corrected(self.x)
        xpc = self.eta_corrected(self.xp)
        return self.emittance(xc, xpc)
    @property
    def vertical_emittance_corrected(self):
        yc = self.eta_corrected(self.y)
        ypc = self.eta_corrected(self.yp)
        return self.emittance(yc, ypc)
    @property
    def beta_x_corrected(self):
        xc = self.eta_corrected(self.x)
        self.twiss['beta_x'] = self.covariance(xc, xc) / self.horizontal_emittance_corrected
        return self.twiss['beta_x']
    @property
    def alpha_x_corrected(self):
        xc = self.eta_corrected(self.x)
        xpc = self.eta_corrected(self.xp)
        self.twiss['alpha_x'] = -1*self.covariance(xc, xpc) / self.horizontal_emittance_corrected
        return self.twiss['alpha_x']
    @property
    def gamma_x_corrected(self):
        xpc = self.eta_corrected(self.xp)
        self.twiss['gamma_x'] = self.covariance(xpc, xpc) / self.horizontal_emittance_corrected
        return self.twiss['gamma_x']
    @property
    def beta_y_corrected(self):
        yc = self.eta_corrected(self.y)
        self.twiss['beta_y'] = self.covariance(yc,yc) / self.vertical_emittance_corrected
        return self.twiss['beta_y']
    @property
    def alpha_y_corrected(self):
        yc = self.eta_corrected(self.y)
        ypc = self.eta_corrected(self.yp)
        self.twiss['alpha_y'] = -1*self.covariance(yc, ypc) / self.vertical_emittance_corrected
        return self.twiss['alpha_y']
    @property
    def gamma_y_corrected(self):
        ypc = self.eta_corrected(self.yp)
        self.twiss['gamma_y'] = self.covariance(ypc,ypc) / self.vertical_emittance_corrected
        return self.twiss['gamma_y']
    @property
    def twiss_analysis_corrected(self):
        return  self.horizontal_emittance_corrected, self.alpha_x_corrected, self.beta_x_corrected, self.gamma_x_corrected, \
                self.vertical_emittance_corrected, self.alpha_y_corrected, self.beta_y_corrected, self.gamma_y_corrected

    @property
    def slice_length(self):
        return self._slicelength

    @slice_length.setter
    def slice_length(self, slicelength):
        self._slicelength = slicelength

    @property
    def slices(self):
        return self._slices

    @slices.setter
    def slices(self, slices):
        twidth = (max(self.t) - min(self.t))
        if twidth == 0:
            t = self.z / (-1 * self.Bz * constants.speed_of_light)
            twidth = (max(t) - min(t))
        if slices == 0:
            slices = int(twidth / 0.1e-12)
        self._slices = slices
        self._slicelength = twidth / self._slices

    def bin_time(self):
        if not hasattr(self,'slice'):
            self.slice = {}
        if not hasattr(self,'_slicelength'):
            self.slice_length = 0.1e-12
            # print("Assuming slice length is 100 fs")
        twidth = (max(self.t) - min(self.t))
        if twidth == 0:
            t = self.z / (-1 * self.Bz * constants.speed_of_light)
            twidth = (max(t) - min(t))
        else:
            t = self.t
        if not self.slice_length > 0.0:
            self.slice_length = twidth / 20.0
        nbins = max([1,int(np.ceil(twidth / self.slice_length))])+2
        self._hist, binst =  np.histogram(t, bins=nbins, range=(min(t)-self.slice_length, max(t)+self.slice_length))
        self.slice['t_Bins'] = binst
        self._t_binned = np.digitize(t, self.slice['t_Bins'])
        self._tfbins = [[self._t_binned == i] for i in range(1, len(binst))]
        self._tbins = [np.array(self.t)[tuple(tbin)] for tbin in self._tfbins]
        self._cpbins = [np.array(self.cp)[tuple(tbin)] for tbin in self._tfbins]

    def bin_momentum(self, width=10**6):
        if not hasattr(self,'slice'):
            self.slice = {}
        pwidth = (max(self.cp) - min(self.cp))
        if width is None:
            self.slice_length_cp = pwidth / self.slices
        else:
            self.slice_length_cp = width
        nbins = max([1,int(np.ceil(pwidth / self.slice_length_cp))])+2
        self._hist, binst =  np.histogram(self.cp, bins=nbins, range=(min(self.cp)-self.slice_length_cp, max(self.cp)+self.slice_length_cp))
        self.slice['cp_Bins'] = binst
        self._cp_binned = np.digitize(self.cp, self.slice['cp_Bins'])
        self._tfbins = [np.array([self._cp_binned == i]) for i in range(1, len(binst))]
        self._cpbins = [self.cp[tuple(cpbin)] for cpbin in self._tfbins]
        self._tbins = [self.t[tuple(cpbin)] for cpbin in self._tfbins]

    @property
    def slice_bins(self):
        if not hasattr(self,'slice'):
            self.bin_time()
        bins = self.slice['t_Bins']
        return (bins[:-1] + bins[1:]) / 2
        # return [t.mean() for t in ]
    @property
    def slice_cpbins(self):
        if not hasattr(self,'slice'):
            self.bin_momentum()
        bins = self.slice['cp_Bins']
        return (bins[:-1] + bins[1:]) / 2
        # return [t.mean() for t in ]
    @property
    def slice_momentum(self):
        if not hasattr(self,'_tbins') or not hasattr(self,'_cpbins'):
            self.bin_time()
        self.slice['Momentum'] = np.array([cpbin.mean() for cpbin in self._cpbins])
        return self.slice['Momentum']
    @property
    def slice_momentum_spread(self):
        if not hasattr(self,'_tbins') or not hasattr(self,'_cpbins'):
            self.bin_time()
        self.slice['Momentum_Spread'] = np.array([cpbin.std() for cpbin in self._cpbins])
        return self.slice['Momentum_Spread']
    @property
    def slice_relative_momentum_spread(self):
        if not hasattr(self,'_tbins') or not hasattr(self,'_cpbins'):
            self.bin_time()
        self.slice['Relative_Momentum_Spread'] = np.array([100*cpbin.std()/cpbin.mean() for cpbin in self._cpbins])
        return self.slice['Relative_Momentum_Spread']

    def slice_data(self, data):
        return [data[tuple(tbin)] for tbin in self._tfbins]

    def emitbins(self, x, y):
        xbins = self.slice_data(x)
        ybins = self.slice_data(y)
        return list(zip(*[xbins, ybins, self._cpbins]))

    @property
    def slice_6D_Volume(self):
        if not hasattr(self,'_tbins') or not hasattr(self,'_cpbins'):
            self.bin_time()
        xbins = self.slice_data(self.x)
        ybins = self.slice_data(self.y)
        zbins = self.slice_data(self.z-np.mean(self.z))
        pxbins = self.slice_data(self.cpx/self.cpz)
        pybins = self.slice_data(self.cpy/self.cpz)
        pzbins = self.slice_data(((self.cpz/np.mean(self.cp)) - 1))
        emitbins = list(zip(xbins, ybins, zbins, pxbins, pybins, pzbins))
        self.slice['6D_Volume'] = np.array([self.volume6D(*a) for a in emitbins])
        return self.slice['6D_Volume']
    @property
    def slice_density(self):
        if not hasattr(self,'_tbins') or not hasattr(self,'_cpbins'):
            self.bin_time()
        xbins = self.slice_data(self.x)
        volume = self.slice_6D_Volume
        self.slice['Density'] = np.array([len(x)/v for x, v in zip(xbins, volume)])
        return self.slice['Density']
    @property
    def slice_horizontal_emittance(self):
        if not hasattr(self,'_tbins') or not hasattr(self,'_cpbins'):
            self.bin_time()
        emitbins = self.emitbins(self.x, self.xp)
        self.slice['Horizontal_Emittance'] = np.array([self.emittance(xbin, xpbin) for xbin, xpbin, cpbin in emitbins])
        return self.slice['Horizontal_Emittance']
    @property
    def slice_vertical_emittance(self):
        if not hasattr(self,'_tbins') or not hasattr(self,'_cpbins'):
            self.bin_time()
        emitbins = self.emitbins(self.y, self.yp)
        self.slice['Vertical_Emittance'] = np.array([self.emittance(ybin, ypbin) for ybin, ypbin, cpbin in emitbins])
        return self.slice['Vertical_Emittance']
    @property
    def slice_normalized_horizontal_emittance(self):
        if not hasattr(self,'_tbins') or not hasattr(self,'_cpbins'):
            self.bin_time()
        emitbins = self.emitbins(self.x, self.xp)
        self.slice['Normalized_Horizontal_Emittance'] = np.array([self.emittance(xbin, xpbin, cpbin) for xbin, xpbin, cpbin in emitbins])
        return self.slice['Normalized_Horizontal_Emittance']
    @property
    def slice_normalized_vertical_emittance(self):
        if not hasattr(self,'_tbins') or not hasattr(self,'_cpbins'):
            self.bin_time()
        emitbins = self.emitbins(self.y, self.yp)
        self.slice['Normalized_Vertical_Emittance'] = np.array([self.emittance(ybin, ypbin, cpbin) for ybin, ypbin, cpbin in emitbins])
        return self.slice['Normalized_Vertical_Emittance']
    @property
    def slice_normalized_mve_horizontal_emittance(self):
        if not hasattr(self,'_tbins') or not hasattr(self,'_cpbins'):
            self.bin_time()
        emitbins = self.emitbins(self.x, self.xp)
        self.slice['Normalized_mve_Horizontal_Emittance'] = np.array([self.mve_emittance(xbin, xpbin, cpbin) for xbin, xpbin, cpbin in emitbins])
        return self.slice['Normalized_mve_Horizontal_Emittance']
    @property
    def slice_normalized_mve_vertical_emittance(self):
        if not hasattr(self,'_tbins') or not hasattr(self,'_cpbins'):
            self.bin_time()
        emitbins = self.emitbins(self.y, self.yp)
        self.slice['Normalized_mve_Vertical_Emittance'] = np.array([self.mve_emittance(ybin, ypbin, cpbin) for ybin, ypbin, cpbin in emitbins])
        return self.slice['Normalized_mve_Vertical_Emittance']
    @property
    def slice_peak_current(self):
        if not hasattr(self,'_hist'):
            self.bin_time()
        f = lambda bin: self.charge / len(self.t) * (len(bin) / (max(bin) - min(bin))) if len(bin) > 1 else 0
        # f = lambda bin: len(bin) if len(bin) > 1 else 0
        self.slice['Peak_Current'] = np.array([f(bin) for bin in self._tbins])
        return abs(self.slice['Peak_Current'])
    @property
    def slice_max_peak_current_slice(self):
        peakI = self.slice_peak_current
        self.slice['Max_Peak_Current_Slice'] = list(abs(peakI)).index(max(abs(peakI)))
        return self.slice['Max_Peak_Current_Slice']

    @property
    def slice_beta_x(self):
        xbins = self.slice_data(self.beam['x'])
        exbins =  self.slice_horizontal_emittance
        emitbins = list(zip(xbins, exbins))
        self.slice['slice_beta_x'] = np.array([self.covariance(x, x)/ex for x, ex in emitbins])
        return self.slice['slice_beta_x']
    @property
    def slice_alpha_x(self):
        xbins = self.slice_data(self.x)
        xpbins = self.slice_data(self.xp)
        exbins =  self.slice_horizontal_emittance
        emitbins = list(zip(xbins, xpbins, exbins))
        self.slice['slice_alpha_x'] = np.array([-1*self.covariance(x, xp)/ex for x, xp, ex in emitbins])
        return self.slice['slice_alpha_x']
    @property
    def slice_gamma_x(self):
        self.twiss['gamma_x'] = self.covariance(self.xp,self.xp) / self.horizontal_emittance
        return self.twiss['gamma_x']
    @property
    def slice_beta_y(self):
        ybins = self.slice_data(self.beam['y'])
        eybins =  self.slice_vertical_emittance
        emitbins = list(zip(ybins, eybins))
        self.slice['slice_beta_y'] = np.array([self.covariance(y, y)/ey for y, ey in emitbins])
        return self.slice['slice_beta_y']
    @property
    def slice_alpha_y(self):
        ybins = self.slice_data(self.y)
        ypbins = self.slice_data(self.yp)
        eybins =  self.slice_vertical_emittance
        emitbins = list(zip(ybins, ypbins, eybins))
        self.slice['slice_alpha_y'] = np.array([-1*self.covariance(y,yp)/ey for y, yp, ey in emitbins])
        return self.twiss['slice_alpha_y']
    @property
    def slice_gamma_y(self):
        self.twiss['gamma_y'] = self.covariance(self.yp,self.yp) / self.vertical_emittance
        return self.twiss['gamma_y']

    def sliceAnalysis(self, density=False):
        self.slice = {}
        self.bin_time()
        peakIPosition = self.slice_max_peak_current_slice
        slice_density = self.slice_density[peakIPosition] if density else 0
        return self.slice_peak_current[peakIPosition], \
            np.std(self.slice_peak_current), \
            self.slice_relative_momentum_spread[peakIPosition], \
            self.slice_normalized_horizontal_emittance[peakIPosition], \
            self.slice_normalized_vertical_emittance[peakIPosition], \
            self.slice_momentum[peakIPosition], \
            self.slice_density[peakIPosition],

    def mvesliceAnalysis(self):
        self.slice = {}
        self.bin_time()
        peakIPosition = self.slice_max_peak_current_slice
        return self.slice_peak_current[peakIPosition], \
            np.std(self.slice_peak_current), \
            self.slice_relative_momentum_spread[peakIPosition], \
            self.slice_normalized_mve_horizontal_emittance[peakIPosition], \
            self.slice_normalized_mve_vertical_emittance[peakIPosition], \
            self.slice_momentum[peakIPosition], \
            self.slice_density[peakIPosition],

    @property
    def chirp(self):
        self.bin_time()
        slice_current_centroid_indices = []
        slice_momentum_centroid = []
        peakIPosition = self.slice_max_peak_current_slice
        peakI = self.slice_peak_current[peakIPosition]
        slicemomentum = self.slice_momentum
        for index, slice_current in enumerate(self.slice_peak_current):
            if abs(peakI - slice_current) < (peakI * 0.75):
                slice_current_centroid_indices.append(index)
        for index in slice_current_centroid_indices:
            slice_momentum_centroid.append(slicemomentum[index])
        chirp = (1e-18 * (slice_momentum_centroid[-1] - slice_momentum_centroid[0]) / (len(slice_momentum_centroid) * self.slice_length))
        return chirp

    @property
    def x(self):
        return self.beam['x']
    @property
    def y(self):
        return self.beam['y']
    @property
    def z(self):
        return self.beam['z']
    @property
    def zn(self):
        return self.beam['z']-np.mean(self.beam['z'])
    @property
    def px(self):
        return self.beam['px']
    @property
    def py(self):
        return self.beam['py']
    @property
    def pz(self):
        return self.beam['pz']
    @property
    def cpx(self):
        return self.beam['px'] / self.q_over_c
    @property
    def cpy(self):
        return self.beam['py'] / self.q_over_c
    @property
    def cpz(self):
        return self.beam['pz'] / self.q_over_c
    @property
    def xp(self):
        return np.arctan(self.px/self.pz)
    @property
    def yp(self):
        return np.arctan(self.py/self.pz)
    @property
    def t(self):
        return self.beam['t']
    @property
    def p(self):
        return self.cp * self.q_over_c
    @property
    def cp(self):
        return np.sqrt(self.cpx**2 + self.cpy**2 + self.cpz**2)
    @property
    def Brho(self):
        return np.mean(self.p) / constants.elementary_charge
    @property
    def gamma(self):
        return np.sqrt(1+(self.cp/self.E0_eV)**2)
    @property
    def BetaGamma(self):
        return self.cp/self.E0_eV
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
    def charge(self):
        return self.beam['total_charge']
    @property
    def sigma_z(self):
        return self.rms(self.Bz*constants.speed_of_light*(self.beam['t'] - np.mean(self.beam['t'])))
    @property
    def momentum_spread(self):
        return self.cp.std()/np.mean(self.cp)
    @property
    def linear_chirp_z(self):
        return -1*self.rms(self.Bz*constants.speed_of_light*self.t)/self.momentum_spread/100

    def computeCorrelations(self, x, y):
        xAve = np.mean(x)
        yAve = np.mean(y)
        C11 = 0
        C12 = 0
        C22 = 0
        for i, ii in enumerate(x):
            dx = x[i] - xAve
            dy = y[i] - yAve
            C11 += dx*dx
            C12 += dx*dy
            C22 += dy*dy
        C11 /= len(x)
        C12 /= len(x)
        C22 /= len(x)
        return C11, C12, C22

    def calculate_etax(self):
        p = self.cp
        pAve = np.mean(p)
        p = [a / pAve - 1 for a in p]
        S11, S16, S66 = self.computeCorrelations(self.x, self.cp)
        eta1 = -pAve * S16/S66 if S66 else 0
        S22, S26, S66 = self.computeCorrelations(self.xp, self.cp)
        etap1 = -pAve * S26/S66 if S66 else 0
        return eta1, etap1, np.mean(self.t)

    def performTransformation(self, x, xp, beta=False, alpha=False, nEmit=False):
        p = self.cp
        pAve = np.mean(p)
        p = [a / pAve - 1 for a in p]
        S11, S16, S66 = self.computeCorrelations(self.x, self.cp)
        eta1 = S16/S66 if S66 else 0
        S22, S26, S66 = self.computeCorrelations(self.xp, self.cp)
        etap1 = S26/S66 if S66 else 0
        for i, ii in enumerate(x):
            x[i] -= p[i] * eta1
            xp[i] -= p[i] * etap1

        S11, S12, S22 = self.computeCorrelations(x, xp)
        emit = np.sqrt(S11*S22 - S12**2)
        beta1 = S11/emit
        alpha1 = -S12/emit
        beta2 = beta if beta is not False else beta1
        alpha2 = alpha if alpha is not False else alpha1
        R11 = beta2/np.sqrt(beta1*beta2)
        R12 = 0
        R21 = (alpha1-alpha2)/np.sqrt(beta1*beta2)
        R22 = beta1/np.sqrt(beta1*beta2)
        if nEmit is not False:
            factor = np.sqrt(nEmit / (emit*pAve))
            R11 *= factor
            R12 *= factor
            R22 *= factor
            R21 *= factor
        for i, ii in enumerate(x):
            x0 = x[i]
            xp0 = xp[i]
            x[i] = R11 * x0 + R12 * xp0
            xp[i] = R21*x0 + R22*xp0
        return x, xp

    def rematchXPlane(self, beta=False, alpha=False, nEmit=False):
        x, xp = self.performTransformation(self.x, self.xp, beta, alpha, nEmit)
        self.beam['x'] = x
        self.beam['xp'] = xp

        cpz = self.cp / np.sqrt(self.beam['xp']**2 + self.yp**2 + 1)
        cpx = self.beam['xp'] * cpz
        cpy = self.yp * cpz
        self.beam['px'] = cpx * self.q_over_c
        self.beam['py'] = cpy * self.q_over_c
        self.beam['pz'] = cpz * self.q_over_c

    def rematchYPlane(self, beta=False, alpha=False, nEmit=False):
        y, yp = self.performTransformation(self.y, self.yp, beta, alpha, nEmit)
        self.beam['y'] = y
        self.beam['yp'] = yp

        cpz = self.cp / np.sqrt(self.xp**2 + self.beam['yp']**2 + 1)
        cpx = self.xp * cpz
        cpy = self.beam['yp'] * cpz
        self.beam['px'] = cpx * self.q_over_c
        self.beam['py'] = cpy * self.q_over_c
        self.beam['pz'] = cpz * self.q_over_c

    def performTransformationPeakISlice(self, xslice, xpslice, x, xp, beta=False, alpha=False, nEmit=False):
        p = self.cp
        pAve = np.mean(p)
        p = [a / pAve - 1 for a in p]
        S11, S16, S66 = self.computeCorrelations(self.x, self.cp)
        eta1 = S16/S66 if S66 else 0
        S22, S26, S66 = self.computeCorrelations(self.xp, self.cp)
        etap1 = S26/S66 if S66 else 0
        for i, ii in enumerate(x):
            x[i] -= p[i] * eta1
            xp[i] -= p[i] * etap1

        S11, S12, S22 = self.computeCorrelations(xslice, xpslice)
        emit = np.sqrt(S11*S22 - S12**2)
        beta1 = S11/emit
        alpha1 = -S12/emit
        beta2 = beta if beta is not False else beta1
        alpha2 = alpha if alpha is not False else alpha1
        R11 = beta2/np.sqrt(beta1*beta2)
        R12 = 0
        R21 = (alpha1-alpha2)/np.sqrt(beta1*beta2)
        R22 = beta1/np.sqrt(beta1*beta2)
        if nEmit is not False:
            factor = np.sqrt(nEmit / (emit*pAve))
            R11 *= factor
            R12 *= factor
            R22 *= factor
            R21 *= factor
        for i, ii in enumerate(x):
            x0 = x[i]
            xp0 = xp[i]
            x[i] = R11 * x0 + R12 * xp0
            xp[i] = R21*x0 + R22*xp0
        return x, xp

    def rematchXPlanePeakISlice(self, beta=False, alpha=False, nEmit=False):
        peakIPosition = self.slice_max_peak_current_slice
        xslice = self.slice_data(self.x)[peakIPosition]
        xpslice = self.slice_data(self.xp)[peakIPosition]
        x, xp = self.performTransformationPeakISlice(xslice, xpslice, self.x, self.xp, beta, alpha, nEmit)
        self.beam['x'] = x
        self.beam['xp'] = xp

        cpz = self.cp / np.sqrt(self.beam['xp']**2 + self.yp**2 + 1)
        cpx = self.beam['xp'] * cpz
        cpy = self.yp * cpz
        self.beam['px'] = cpx * self.q_over_c
        self.beam['py'] = cpy * self.q_over_c
        self.beam['pz'] = cpz * self.q_over_c

    def rematchYPlanePeakISlice(self, beta=False, alpha=False, nEmit=False):
        peakIPosition = self.slice_max_peak_current_slice
        yslice = self.slice_data(self.y)[peakIPosition]
        ypslice = self.slice_data(self.yp)[peakIPosition]
        y, yp = self.performTransformationPeakISlice(yslice, ypslice, self.y, self.yp, beta, alpha, nEmit)
        self.beam['y'] = y
        self.beam['yp'] = yp

        cpz = self.cp / np.sqrt(self.xp**2 + self.beam['yp']**2 + 1)
        cpx = self.xp * cpz
        cpy = self.beam['yp'] * cpz
        self.beam['px'] = cpx * self.q_over_c
        self.beam['py'] = cpy * self.q_over_c
        self.beam['pz'] = cpz * self.q_over_c


    @property
    def Sx(self):
        return np.sqrt(self.covariance(self.x,self.x))
    @property
    def Sy(self):
        return np.sqrt(self.covariance(self.y,self.y))

    @property
    def Sz(self):
        return np.sqrt(self.covariance(self.z,self.z))

    # These are also useful for some calculations
    @property
    def Mx(self):
        return np.mean(self.beam['x'])

    @property
    def My(self):
        return np.mean(self.beam['y'])
    @property
    def Mz(self):
        return np.mean(self.beam['z'])

    @property
    def Mzn(self):
        return 0.0


    # Added to use with dwa_beam_tools
    @property
    def nMacros(self):
        return len(self.x)

    # Could cause a problem with a variable weight macro file...
    @property
    def charge_per_macro(self):
        return self.charge / self.nMacros

    @property
    def code(self):
        if 'code' in self.beam:
            return self.beam['code']
        else:
            return 'unknown'

    def driftBeam(self, L):
        self.beam['x'] = self.x + L * self.xp
        self.beam['y'] = self.y + L * self.yp
        self.beam['z'] = self.z + L
        # xp,yp,pz doesn't change