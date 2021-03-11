import os
import numpy as np
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

class ParticleWatcherAnalyzer:
    def __init__(self, analytic_particle_watcher, path, do_clear_data = True):
        self.particle_watcher = analytic_particle_watcher

    # call everytime post results are generated
    def SetNodalMaxImpactVelocities(self, analytic_model_part):
        self.analytic_particle_watcher.SetNodalMaxImpactVelocities(analytic_model_part)

    def SetNodalMaxFaceImpactVelocities(self, analytic_model_part):
        self.analytic_particle_watcher.SetNodalMaxFaceImpactVelocities(analytic_model_part)


class FaceWatcherAnalyzer:
    def __init__(self, name, analytic_face_watcher):
        self.face_watcher = analytic_face_watcher
        self.face_watcher_name = name
        # the following objects are useful if data is chunked into several databases
        self.inlet = None
        self.n_particles_accumulated = 0
        self.mass_accumulated = 0.0

    def MakeReading(self):
        times, number_flux, mass_flux, vel_nr_mass, vel_tg_mass = [], [], [], [], []
        self.face_watcher.GetTotalFlux(times, number_flux, mass_flux, vel_nr_mass, vel_tg_mass)
        lists = [times, number_flux, mass_flux, vel_nr_mass, vel_tg_mass]
        times, number_flux, mass_flux, vel_nr_mass, vel_tg_mass = [l for l in lists]
        length = len(times)
        assert length == len(number_flux) == len(mass_flux)
        shape = (length, )
        number_flux, mass_flux = self.CalculateAccumulatedVectors(length, number_flux, mass_flux)
        if self.inlet is not None:
            self.inlet_accumulated_mass.append(self.inlet.GetMassInjectedSoFar())
            self.inlet_accumulated_number_of_particles.append(self.inlet.GetNumberOfParticlesInjectedSoFar())

        return shape, times, number_flux, mass_flux, vel_nr_mass, vel_tg_mass

    def CalculateAccumulatedVectors(self, length, number_flux, mass_flux):
        acc_number_flux = self.CalculateAccumulated(original_list = number_flux, old_accumulated = self.n_particles_accumulated)
        acc_mass_flux = self.CalculateAccumulated(original_list = mass_flux, old_accumulated = self.mass_accumulated)

        return acc_number_flux, acc_mass_flux

    def CalculateAccumulated(self, original_list, old_accumulated = 0):
        new_accumulated = np.cumsum(np.array(original_list)) + old_accumulated
        return new_accumulated

    def UpdateData(self, time):
        shape, time, n_particles, mass, vel_nr_mass = self.MakeReading()[:-1]  # initial with 1 for each surface, should be one for each condition in each surface
        total_mass = sum(mass)
        if total_mass:
            avg_vel_nr = vel_nr_mass / total_mass  # sum (normal vel * particle_mass) / total mass flux of that timestep
        else:
            avg_vel_nr = [0.] * len(mass)
        return shape, time, n_particles, mass, avg_vel_nr

    def MakeInletMassPlot(self):
        self.MakeInletReading()

    def SetInlet(self, inlet):
        self.inlet = inlet

    def UpdateVariables(self, n_particles_old, n_mass_old):
        self.n_particles_accumulated = n_particles_old
        self.mass_accumulated = n_mass_old

    def ClearData(self):
        self.face_watcher.ClearData()

class FaceAnalyzerClass:
    def __init__(self, model_part, main_path, do_clear_data = True):
        self.model_part = model_part
        self.main_path = main_path
        self.face_watcher_dict = dict()
        self.face_watcher_analysers = dict()
        self.do_clear_data = do_clear_data
        self.times_data_base_names = []
        self.n_particles_data_base_names = []
        self.mass_data_base_names = []
        self.new_path = self.main_path + '/flux_data_new.hdf5'
        self.old_path = self.new_path.replace('_new.hdf5', '.hdf5')
        self.name_n_particles = 'n_accum'
        self.name_mass = 'm_accum'
        self.name_avg_vel_nr = 'mass_avg_normal_vel'

        for sub_part in self.model_part:
            if sub_part[IS_GHOST] == True:
                self.face_watcher_dict[sub_part.Name] = AnalyticFaceWatcher(sub_part)
                self.face_watcher_analysers[sub_part.Name] = FaceWatcherAnalyzer(name=sub_part.Name, analytic_face_watcher=self.face_watcher_dict[sub_part.Name])

        self.RemoveFiles()

    def MakeAnalyticsMeasurements(self):
        for face_watcher in self.face_watcher_dict.values():
            face_watcher.MakeMeasurements()

    def MakeAnalyticsPipeLine(self, time):
        self.CreateNewFile()
        self.UpdateDataBases(time)
        self.RemoveOldFile()

    def CreateNewFile(self):
        import h5py
        if os.path.exists(self.new_path):
            os.rename(self.new_path, self.old_path)

        h5py.File(self.new_path, 'a')

    def RemoveOldFile(self):
        if os.path.exists(self.old_path):
            os.remove(self.old_path)

    def RemoveFiles(self):
        for path in (p for p in [self.new_path, self.old_path] if os.path.exists(p)):
            os.remove(path)

    def UpdateDataBases(self, time):
        if self.OldFileExists():
            self.UpdateDataFile(time)
        else:
            self.CreateDataFile(time)

    def UpdateDataFile(self, time):

        # how to create subgrouped datasets with variable name:
        # group2 = f.create_group('group2/subfolder')
        # group2.create_dataset('data',data=d)
        import h5py
        with h5py.File(self.new_path, 'a') as f, h5py.File(self.old_path, 'r') as f_old:
            for sub_part in self.model_part:
                if sub_part[IS_GHOST]:
                    shape, time, n_particles, mass, avg_vel_nr = self.face_watcher_analysers[sub_part.Name].UpdateData(time)
                    shape_old = f_old['/' + sub_part.Name + '/time'].shape
                    current_shape = (shape_old[0] + shape[0], )
                    time_db, n_particles_db, mass_db, avg_vel_nr_db = self.CreateDataSets(f, current_shape, sub_part.Name)

                    time_db[:shape_old[0]] = f_old['/' + sub_part.Name + '/time'][:]
                    time_db[shape_old[0]:] = time[:]
                    n_particles_db[:shape_old[0]] = f_old['/' + sub_part.Name + '/' + self.name_n_particles][:]
                    n_particles_db[shape_old[0]:] = n_particles[:]
                    mass_db[:shape_old[0]] = f_old['/' + sub_part.Name + '/' + self.name_mass][:]
                    mass_db[shape_old[0]:] = mass[:]
                    avg_vel_nr_db[:shape_old[0]] = f_old['/' + sub_part.Name + '/' + self.name_avg_vel_nr][:]
                    avg_vel_nr_db[shape_old[0]:] = avg_vel_nr[:]
                    if self.do_clear_data:
                        if len(n_particles):
                            self.face_watcher_analysers[sub_part.Name].UpdateVariables(n_particles[-1], mass[-1])
                        self.face_watcher_analysers[sub_part.Name].ClearData()

        # how to extract data from h5 subgrouped datasets:
        #input_data = h5py.File('Cemib_P660_SpreadPattern.dat.hdf5','r')
        #x_h5 = input_data.get('/patch/X')
        #x = np.array(x_h5)

    def CreateDataFile(self, time):

        # how to create subgrouped datasets with variable name:
        # group2 = f.create_group('group2/subfolder')
        # group2.create_dataset('data',data=d)
        import h5py
        with h5py.File(self.new_path, 'a') as f:
            for sub_part in self.model_part:
                if sub_part[IS_GHOST]:
                    shape, time, n_particles, mass, avg_vel_nr = self.face_watcher_analysers[sub_part.Name].UpdateData(time)
                    time_db, n_particles_db, mass_db, avg_vel_nr_db = self.CreateDataSets(f, shape, sub_part.Name)
                    time_db[:] = time[:]
                    n_particles_db[:] = n_particles[:]
                    mass_db[:] = mass[:]
                    avg_vel_nr_db[:] = avg_vel_nr[:]
                    if self.do_clear_data:
                        if len(n_particles):
                            self.face_watcher_analysers[sub_part.Name].UpdateVariables(n_particles[-1], mass[-1])
                        self.face_watcher_analysers[sub_part.Name].ClearData()

        #input_data = h5py.File('Cemib_P660_SpreadPattern.dat.hdf5','r')
        #x_h5 = input_data.get('/patch/X')
        #x = np.array(x_h5)

    def CreateDataSets(self, f, current_shape, sp_name):
        surface_data = f.require_group(sp_name)
        surface_data.attrs['Surface Identifier'] = sp_name

        time_db = surface_data.require_dataset('time', shape = current_shape, dtype = np.float64)
        n_particles_db = surface_data.require_dataset(self.name_n_particles, shape = current_shape, dtype = np.int64)
        mass_db = surface_data.require_dataset(self.name_mass, shape = current_shape, dtype = np.float64)
        avg_vel_nr_db = surface_data.require_dataset(self.name_avg_vel_nr, shape = current_shape, dtype = np.float64)

        return time_db, n_particles_db, mass_db, avg_vel_nr_db

    def OldFileExists(self):
        return os.path.exists(self.old_path)

    def GetJointData(self, data_base_names):
        data_list = []
        import h5py
        with h5py.File(self.new_path, 'r') as f:
            if self.do_clear_data: # join all databases
                for name in data_base_names:
                    data_list.append(f['/' + name].value)
                joint_list = np.concatenate(data_list, axis = 0)
            else: # get the latest
                joint_list = f['/' + data_base_names[-1]].value

        return joint_list

    def GetTimes(self):
        return self.GetJointData(self.times_data_base_names)

    def GetNumberOfParticlesFlux(self):
        return self.GetJointData(self.n_particles_data_base_names)

    def GetMassFlux(self):
        return self.GetJointData(self.mass_data_base_names)

    def MakeTotalFluxPlot(self):
        import matplotlib.pyplot as plt
        import h5py
        with h5py.File(self.file_path) as f:
            times = f['/' + self.face_watcher_name + '/' + '/time'].value
            mass_flux = f['/' + self.face_watcher_name + '/' + '/m_accum'].value
        plt.xlabel('time')
        plt.ylabel('accumulated mass throughput')
        plt.plot(times, mass_flux)
        plt.savefig(self.main_path + '/mass_throughput.pdf', bbox_inches='tight')

    def MakeFluxOfNumberOfParticlesPlot(self):
        import matplotlib.pyplot as plt
        self.MakeReading()
        times = self.GetTimes()
        flux = self.GetNumberOfParticlesFlux()
        '''
        plt.xlabel('time')
        plt.ylabel('accumulated number of particles through surface')
        plt.plot(times, flux)
        plt.savefig(self.main_path + '/throughput.svg')
        plt.clf()
        '''
