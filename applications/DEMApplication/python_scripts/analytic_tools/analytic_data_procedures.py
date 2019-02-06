import os

class ParticleWatcherAnalyzer:
    def __init__(self, analytic_particle_watcher, path, do_clear_data = True):
        self.particle_watcher = analytic_particle_watcher

    # call everytime post results are generated
    def SetNodalMaxImpactVelocities(self, analytic_model_part):
        self.analytic_particle_watcher.SetNodalMaxImpactVelocities(analytic_model_part)

    def SetNodalMaxFaceImpactVelocities(self, analytic_model_part):
        self.analytic_particle_watcher.SetNodalMaxFaceImpactVelocities(analytic_model_part)


class FaceWatcherAnalyzer:
    def __init__(self, name, analytic_face_watcher, path, do_clear_data = True):
        self.face_watcher = analytic_face_watcher
        self.face_watcher_name = name
        self.do_clear_data = do_clear_data
        # the following objects are useful if data is chunked into several databases
        self.times_data_base_names = []
        self.n_particles_data_base_names = []
        self.mass_data_base_names = []
        self.n_particles_accumulated = 0
        self.mass_accumulated = 0.0
        self.inlet = None

        self.folder_path = path
        self.file_path = path + '/flux_data_new.hdf5'
        self.old_file_path = self.file_path.replace('_new.hdf5', '.hdf5')

        FaceWatcherAnalyzer.file_path = self.file_path
        FaceWatcherAnalyzer.file_path_old =self.old_file_path
        FaceWatcherAnalyzer.RemoveFiles()

    @staticmethod
    def RemoveFiles():
        new_path = FaceWatcherAnalyzer.file_path
        old_path = FaceWatcherAnalyzer.file_path_old
        for path in (p for p in [new_path, old_path] if os.path.exists(p)):
            os.remove(path)

    @staticmethod
    def CreateNewFile():
        new_path = FaceWatcherAnalyzer.file_path
        old_path = FaceWatcherAnalyzer.file_path_old
        if os.path.exists(new_path):
            os.rename(new_path, old_path)

        import h5py
        h5py.File(new_path)

    @staticmethod
    def RemoveOldFile():
        old_path = FaceWatcherAnalyzer.file_path_old
        if os.path.exists(old_path):
            os.remove(old_path)

    def MakeReading(self):
        times, number_flux, mass_flux, vel_nr_mass, vel_tg_mass = [], [], [], [], []
        self.face_watcher.GetTotalFlux(times, number_flux, mass_flux, vel_nr_mass, vel_tg_mass)
        lists = [times, number_flux, mass_flux, vel_nr_mass, vel_tg_mass]
        #times, number_flux, mass_flux, vel_nr_mass, vel_tg_mass = [np.array(l) for l in lists]
        times, number_flux, mass_flux, vel_nr_mass, vel_tg_mass = [l for l in lists]
        length = len(times)
        assert length == len(number_flux) == len(mass_flux)
        shape = (length, )
        number_flux, mass_flux = self.CalculateAccumulatedVectors(length, number_flux, mass_flux)
        if self.inlet is not None:
            self.inlet_accumulated_mass.append(self.inlet.GetMassInjectedSoFar())
            self.inlet_accumulated_number_of_particles.append(self.inlet.GetNumberOfParticlesInjectedSoFar())

        return shape, times, number_flux, mass_flux, vel_nr_mass, vel_tg_mass

    def GetTimes(self):
        return self.GetJointData(self.times_data_base_names)

    def GetNumberOfParticlesFlux(self):
        return self.GetJointData(self.n_particles_data_base_names)

    def GetMassFlux(self):
        return self.GetJointData(self.mass_data_base_names)

    def GetJointData(self, data_base_names):
        import numpy as np
        import h5py
        data_list = []

        with h5py.File(self.file_path, 'r') as f:
            if self.do_clear_data: # join all databases
                for name in data_base_names:
                    data_list.append(f['/' + name].value)
                joint_list = np.concatenate(data_list, axis = 0)
            else: # get the latest
                joint_list = f['/' + data_base_names[-1]].value

        return joint_list

    def CalculateAccumulatedVectors(self, length, number_flux, mass_flux):
        acc_number_flux = self.CalculateAccumulated(original_list = number_flux, old_accumulated = self.n_particles_accumulated)
        acc_mass_flux = self.CalculateAccumulated(original_list = mass_flux, old_accumulated = self.mass_accumulated)

        return acc_number_flux, acc_mass_flux

    def CalculateAccumulated(self, original_list, old_accumulated = 0):
        import numpy as np
        new_accumulated = np.cumsum(np.array(original_list)) + old_accumulated
        return new_accumulated

    def OldFileExists(self):
        return os.path.exists(self.old_file_path)

    def UpdateDataFiles(self, time):
        shape, time, n_particles, mass, vel_nr_mass = self.MakeReading()[:-1]  # initial with 1 for each surface, should be one for each condition in each surface
        total_mass = sum(mass)
        #total_mass = np.sum(mass)
        if total_mass:
            avg_vel_nr = vel_nr_mass / total_mass  # sum (normal vel * particle_mass) / total mass flux of that timestep
            #avg_vel_tg = vel_tg_mass / total_mass
        else:
            #avg_vel_nr = np.zeros(mass.size)
            avg_vel_nr = [0.] * len(mass)
            #avg_vel_tg = np.zeros(mass.size)
        name_n_particles = 'n_accum'
        name_mass = 'm_accum'
        name_avg_vel_nr = 'mass_avg_normal_vel'

        # how to create subgrouped datasets with variable name:
        # group2 = f.create_group('group2/subfolder')
        # group2.create_dataset('data',data=d)

        def CreateDataSets(f, current_shape):
            surface_data = f.require_group(self.face_watcher_name)
            surface_data.attrs['Surface Identifier'] = self.face_watcher_name

            #time_db = surface_data.require_dataset('time', shape = current_shape, dtype = np.float64)
            #n_particles_db = surface_data.require_dataset(name_n_particles, shape = current_shape, dtype = np.int64)
            #mass_db = surface_data.require_dataset(name_mass, shape = current_shape, dtype = np.float64)
            #avg_vel_nr_db = surface_data.require_dataset(name_avg_vel_nr, shape = current_shape, dtype = np.float64)

            time_db = surface_data.require_dataset('time', shape = current_shape, dtype = float)
            n_particles_db = surface_data.require_dataset(name_n_particles, shape = current_shape, dtype = int)
            mass_db = surface_data.require_dataset(name_mass, shape = current_shape, dtype = float)
            avg_vel_nr_db = surface_data.require_dataset(name_avg_vel_nr, shape = current_shape, dtype = float)
            return time_db, n_particles_db, mass_db, avg_vel_nr_db

        import h5py
        if self.OldFileExists():

            with h5py.File(self.file_path) as f, h5py.File(self.old_file_path, 'r') as f_old:
                shape_old = f_old['/' + self.face_watcher_name + '/time'].shape
                current_shape = (shape_old[0] + shape[0], )
                time_db, n_particles_db, mass_db, avg_vel_nr_db = CreateDataSets(f, current_shape)

                time_db[:shape_old[0]] = f_old['/' + self.face_watcher_name + '/time'][:]
                time_db[shape_old[0]:] = time[:]
                n_particles_db[:shape_old[0]] = f_old['/' + self.face_watcher_name + '/' + name_n_particles][:]
                n_particles_db[shape_old[0]:] = n_particles[:]
                mass_db[:shape_old[0]] = f_old['/' + self.face_watcher_name + '/' + name_mass][:]
                mass_db[shape_old[0]:] = mass[:]
                avg_vel_nr_db[:shape_old[0]] = f_old['/' + self.face_watcher_name + '/' + name_avg_vel_nr][:]
                avg_vel_nr_db[shape_old[0]:] = avg_vel_nr[:]

        else:
            with h5py.File(self.file_path) as f:
                time_db, n_particles_db, mass_db, avg_vel_nr_db = CreateDataSets(f, shape)
                time_db[:] = time[:]
                n_particles_db[:] = n_particles[:]
                mass_db[:] = mass[:]
                avg_vel_nr_db[:] = avg_vel_nr[:]

        # how to extract data from h5 subgrouped datasets:
        #input_data = h5py.File('Cemib_P660_SpreadPattern.dat.hdf5','r')
        #x_h5 = input_data.get('/patch/X')
        #x = np.array(x_h5)

        if self.do_clear_data:
            if len(n_particles):
                self.n_particles_accumulated = n_particles[-1]
                self.mass_accumulated = mass[-1]
            self.face_watcher.ClearData()

    def MakeTotalFluxPlot(self):
        import matplotlib.pyplot as plt
        import h5py
        with h5py.File(self.file_path) as f:
            times = f['/' + self.face_watcher_name + '/' + '/time'].value
            mass_flux = f['/' + self.face_watcher_name + '/' + '/m_accum'].value
        plt.xlabel('time')
        plt.ylabel('accumulated mass throughput')
        plt.plot(times, mass_flux)
        plt.savefig(self.folder_path + '/mass_throughput.pdf', bbox_inches='tight')

    def MakeFluxOfNumberOfParticlesPlot(self):
        import matplotlib.pyplot as plt
        self.MakeReading()
        times = self.GetTimes()
        flux = self.GetNumberOfParticlesFlux()
        '''
        plt.xlabel('time')
        plt.ylabel('accumulated number of particles through surface')
        plt.plot(times, flux)
        plt.savefig(self.folder_path + '/throughput.svg')
        plt.clf()
        '''

    def MakeInletMassPlot(self):
        self.MakeInletReading()

    def SetInlet(self, inlet):
        self.inlet = inlet
