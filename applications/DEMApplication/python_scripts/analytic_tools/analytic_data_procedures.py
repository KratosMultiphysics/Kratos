import os
import numpy as np
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from datetime import datetime
import h5py

class SurfaceAnalyzer:
    def __init__(self, smp):
        self.inlet = None
        self.n_particles_accumulated = 0
        self.mass_accumulated = 0.0
        self.radii_accumulated = []
        self.smp_name = smp.Name
        self.face_watcher = AnalyticFaceWatcher(smp)

    def MakeMeasurements(self):
        '''
        From python to c++'''
        self.face_watcher.MakeMeasurements()

    def MakeReading(self):
        times, number_flux, mass_flux, radii_flux, vel_nr_mass, vel_tg_mass = [], [], [], [], [], []
        self.face_watcher.GetTotalFlux(times, number_flux, mass_flux, radii_flux, vel_nr_mass, vel_tg_mass)
        lists = [times, number_flux, mass_flux, radii_flux, vel_nr_mass, vel_tg_mass]
        times, number_flux, mass_flux, radii_flux, vel_nr_mass, vel_tg_mass = [l for l in lists]
        length = len(times)
        assert length == len(number_flux) == len(mass_flux)
        shape = (length, )
        number_flux, mass_flux = self.CalculateAccumulatedVectors(length, number_flux, mass_flux)
        radii_flux = self.CalculateAccumulatedLists(radii_flux)
        if self.inlet is not None:
            self.inlet_accumulated_mass.append(self.inlet.GetMassInjectedSoFar())
            self.inlet_accumulated_number_of_particles.append(self.inlet.GetNumberOfParticlesInjectedSoFar())

        return shape, times, number_flux, mass_flux, radii_flux, vel_nr_mass, vel_tg_mass

    def CalculateAccumulatedVectors(self, length, number_flux, mass_flux):
        acc_number_flux = self.CalculateAccumulated(original_list = number_flux, old_accumulated = self.n_particles_accumulated)
        acc_mass_flux = self.CalculateAccumulated(original_list = mass_flux, old_accumulated = self.mass_accumulated)

        return acc_number_flux, acc_mass_flux

    def CalculateAccumulatedLists(self, radii_flux):
        acc_radii_flux = radii_flux

        return acc_radii_flux

    def CalculateAccumulated(self, original_list, old_accumulated = 0):
        new_accumulated = np.cumsum(np.array(original_list)) + old_accumulated
        return new_accumulated

    def UpdateData(self, time):
        shape, time, n_particles, mass, radii, vel_nr_mass = self.MakeReading()[:-1]
        # initial with 1 for each surface, should be one for each condition in each surface
        total_mass = sum(mass)
        if total_mass:
            avg_vel_nr = vel_nr_mass / total_mass
            # sum (normal vel * particle_mass) / total mass flux of that timestep
        else:
            avg_vel_nr = [0.] * len(mass)
        return shape, time, n_particles, mass, radii, avg_vel_nr


    def MakeInletMassPlot(self):
        self.MakeInletReading()

    def SetInlet(self, inlet):
        self.inlet = inlet

    def UpdateVariables(self, n_particles_old, n_mass_old):
        self.n_particles_accumulated = n_particles_old
        self.mass_accumulated = n_mass_old

    def ClearData(self):
        self.face_watcher.ClearData()

class ParticlesAnalyzerClass:
    def __init__(self, model_part):
        self.analytic_model_part = model_part
        self.particle_analyzer = AnalyticParticleWatcher()

    def MakeAnalyticsMeasurements(self):
       self.particle_analyzer.MakeMeasurements(self.analytic_model_part)

    def SetNodalMaxImpactVelocities(self):
        self.particle_analyzer.SetNodalMaxImpactVelocities(self.analytic_model_part)

    def SetNodalMaxFaceImpactVelocities(self):
        self.particle_analyzer.SetNodalMaxFaceImpactVelocities(self.analytic_model_part)


class SurfacesAnalyzerClass:
    def __init__(self, sub_model_parts, main_path, do_clear_data = True):
        self.sub_model_parts = sub_model_parts
        self.main_path = main_path

        self.surface_analyzers_list = []
        self.do_clear_data = do_clear_data
        self.times_data_base_names = []
        self.n_particles_data_base_names = []
        self.mass_data_base_names = []
        self.new_path = self.main_path + '/flux_data_new.hdf5'
        self.old_path = self.new_path.replace('_new.hdf5', '.hdf5')
        self.name_n_particles = 'n_accum'
        self.name_mass = 'm_accum'
        self.name_radii = 'particles_radii'
        self.name_avg_vel_nr = 'mass_avg_normal_vel'
        self.ghost_smp_detected = False

        for smp in sub_model_parts:
            if smp[IS_GHOST] == True:
                self.ghost_smp_detected = True
                self.surface_analyzers_list.append(SurfaceAnalyzer(smp))

        self.RemoveFiles()


    def MakeAnalyticsMeasurements(self):
        '''
        This function is used as interface to reach the SurfaceAnalyzer MakeMeasurements
        and from there, the cpp function.'''
        for analyzer in self.surface_analyzers_list:
            analyzer.MakeMeasurements()

    def MakeAnalyticsPipeLine(self, time):
        if self.ghost_smp_detected:
            self.CreateNewFile()
            self.UpdateDataBases(time)
            self.RemoveOldFile()

    def CreateNewFile(self):
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
        with h5py.File(self.new_path, 'a') as f, h5py.File(self.old_path, 'r') as f_old:
            for analyzer in self.surface_analyzers_list:
                shape, time, n_particles, mass, radii, avg_vel_nr = analyzer.UpdateData(time)
                shape_old = f_old['/' + analyzer.smp_name + '/time'].shape
                current_shape = (shape_old[0] + shape[0], )
                time_db, n_particles_db, mass_db, avg_vel_nr_db = self.CreateDataSets(f, current_shape, analyzer.smp_name)
                radii_db = self.CreateDataSetOfRadii(f, current_shape[0], analyzer.smp_name)
                surface_data = f.require_group(analyzer.smp_name)
                surface_data.attrs['Surface Identifier'] = analyzer.smp_name

                time_db[:shape_old[0]] = f_old['/' + analyzer.smp_name + '/time'][:]
                time_db[shape_old[0]:] = time
                n_particles_db[:shape_old[0]] = f_old['/' + analyzer.smp_name + '/' + self.name_n_particles][:]
                n_particles_db[shape_old[0]:] = n_particles
                mass_db[:shape_old[0]] = f_old['/' + analyzer.smp_name + '/' + self.name_mass][:]
                mass_db[shape_old[0]:] = mass

                start_time = datetime.now()
                radii_db[:shape_old[0]] = f_old['/' + analyzer.smp_name + '/' + self.name_radii][:]
                end_time = datetime.now()
                print('Duration: {}'.format(end_time - start_time))


                radii_db[shape_old[0]:,0] = radii

                avg_vel_nr_db[:shape_old[0]] = f_old['/' + analyzer.smp_name + '/' + self.name_avg_vel_nr][:]
                avg_vel_nr_db[shape_old[0]:] = avg_vel_nr
                if self.do_clear_data:
                    if len(n_particles):
                        (analyzer).UpdateVariables(n_particles[-1], mass[-1])
                    (analyzer).ClearData()

        # how to extract data from h5 subgrouped datasets:
        #input_data = h5py.File('Cemib_P660_SpreadPattern.dat.hdf5','r')
        #x_h5 = input_data.get('/patch/X')
        #x = np.array(x_h5)

    def CreateDataFile(self, time):

        # how to create subgrouped datasets with variable name:
        # group2 = f.create_group('group2/subfolder')
        # group2.create_dataset('data',data=d)
        with h5py.File(self.new_path, 'a') as f:
            for analyzer in self.surface_analyzers_list:
                    shape, time, n_particles, mass, radii, avg_vel_nr = analyzer.UpdateData(time)
                    time_db, n_particles_db, mass_db, avg_vel_nr_db = self.CreateDataSets(f, shape, analyzer.smp_name)
                    radii_db = self.CreateDataSetOfRadii(f, len(radii), analyzer.smp_name)
                    time_db[:] = time[:]
                    n_particles_db[:] = n_particles[:]
                    mass_db[:] = mass[:]
                    avg_vel_nr_db[:] = avg_vel_nr[:]
                    radii_db[:,0] = radii[:]
                    if self.do_clear_data:
                        if len(n_particles):
                            analyzer.UpdateVariables(n_particles[-1], mass[-1])
                        analyzer.ClearData()

    def CreateDataSets(self, f, current_shape, sp_name):
        surface_data = f.require_group(sp_name)
        surface_data.attrs['Surface Identifier'] = sp_name
        time_db = surface_data.require_dataset('time', shape = current_shape, dtype = np.float64)
        n_particles_db = surface_data.require_dataset(self.name_n_particles, shape = current_shape, dtype = np.int64)
        mass_db = surface_data.require_dataset(self.name_mass, shape = current_shape, dtype = np.float64)
        avg_vel_nr_db = surface_data.require_dataset(self.name_avg_vel_nr, shape = current_shape, dtype = np.float64)

        return time_db, n_particles_db, mass_db, avg_vel_nr_db

    def CreateDataSetOfRadii(self, f, shape, sp_name):
        dt = h5py.special_dtype(vlen=np.dtype('float64'))
        surface_data = f.require_group(sp_name)
        radii_db = surface_data.create_dataset(self.name_radii, shape=(shape,1) , maxshape = (None,1), dtype = dt)
        return radii_db


    def OldFileExists(self):
        return os.path.exists(self.old_path)

    #TODO: Decide what to do with these unused.
    # Currently not being used
    def GetJointData(self, data_base_names):
        data_list = []
        with h5py.File(self.new_path, 'r') as f:
            if self.do_clear_data: # join all databases
                for name in data_base_names:
                    data_list.append(f['/' + name].value)
                joint_list = np.concatenate(data_list, axis = 0)
            else: # get the latest
                joint_list = f['/' + data_base_names[-1]].value

        return joint_list

    # Currently not being used
    def GetMassFlux(self):
        return self.GetJointData(self.mass_data_base_names)

    # Currently not being used
    def MakeTotalFluxPlot(self):
        import matplotlib.pyplot as plt
        with h5py.File(self.file_path) as f:
            times = f['/' + self.face_watcher_name + '/' + '/time'].value
            mass_flux = f['/' + self.face_watcher_name + '/' + '/m_accum'].value
        plt.xlabel('time')
        plt.ylabel('accumulated mass throughput')
        plt.plot(times, mass_flux)
        plt.savefig(self.main_path + '/mass_throughput.pdf', bbox_inches='tight')

    # Currently not being used
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
    # Currently not being used
    def GetNumberOfParticlesFlux(self):
        return self.GetJointData(self.n_particles_data_base_names)

    # Currently not being used
    def GetTimes(self):
        return self.GetJointData(self.times_data_base_names)
