import matplotlib.pyplot as plt
import h5py
import numpy as np
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
    def __init__(self, analytic_face_watcher, path, do_clear_data = True):
        self.face_watcher = analytic_face_watcher
        self.dtype = np.float64
        self.do_clear_data = do_clear_data
        # the following objects are useful if data is chunked into several databases
        self.times_data_base_names = []
        self.n_particles_data_base_names = []
        self.mass_data_base_names = []
        self.n_particles_accumulated = 0
        self.mass_accumulated = 0.0
        self.folder_path = path
        self.file_path = path + '/flux_data.hdf5'
        self.inlet = None

        if os.path.exists(self.file_path):
            os.remove(self.file_path)

    def MakeReading(self):
        times, number_flux, mass_flux = [], [], []
        self.face_watcher.GetTotalFlux(times, number_flux, mass_flux)
        length = len(times)
        assert length == len(number_flux) == len(mass_flux)
        shape = (length, )
        number_flux, mass_flux = self.CalculateAccumulatedVectors(length, number_flux, mass_flux)
        if self.inlet is not None:
            self.inlet_accumulated_mass.append(self.inlet.GetMassInjectedSoFar())
            self.inlet_accumulated_number_of_particles.append(self.inlet.GetNumberOfParticlesInjectedSoFar())

        return shape, times, number_flux, mass_flux

    def GetTimes(self):
        return self.GetJointData(self.times_data_base_names)

    def GetNumberOfParticlesFlux(self):
        return self.GetJointData(self.n_particles_data_base_names)

    def GetMassFlux(self):
        return self.GetJointData(self.mass_data_base_names)

    def GetJointData(self, data_base_names):
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
        new_accumulated = np.cumsum(np.array(original_list)) + old_accumulated
        return new_accumulated

    def UpdateDataFile(self, time):
        shape, times, n_particles_data, mass_data = self.MakeReading()
        label = str(len(self.times_data_base_names))
        name_times = 'times ' + label
        name_n_particles = 'accumulated throughput (number of particles) ' + label
        name_mass = 'accumulated mass throughput ' + label
        self.times_data_base_names.append(name_times)
        self.n_particles_data_base_names.append(name_n_particles)
        self.mass_data_base_names.append(name_mass)

        with h5py.File(self.file_path) as f:
            f.require_dataset(name_times, data = times, shape = shape, dtype = np.float64)
            f.require_dataset(name_n_particles, data = n_particles_data, shape = shape, dtype = np.int64)
            f.require_dataset(name_mass, data = mass_data, shape = shape, dtype = np.float64)

        if self.do_clear_data:
            if len(n_particles_data):
                self.n_particles_accumulated = n_particles_data[-1]
                self.mass_accumulated = mass_data[-1]
            self.face_watcher.ClearData()

    def MakeTotalFluxPlot(self):
        self.MakeReading()
        times = self.GetTimes()
        mass_flux = self.GetMassFlux()
        plt.xlabel('time')
        plt.ylabel('accumulated mass throughput')
        plt.plot(times, mass_flux)
        plt.savefig(self.folder_path + '/mass_throughput.svg')

    def MakeFluxOfNumberOfParticlesPlot(self):
        self.MakeReading()
        times = self.GetTimes()
        flux = self.GetNumberOfParticlesFlux()
        plt.xlabel('time')
        plt.ylabel('accumulated number of particles through surface')
        plt.plot(times, flux)
        plt.savefig(self.folder_path + '/throughput.svg')
        plt.clf()

    def MakeInletMassPlot(self):
        self.MakeInletReading()

    def SetInlet(self, inlet):
        self.inlet = inlet
