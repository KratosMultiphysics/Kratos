import bisect as bi
import numpy as np
import h5py
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *
from DEM_procedures import KratosPrint as Say
import json

def DeleteDataSet(file_or_group, dset_name):
    if dset_name in file_or_group:
        file_or_group.__delitem__(dset_name)

def CreateDataset(file_or_group, name, data):
    if name in file_or_group:
        file_or_group.__delitem__(name)

    dtype = float

    if len(data):
        dtype = type(data[0])

    file_or_group.create_dataset(dtype = dtype, name = name, data = data)

def CreateGroup(file_or_group, name):
    if name in file_or_group:
        file_or_group['/'].__delitem__(name)

    return file_or_group.create_group(name)

def WriteDataToFile(file_or_group, names, data):
    for name, datum in zip(names, data):
        DeleteDataSet(file_or_group, name)
        file_or_group.create_dataset(name = name, data = datum)

def Index():
    index = 0
    while True:
        yield index
        index += 1


class FluidHDF5Loader:

    def __init__(self, fluid_model_part, particles_model_part, pp, main_path):
        self.n_nodes = len(fluid_model_part.Nodes)
        self.shape = (self.n_nodes,)
        self.store_pressure = pp.CFD_DEM["store_fluid_pressure_option"].GetBool()
        self.store_gradient = pp.CFD_DEM["store_full_gradient_option"].GetBool()
        self.there_are_more_steps_to_load = True
        self.main_path = main_path
        self.pp = pp
        self.fluid_model_part = fluid_model_part
        self.disperse_phase_model_part = particles_model_part

        number_of_variables = 3

        if pp.CFD_DEM["store_fluid_pressure_option"].GetBool():
            number_of_variables += 1
        if pp.CFD_DEM["load_derivatives"].GetBool():
            number_of_variables += 9

        self.extended_shape = self.shape + (number_of_variables, )
        self.file_name = self.pp.CFD_DEM.AddEmptyValue("prerun_fluid_file_name").GetString()
        self.file_path = main_path + self.file_name

        if pp.CFD_DEM["fluid_already_calculated"].GetBool():

            with h5py.File(self.file_path, 'r') as f:
                self.times_str = list([str(key) for key in f.keys() if key not in {'nodes'}])
                nodes_ids = np.array([node_id for node_id in f['nodes'][:, 0]])
                self.permutations = np.array(range(len(nodes_ids)))
                # obtaining the vector of permutations by ordering [0, 1, ..., n_nodes] as nodes_ids, by increasing order of id.
                self.permutations = np.array([x for (y, x) in sorted(zip(nodes_ids, self.permutations))])
                self.times     = np.array([float(f[key].attrs['time']) for key in self.times_str])
                self.times_str = np.array([x for (y, x) in sorted(zip(self.times, self.times_str))])
                self.times = sorted(self.times)
                self.dt = self.times[-1] - self.times[-2]

            self.old_data_array = np.zeros(self.extended_shape)
            self.future_data_array = np.zeros(self.extended_shape)
            self.old_time_index = 0
            self.future_time_index = 1
            viscosity = 1e-6
            density = 1000.
            for node in self.fluid_model_part.Nodes:
                node.SetSolutionStepValue(VISCOSITY, viscosity)
                node.SetSolutionStepValue(DENSITY, density)
        else:
            self.compression_type = 'gzip'
            for node in self.fluid_model_part.Nodes:
                viscosity = node.GetSolutionStepValue(VISCOSITY)
                density = node.GetSolutionStepValue(DENSITY)
                break
            with h5py.File(self.file_path, 'w') as f:
                f.attrs['kinematic viscosity'] = viscosity
                f.attrs['time step'] = pp.Dt
                f.attrs['density'] = density
                f.attrs['solver type'] = pp.FluidSolverConfiguration.solver_type
                f.attrs['linear system solver type'] = pp.FluidSolverConfiguration.linear_solver_config.solver_type
                f.attrs['use orthogonal subscales'] = bool(pp.FluidSolverConfiguration.oss_switch)
                self.dtype = np.float64
                nodes = np.array([(node.Id, node.X, node.Y, node.Z) for node in fluid_model_part.Nodes])
                f.create_dataset(name = 'nodes', compression = self.compression_type, data = nodes, dtype = np.float64)

            if pp.CFD_DEM["store_fluid_in_single_precision"].GetBool():
                self.dtype = np.float32

        self.last_time = 0.0

        self.current_data_array = np.zeros(self.extended_shape)

    def GetOldTimeIndicesAndWeights(self, current_time, times_array, fluid_dt):
        old_index = bi.bisect(times_array, current_time)
        future_index = old_index + 1
        old_time =  times_array[old_index]

        if future_index >= len(times_array):
            alpha_old = 1
            future_index = old_index
            alpha_future = 0
            self.there_are_more_steps_to_load = False
        else:
            alpha_old = max(0, (current_time - old_time) / fluid_dt)
            alpha_future = 1.0 - alpha_old

        return old_index, alpha_old, future_index, alpha_future

    def CanLoadMoreSteps(self):
        return self.there_are_more_steps_to_load

    def FillUpSingleDataset(self, name, variable, variable_index_in_temp_array):
        with h5py.File(self.file_path, 'r+') as f:
            f.create_dataset(name, compression = self.compression_type, shape = self.shape, dtype = self.dtype)
            for i_node, node in enumerate(self.fluid_model_part.Nodes):
                self.current_data_array[i_node, variable_index_in_temp_array] = node.GetSolutionStepValue(variable)
            f[name][:] = self.current_data_array[:, variable_index_in_temp_array]

    def FillFluidDataStep(self):
        time = self.fluid_model_part.ProcessInfo[TIME]
        name = str(time)
        with h5py.File(self.file_path) as f:
            f.create_group(name = name)
            f[name].attrs['time'] = time

        indices = Index()
        if not self.last_time == time:
            self.FillUpSingleDataset(name + '/vx', VELOCITY_X, next(indices))
            self.FillUpSingleDataset(name + '/vy', VELOCITY_Y, next(indices))
            self.FillUpSingleDataset(name + '/vz', VELOCITY_Z, next(indices))

            if self.store_pressure:
                self.FillUpSingleDataset(name + '/p', PRESSURE, next(indices))

        if self.store_gradient:
            self.FillUpSingleDataset(name + '/dvxx', VELOCITY_X_GRADIENT_X, next(indices))
            self.FillUpSingleDataset(name + '/dvxy', VELOCITY_X_GRADIENT_Y, next(indices))
            self.FillUpSingleDataset(name + '/dvxz', VELOCITY_X_GRADIENT_Z, next(indices))
            self.FillUpSingleDataset(name + '/dvyx', VELOCITY_Y_GRADIENT_X, next(indices))
            self.FillUpSingleDataset(name + '/dvyy', VELOCITY_Y_GRADIENT_Y, next(indices))
            self.FillUpSingleDataset(name + '/dvyz', VELOCITY_Y_GRADIENT_Z, next(indices))
            self.FillUpSingleDataset(name + '/dvzx', VELOCITY_Z_GRADIENT_X, next(indices))
            self.FillUpSingleDataset(name + '/dvzy', VELOCITY_Z_GRADIENT_Y, next(indices))
            self.FillUpSingleDataset(name + '/dvzz', VELOCITY_Z_GRADIENT_Z, next(indices))

        self.last_time = time

    def ConvertComponent(self, f, component_name):
        if '/vx' in component_name:
            read_values = f[component_name][:,]
        elif '/vy' in component_name:
            read_values = f[component_name][:,]
        elif '/vz' in component_name:
            read_values = f[component_name][:,]
        else:
            read_values = f[component_name][:,]

        return read_values[self.permutations]

    def UpdateFluidVariable(self, name, variable, variable_index_in_temp_array, must_load_future_values_from_database, alpha_old, alpha_future):
        if must_load_future_values_from_database:
            with h5py.File(self.file_path, 'r') as f:
                self.future_data_array[:, variable_index_in_temp_array] = self.ConvertComponent(f, name)

        self.current_data_array[:, variable_index_in_temp_array] = alpha_old * self.old_data_array[:, variable_index_in_temp_array] + alpha_future * self.future_data_array[:, variable_index_in_temp_array]

        for i_node, node in enumerate(self.fluid_model_part.Nodes):
            node.SetSolutionStepValue(variable, self.current_data_array[i_node, variable_index_in_temp_array])

    def LoadFluid(self, DEM_time):
        Say('\nLoading fluid from hdf5 file...')
        # getting time indices and weights (identifyint the two fluid time steps surrounding the current DEM step and assigning correspnding weights)
        old_time_index, alpha_old, future_time_index, alpha_future = self.GetOldTimeIndicesAndWeights(DEM_time, self.times, self.dt)
        old_step_dataset_name    = self.times_str[old_time_index]
        future_step_dataset_name = self.times_str[future_time_index]
        must_load_from_database = not self.old_time_index == old_time_index # old and future time steps must be updated

        if must_load_from_database:
            # new old becomes old future
            self.old_data_array, self.future_data_array = self.future_data_array, self.old_data_array

        indices = Index()
        self.UpdateFluidVariable(future_step_dataset_name + '/vx', VELOCITY_X, next(indices), must_load_from_database, alpha_old, alpha_future)
        self.UpdateFluidVariable(future_step_dataset_name + '/vy', VELOCITY_Y, next(indices), must_load_from_database, alpha_old, alpha_future)
        self.UpdateFluidVariable(future_step_dataset_name + '/vz', VELOCITY_Z, next(indices), must_load_from_database, alpha_old, alpha_future)

        if self.store_pressure:
            self.UpdateFluidVariable(future_step_dataset_name + '/p', PRESSURE, next(indices), must_load_from_database, alpha_old, alpha_future)

        if self.store_gradient:
            self.UpdateFluidVariable(future_step_dataset_name + '/dvxx', VELOCITY_X_GRADIENT_X, next(indices), must_load_from_database, alpha_old, alpha_future)
            self.UpdateFluidVariable(future_step_dataset_name + '/dvxy', VELOCITY_X_GRADIENT_Y, next(indices), must_load_from_database, alpha_old, alpha_future)
            self.UpdateFluidVariable(future_step_dataset_name + '/dvxz', VELOCITY_X_GRADIENT_Z, next(indices), must_load_from_database, alpha_old, alpha_future)
            self.UpdateFluidVariable(future_step_dataset_name + '/dvyx', VELOCITY_Y_GRADIENT_X, next(indices), must_load_from_database, alpha_old, alpha_future)
            self.UpdateFluidVariable(future_step_dataset_name + '/dvyy', VELOCITY_Y_GRADIENT_Y, next(indices), must_load_from_database, alpha_old, alpha_future)
            self.UpdateFluidVariable(future_step_dataset_name + '/dvyz', VELOCITY_Y_GRADIENT_Z, next(indices), must_load_from_database, alpha_old, alpha_future)
            self.UpdateFluidVariable(future_step_dataset_name + '/dvzx', VELOCITY_Z_GRADIENT_X, next(indices), must_load_from_database, alpha_old, alpha_future)
            self.UpdateFluidVariable(future_step_dataset_name + '/dvzy', VELOCITY_Z_GRADIENT_Y, next(indices), must_load_from_database, alpha_old, alpha_future)
            self.UpdateFluidVariable(future_step_dataset_name + '/dvzz', VELOCITY_Z_GRADIENT_Z, next(indices), must_load_from_database, alpha_old, alpha_future)
        Say('Finished loading fluid from hdf5 file.\n')

class ParticleHistoryLoader:
    def __init__(self, particles_model_part, particle_watcher, pp, main_path):
        self.pp = pp
        self.model_part = particles_model_part
        self.particle_watcher = particle_watcher
        self.main_path = main_path
        self.particles_list_file_name = self.main_path + '/all_particles.hdf5'
        self.prerun_fluid_file_name = pp.CFD_DEM.AddEmptyValue("prerun_fluid_file_name").GetString()
        self.CreateAllParticlesFileIfNecessary()

    def CreateAllParticlesFileIfNecessary(self):
        if not self.pp.CFD_DEM["full_particle_history_watcher"].GetString() == 'Empty':
            nodes = [node for node in self.model_part.Nodes if node.IsNot(BLOCKED)]
            Ids = np.array([node.Id for node in nodes], dtype = int)
            X0s = np.array([node.X0 for node in nodes])
            Y0s = np.array([node.Y0 for node in nodes])
            Z0s = np.array([node.Z0 for node in nodes])
            radii = np.array([node.GetSolutionStepValue(RADIUS) for node in nodes])
            times = np.array([0.0 for node in nodes])

            with h5py.File(self.particles_list_file_name) as f:
                WriteDataToFile(file_or_group = f,
                                names = ['Id', 'X0', 'Y0', 'Z0', 'RADIUS', 'TIME'],
                                data = [Ids, X0s, Y0s, Z0s, radii, times])

    def UpdateListOfAllParticles(self):
        Ids, X0s, Y0s, Z0s, radii, times = [], [], [], [], [], []
        self.particle_watcher.GetNewParticlesData(Ids, X0s, Y0s, Z0s, radii, times)
        names = ['Id', 'X0', 'Y0', 'Z0', 'RADIUS', 'TIME']
        data = [Ids, X0s, Y0s, Z0s, radii, times]
        new_data = []

        with h5py.File(self.particles_list_file_name) as f:
            for name, datum in zip(names, data):
                old_datum = f['/' + name][:]
                new_datum = np.concatenate((old_datum, datum))
                new_data.append(new_datum)

            WriteDataToFile(file_or_group = f, names = names, data = new_data)

    # This function records the particles initial positions in an hdf5 file.
    # It records both the particles that fall within an input bounding box
    # and all the particles in the model part.
    def RecordParticlesInBox(self, bounding_box = BoundingBoxRule()):
        self.bounding_box = bounding_box
        time = self.model_part.ProcessInfo[TIME]

        def IsInside(node):
            is_a_particle = node.IsNot(BLOCKED)
            is_inside = self.bounding_box.CheckIfRuleIsMet(time, node.X, node.Y, node.Z)
            return is_a_particle and is_inside

        nodes_inside = [node for node in self.model_part.Nodes if IsInside(node)]
        Ids_inside = np.array([node.Id for node in nodes_inside])
        X0s_inside = np.array([node.X0 for node in nodes_inside])
        Y0s_inside = np.array([node.Y0 for node in nodes_inside])
        Z0s_inside = np.array([node.Z0 for node in nodes_inside])
        radii_inside = np.array([node.GetSolutionStepValue(RADIUS) for node in nodes_inside])

        if len(radii_inside):
            mean_radius = sum(radii_inside) / len(radii_inside)
        else:
            mean_radius = 1.0

        with h5py.File('particles_snapshot.hdf5') as f:
            prerun_fluid_file_name = self.prerun_fluid_file_name.split('/')[- 1]
            current_fluid = CreateGroup(f, prerun_fluid_file_name)
            # storing the input parameters for this run, the one corresponding
            # to the current pre-calculated fluid
            for k, v in ((k, v) for k, v in json.loads(self.pp.CFD_DEM.WriteJsonString()).items() if 'comment' not in k):
                current_fluid.attrs[k] = v

            snapshot_name = 't=' + str(round(time, 3)) + '_RADIUS=' + str(round(mean_radius, 4)) + '_in_box'
            snapshot = CreateGroup(current_fluid, snapshot_name)
            snapshot.attrs['time'] = time

            names = ['Id', 'X0', 'Y0', 'Z0', 'RADIUS']
            data = [Ids_inside, X0s_inside, Y0s_inside, Z0s_inside, radii_inside]

            for dset_name, datum in zip(names, data):
                CreateDataset(snapshot, dset_name, datum)
