import bisect as bi
import numpy as np
import h5py
from KratosMultiphysics import *


def GetOldTimeIndicesAndWeights(current_time, times_array):
    old_index = bi.bisect(times_array, current_time)
    future_index = old_index + 1
    old_time =  times_array[old_index]

    if future_index >= len(times_array):
        alpha_old = 1
        future_index = old_index
        alpha_future = 0
    else:
        fluid_dt = times_array[future_index] - old_time
        alpha_old = (current_time - old_time) / fluid_dt
        alpha_future = 1.0 - alpha_old

    return old_index, alpha_old, future_index, alpha_future

class FluidHDF5Loader:
    def __init__(self, fluid_model_part, pp, main_path):
        self.n_nodes = len(fluid_model_part.Nodes)
        self.shape = (self.n_nodes,)
        self.store_pressure = pp.CFD_DEM.store_fluid_pressure
        self.store_gradient = pp.CFD_DEM.store_full_gradient

        number_of_variables = 3

        if pp.CFD_DEM.store_fluid_pressure:
            number_of_variables += 1
        if pp.CFD_DEM.load_derivatives:
            number_of_variables += 9

        self.extended_shape = self.shape + (number_of_variables,)

        self.file_name = main_path + '/fluid_results.hdf5'
        self.fluid_model_part = fluid_model_part

        if pp.CFD_DEM.fluid_already_calculated:

            with h5py.File(self.file_name, 'r') as f:
                self.times_str = np.array([key for key in f.keys() if key not in {'density', 'viscosity', 'nodes', 'number of nodes'}])
                nodes_ids = np.array([node_id for node_id in f['nodes'][:, 0]])
                self.permutations = np.array(range(len(nodes_ids)))
                # obtaining the vector of permutations by ordering [0, 1, ..., n_nodes] as nodes_ids, by increasing order of id.
                self.permutations = np.array([x for (y, x) in sorted(zip(nodes_ids, self.permutations))])

            self.times     = np.array([float(key) for key in self.times_str])
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
            with h5py.File(self.file_name, 'w') as f:
                f.attrs['kinematic viscosity'] = viscosity
                f.attrs['time step'] = pp.Dt
                f.attrs['density'] = density
                f.attrs['solver type'] = pp.FluidSolverConfiguration.solver_type
                f.attrs['linear system solver type'] = pp.FluidSolverConfiguration.linear_solver_config.solver_type
                f.attrs['use orthogonal subscales'] = bool(pp.FluidSolverConfiguration.oss_switch)
                self.dtype = np.float64

            if pp.CFD_DEM.store_fluid_in_single_precision:
                self.dtype = np.float32

        self.last_time = 0.0

        self.current_data_array = np.zeros(self.extended_shape)

    def Index(self):
        index = 0
        while True:
            yield index
            index += 1

    def FillUpSingleDataset(self, name, variable, variable_index_in_temp_array):
        with h5py.File(self.file_name, 'r+') as f:
            f.create_dataset(name, compression = self.compression_type, shape = self.shape, dtype = self.dtype)
            for i_node, node in enumerate(self.fluid_model_part.Nodes):
                self.current_data_array[i_node, variable_index_in_temp_array] = node.GetSolutionStepValue(variable)
            f[name][:] = self.current_data_array[:, variable_index_in_temp_array]

    def FillFluidDataStep(self):
        time = self.fluid_model_part.ProcessInfo[TIME]
        name = str(time)

        indices = self.Index()
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
            return np.array(f[component_name.replace('/vx', '/VELOCITY')][:, 0])
        elif '/vy' in component_name:
            return np.array(f[component_name.replace('/vy', '/VELOCITY')][:, 1])
        elif '/vz' in component_name:
            return np.array(f[component_name.replace('/vz', '/VELOCITY')][:, 2])
        else:
            return np.transpose(f[component_name.replace('/p', '/PRESSURE')])

    def UpdateFluidVariable(self, name, variable, variable_index_in_temp_array, must_load_future_values_from_database, alpha_old, alpha_future):
        if must_load_future_values_from_database:
            with h5py.File(self.file_name, 'r') as f:
                self.future_data_array[:, variable_index_in_temp_array] = self.ConvertComponent(f, name)[:]
                for i, j in enumerate(self.permutations):
                    self.future_data_array[i, variable_index_in_temp_array] = self.future_data_array[j, variable_index_in_temp_array]

        self.current_data_array[:, variable_index_in_temp_array] = alpha_old * self.old_data_array[:, variable_index_in_temp_array] + alpha_future * self.future_data_array[:, variable_index_in_temp_array]

        for i_node, node in enumerate(self.fluid_model_part.Nodes):
            node.SetSolutionStepValue(variable, self.current_data_array[i_node, variable_index_in_temp_array])

    def LoadFluid(self, DEM_time):
        # getting time indices and weights (identifyint the two fluid time steps surrounding the current DEM step and assigning correspnding weights)
        old_time_index, alpha_old, future_time_index, alpha_future = GetOldTimeIndicesAndWeights(DEM_time, self.times)
        old_step_dataset_name    = self.times_str[old_time_index]
        future_step_dataset_name = self.times_str[future_time_index]
        must_load_from_database = not self.old_time_index == old_time_index # old and future time steps must be updated

        if must_load_from_database:
            # new old becomes old future
            self.old_data_array, self.future_data_array = self.future_data_array, self.old_data_array

        indices = self.Index()
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
