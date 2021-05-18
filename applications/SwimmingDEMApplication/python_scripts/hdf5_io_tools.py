import bisect as bi
import numpy as np
import h5py
import KratosMultiphysics as Kratos
import KratosMultiphysics.SwimmingDEMApplication as SDEM
from KratosMultiphysics.DEMApplication.DEM_procedures import KratosPrintInfo as Say
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

def CreateGroup(file_or_group, name, overwrite_previous = True):
    if name in file_or_group:
        if overwrite_previous:
            file_or_group['/'].__delitem__(name)
        else:
            return file_or_group['/' + name]

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

    def __init__(self, project_parameters, fluid_model_part, particles_model_part, main_path):
        self.n_nodes = len(fluid_model_part.Nodes)
        self.shape = (self.n_nodes,)
        self.store_pressure = project_parameters["store_fluid_pressure_option"].GetBool()
        self.store_gradient = project_parameters["store_full_gradient_option"].GetBool()
        self.load_derivatives = project_parameters["load_derivatives"].GetBool()
        self.there_are_more_steps_to_load = True
        self.main_path = main_path
        self.fluid_model_part = fluid_model_part
        self.disperse_phase_model_part = particles_model_part

        number_of_variables = 3

        if project_parameters["store_fluid_pressure_option"].GetBool():
            number_of_variables += 1
        if self.load_derivatives or self.store_gradient:
            number_of_variables += 9

        self.extended_shape = self.shape + (number_of_variables, )
        self.file_name = self.GetFileName()
        self.file_path = main_path + self.file_name

        if project_parameters["custom_fluid"]["fluid_already_calculated"].GetBool():

            with h5py.File(self.file_path, 'r') as f:
                nodes_ids = np.array([node_id for node_id in f['nodes'][:, 0]])
                self.permutations = np.array(range(len(nodes_ids)))
                # obtaining the vector of permutations by ordering [0, 1, ..., n_nodes] as nodes_ids, by increasing order of id.
                self.permutations = np.array([x for (y, x) in sorted(zip(nodes_ids, self.permutations))])
                self.CheckTimes(f)

            self.data_array_past = np.zeros(self.extended_shape)
            self.data_array_future = np.zeros(self.extended_shape)
            self.time_index_past = - 1 # it starts at an absurd value
            self.time_index_future = - 1
            viscosity = 1e-6
            density = 1000. # BIG TODO: READ THIS FROM NODES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            for node in self.fluid_model_part.Nodes:
                node.SetSolutionStepValue(Kratos.VISCOSITY, viscosity)
                node.SetSolutionStepValue(Kratos.DENSITY, density)
        else:
            self.dtype = np.float64
            if project_parameters["store_fluid_in_single_precision"].GetBool():
                self.dtype = np.float32

            self.compression_type = 'gzip'
            for node in self.fluid_model_part.Nodes:
                viscosity = node.GetSolutionStepValue(Kratos.VISCOSITY)
                density = node.GetSolutionStepValue(Kratos.DENSITY)
                break

            with h5py.File(self.file_path, 'w') as f:
                f.attrs['kinematic viscosity'] = viscosity
                f.attrs['time step'] = project_parameters.Dt
                f.attrs['density'] = density
                solver_settings = project_parameters["fluid_parameters"]["solver_settings"]
                f.attrs['solver type'] = solver_settings["solver_type"].GetString()
                f.attrs['linear system solver type'] = solver_settings["linear_solver_settings"]["solver_type"].GetString()
                f.attrs['use orthogonal subscales'] = False
                nodes = np.array([(node.Id, node.X, node.Y, node.Z) for node in fluid_model_part.Nodes])
                f.create_dataset(name = 'nodes', compression = self.compression_type, data = nodes, dtype = np.float64)

        self.last_time = 0.0

        self.current_data_array = np.zeros(self.extended_shape)

    def GetFileName(self):
        return self.project_parameters.AddEmptyValue("prerun_fluid_file_name").GetString()

    def CheckTimes(self, hdf5_file):
        self.times_str = list([str(key) for key in hdf5_file.keys() if 'time' in hdf5_file['/' + key].attrs])
        self.times = np.array([float(hdf5_file[key].attrs['time']) for key in self.times_str])

        if len(self.times) < 2:
            raise ValueError("\nThere are only " + str(len(self.times)) + ' time steps stored in the hdf5 file. At least two are needed.\n')

        self.times_str = np.array([x for (y, x) in sorted(zip(self.times, self.times_str))])
        self.times = sorted(self.times)
        self.dt = self.times[-1] - self.times[-2]

    def GetTimeIndicesAndWeights(self, current_time):
        index_future = bi.bisect(self.times, current_time)
        index_past = max(0, index_future - 1)
        time_past = self.times[index_past]

        if index_future == len(self.times): # we are beyond the last time
            alpha_past = 0
            alpha_future = 1
            index_future = index_past
            self.there_are_more_steps_to_load = False
        else:
            alpha_future = max(0, (current_time - time_past) / self.dt)
            alpha_past = 1.0 - alpha_future

        return index_past, alpha_past, index_future, alpha_future

    def CanLoadMoreSteps(self):
        return self.there_are_more_steps_to_load

    def FillUpSingleDataset(self, name, variable, variable_index_in_temp_array):
        with h5py.File(self.file_path, 'r+') as f:
            f.create_dataset(name, compression = self.compression_type, shape = self.shape, dtype = self.dtype)
            for i_node, node in enumerate(self.fluid_model_part.Nodes):
                self.current_data_array[i_node, variable_index_in_temp_array] = node.GetSolutionStepValue(Kratos.variable)
            f[name][:] = self.current_data_array[:, variable_index_in_temp_array]

    def FillFluidDataStep(self):
        time = self.fluid_model_part.ProcessInfo[Kratos.TIME]
        name = str(time)
        with h5py.File(self.file_path) as f:
            f.create_group(name = name)
            f[name].attrs['time'] = time

        index = Index()
        if not self.last_time == time:
            self.FillUpSingleDataset(name + '/vx', Kratos.VELOCITY_X, next(index))
            self.FillUpSingleDataset(name + '/vy', Kratos.VELOCITY_Y, next(index))
            self.FillUpSingleDataset(name + '/vz', Kratos.VELOCITY_Z, next(index))

            if self.store_pressure:
                self.FillUpSingleDataset(name + '/p', Kratos.PRESSURE, next(index))

        if self.store_gradient:
            self.FillUpSingleDataset(name + '/dvxx', Kratos.VELOCITY_X_GRADIENT_X, next(index))
            self.FillUpSingleDataset(name + '/dvxy', Kratos.VELOCITY_X_GRADIENT_Y, next(index))
            self.FillUpSingleDataset(name + '/dvxz', Kratos.VELOCITY_X_GRADIENT_Z, next(index))
            self.FillUpSingleDataset(name + '/dvyx', Kratos.VELOCITY_Y_GRADIENT_X, next(index))
            self.FillUpSingleDataset(name + '/dvyy', Kratos.VELOCITY_Y_GRADIENT_Y, next(index))
            self.FillUpSingleDataset(name + '/dvyz', Kratos.VELOCITY_Y_GRADIENT_Z, next(index))
            self.FillUpSingleDataset(name + '/dvzx', Kratos.VELOCITY_Z_GRADIENT_X, next(index))
            self.FillUpSingleDataset(name + '/dvzy', Kratos.VELOCITY_Z_GRADIENT_Y, next(index))
            self.FillUpSingleDataset(name + '/dvzz', Kratos.VELOCITY_Z_GRADIENT_Z, next(index))

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

    def UpdateFluidVariable(self,
                            name,
                            variable,
                            variable_index,
                            must_load_future_values_from_database,
                            alpha_past,
                            alpha_future):
        if must_load_future_values_from_database:
            with h5py.File(self.file_path, 'r') as f:
                self.data_array_future[:, variable_index] = self.ConvertComponent(f, name)

        self.current_data_array[:, variable_index] = (
            alpha_past * self.data_array_past[:, variable_index]
          + alpha_future * self.data_array_future[:, variable_index])

        for i, node in enumerate(self.fluid_model_part.Nodes):
            node.SetSolutionStepValue(Kratos.variable, self.current_data_array[i, variable_index])

    def GetDatasetName(self, time_index_future):
        return self.times_str[time_index_future]

    def LoadFluid(self, fluid_time):
        Say('\nLoading fluid from hdf5 file...')
        # getting time indices and weights (identifying the two fluid time steps surrounding the current DEM step and assigning correspnding weights)
        time_index_past, alpha_past, time_index_future, alpha_future = self.GetTimeIndicesAndWeights(fluid_time)
        future_step_dataset_name = self.GetDatasetName(time_index_future)
        must_load_from_database = self.time_index_past != time_index_past or self.time_index_future != time_index_future# old and future time steps must be updated

        if must_load_from_database: # the current time is not between the two already loaded time steps
            # the old future becomes the new past
            self.data_array_past, self.data_array_future = self.data_array_future, self.data_array_past

        index = Index()
        self.UpdateFluidVariable(future_step_dataset_name + '/vx', Kratos.VELOCITY_X, next(index), must_load_from_database, alpha_past, alpha_future)
        self.UpdateFluidVariable(future_step_dataset_name + '/vy', Kratos.VELOCITY_Y, next(index), must_load_from_database, alpha_past, alpha_future)
        self.UpdateFluidVariable(future_step_dataset_name + '/vz', Kratos.VELOCITY_Z, next(index), must_load_from_database, alpha_past, alpha_future)

        if self.store_pressure:
            self.UpdateFluidVariable(future_step_dataset_name + '/p', Kratos.PRESSURE, next(index), must_load_from_database, alpha_past, alpha_future)

        if self.load_derivatives:
            self.UpdateFluidVariable(future_step_dataset_name + '/dvxx', Kratos.VELOCITY_X_GRADIENT_X, next(index), must_load_from_database, alpha_past, alpha_future)
            self.UpdateFluidVariable(future_step_dataset_name + '/dvxy', Kratos.VELOCITY_X_GRADIENT_Y, next(index), must_load_from_database, alpha_past, alpha_future)
            self.UpdateFluidVariable(future_step_dataset_name + '/dvxz', Kratos.VELOCITY_X_GRADIENT_Z, next(index), must_load_from_database, alpha_past, alpha_future)
            self.UpdateFluidVariable(future_step_dataset_name + '/dvyx', Kratos.VELOCITY_Y_GRADIENT_X, next(index), must_load_from_database, alpha_past, alpha_future)
            self.UpdateFluidVariable(future_step_dataset_name + '/dvyy', Kratos.VELOCITY_Y_GRADIENT_Y, next(index), must_load_from_database, alpha_past, alpha_future)
            self.UpdateFluidVariable(future_step_dataset_name + '/dvyz', Kratos.VELOCITY_Y_GRADIENT_Z, next(index), must_load_from_database, alpha_past, alpha_future)
            self.UpdateFluidVariable(future_step_dataset_name + '/dvzx', Kratos.VELOCITY_Z_GRADIENT_X, next(index), must_load_from_database, alpha_past, alpha_future)
            self.UpdateFluidVariable(future_step_dataset_name + '/dvzy', Kratos.VELOCITY_Z_GRADIENT_Y, next(index), must_load_from_database, alpha_past, alpha_future)
            self.UpdateFluidVariable(future_step_dataset_name + '/dvzz', Kratos.VELOCITY_Z_GRADIENT_Z, next(index), must_load_from_database, alpha_past, alpha_future)

        if self.time_index_past == - 1: # it is the first upload
            self.data_array_past[:] = self.data_array_future[:]

        self.time_index_past = time_index_past
        self.time_index_future = time_index_future

        Say('Finished loading fluid from hdf5 file.\n')

class ParticleHistoryLoader:
    def __init__(self, parameters, particles_model_part, particle_watcher, main_path):
        self.parameters = parameters
        self.model_part = particles_model_part
        self.particle_watcher = particle_watcher
        self.main_path = main_path
        self.particles_list_file_name = self.main_path + '/all_particles.hdf5'
        self.prerun_fluid_file_name = parameters.AddEmptyValue("prerun_fluid_file_name").GetString()

        self.CreateAllParticlesFileIfNecessary()
        self.run_code = None

    def CreateAllParticlesFileIfNecessary(self):
        if not self.parameters["full_particle_history_watcher"].GetString() == 'Empty':
            nodes = [node for node in self.model_part.Nodes if node.IsNot(Kratos.BLOCKED)]
            Ids = np.array([node.Id for node in nodes], dtype = int)
            X0s = np.array([node.X0 for node in nodes])
            Y0s = np.array([node.Y0 for node in nodes])
            Z0s = np.array([node.Z0 for node in nodes])
            radii = np.array([node.GetSolutionStepValue(Kratos.RADIUS) for node in nodes])
            times = np.array([0.0 for node in nodes])

            with h5py.File(self.particles_list_file_name, 'w') as f:
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
    def RecordParticlesInBox(self, bounding_box = SDEM.BoundingBoxRule()):
        self.bounding_box = bounding_box
        time = self.model_part.ProcessInfo[Kratos.TIME]

        def IsInside(node):
            is_a_particle = node.IsNot(Kratos.BLOCKED)
            is_inside = self.bounding_box.CheckIfRuleIsMet(time, node.X, node.Y, node.Z)
            return is_a_particle and is_inside

        nodes_inside = [node for node in self.model_part.Nodes if IsInside(node)]
        Ids_inside = np.array([node.Id for node in nodes_inside])
        X0s_inside = np.array([node.X0 for node in nodes_inside])
        Y0s_inside = np.array([node.Y0 for node in nodes_inside])
        Z0s_inside = np.array([node.Z0 for node in nodes_inside])
        radii_inside = np.array([node.GetSolutionStepValue(Kratos.RADIUS) for node in nodes_inside])

        if len(radii_inside):
            mean_radius = sum(radii_inside) / len(radii_inside)
        else:
            mean_radius = 1.0

        with h5py.File(self.main_path + '/particles_snapshots.hdf5') as f:
            prerun_fluid_file_name = self.prerun_fluid_file_name.split('/')[- 1]
            current_fluid = CreateGroup(f, prerun_fluid_file_name, overwrite_previous = False)

            # snapshot_name = 't=' + str(round(time, 3)) + '_RADIUS=' + str(round(mean_radius, 4)) + '_in_box'
            snapshot_name = str(len(current_fluid.items()) + 1)
            self.run_code = prerun_fluid_file_name.strip('.hdf5') + '_' + snapshot_name

            snapshot = CreateGroup(current_fluid, snapshot_name)
            snapshot.attrs['time'] = time
            snapshot.attrs['particles_nondimensional_radius'] = mean_radius
            # storing the input parameters for this run, the one corresponding
            # to the current pre-calculated fluid
            for k, v in ((k, v) for k, v in json.loads(self.parameters.WriteJsonString()).items() if 'comment' not in k):
                snapshot.attrs[k] = v

            names = ['Id', 'X0', 'Y0', 'Z0', 'RADIUS']
            data = [Ids_inside, X0s_inside, Y0s_inside, Z0s_inside, radii_inside]

            for dset_name, datum in zip(names, data):
                CreateDataset(snapshot, dset_name, datum)

    def GetRunCode(self):
        if self.run_code == None:
            raise RuntimeError('No run code has been generated so far, because no snapshot has been performed yet')
        else:
            return self.run_code
