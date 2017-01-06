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
    def __init__(self, fluid_model_part, pp):
        self.n_nodes = len(fluid_model_part.Nodes)
        self.shape = (self.n_nodes,)
        self.store_pressure = pp.CFD_DEM.store_fluid_pressure
        
        if pp.CFD_DEM.store_fluid_pressure:
            self.extended_shape = self.shape + (4,)
        else:
            self.extended_shape = self.shape + (3,)
            
        self.file_name = '../fluid_results.hdf5'
        self.fluid_model_part = fluid_model_part
        
        if pp.CFD_DEM.fluid_already_calculated:
            with h5py.File(self.file_name, 'r') as f:
                self.times_str = np.array([key for key in list(f.keys())])
            self.times     = np.array([float(key) for key in self.times_str])
            self.old_data_array = np.zeros(self.extended_shape)
            self.future_data_array = np.zeros(self.extended_shape)
            self.old_time_index = 0
            self.future_time_index = 1
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

    def FillUpSingleDataset(self, name, variable, variable_index_in_temp_array):
        with h5py.File(self.file_name, 'r+') as f:
            f.create_dataset(name, compression = self.compression_type, shape = self.shape, dtype = self.dtype)
            i_node = 0
            for node in self.fluid_model_part.Nodes:
                self.current_data_array[i_node, variable_index_in_temp_array] = node.GetSolutionStepValue(variable)
                i_node += 1
            f[name][:] = self.current_data_array[:, variable_index_in_temp_array]
        
    def FillFluidDataStep(self):
        time = self.fluid_model_part.ProcessInfo[TIME]
        name = str(time)

        if not self.last_time == time:
            self.FillUpSingleDataset(name + '/vx', VELOCITY_X, 0)
            self.FillUpSingleDataset(name + '/vy', VELOCITY_Y, 1)
            self.FillUpSingleDataset(name + '/vz', VELOCITY_Z, 2)
            
            if self.store_pressure:
                self.FillUpSingleDataset(name + '/p', PRESSURE, 3)
                
            self.last_time = time

    def UpdateFluidVariable(self, name, variable, variable_index_in_temp_array, must_load_future_values_from_database, alpha_old, alpha_future):        
        if must_load_future_values_from_database:
            with h5py.File(self.file_name, 'r') as f:
                self.future_data_array[:, variable_index_in_temp_array] = f[name]
        
        self.current_data_array[:, variable_index_in_temp_array] = alpha_old * self.old_data_array[:, variable_index_in_temp_array] + alpha_future * self.future_data_array[:, variable_index_in_temp_array]
        
        i_node = 0

        for node in self.fluid_model_part.Nodes:
            node.SetSolutionStepValue(variable, self.current_data_array[i_node, variable_index_in_temp_array])
            i_node += 1           
           
    def LoadFluid(self, DEM_time):   
        # getting time indices and weights (identifyint the two fluid time steps surrounding the current DEM step and assigning correspnding weights) 
        old_time_index, alpha_old, future_time_index, alpha_future = GetOldTimeIndicesAndWeights(DEM_time, self.times)
        old_step_dataset_name    = self.times_str[old_time_index]
        future_step_dataset_name = self.times_str[future_time_index]
        must_load_from_database = not self.old_time_index == old_time_index # old and future time steps must be updated
        
        if must_load_from_database: 
            # new old becomes old future
            self.old_data_array, self.future_data_array = self.future_data_array, self.old_data_array
            
        self.UpdateFluidVariable(future_step_dataset_name + '/vx', VELOCITY_X, 0, must_load_from_database, alpha_old, alpha_future)
        self.UpdateFluidVariable(future_step_dataset_name + '/vy', VELOCITY_Y, 1, must_load_from_database, alpha_old, alpha_future)  
        self.UpdateFluidVariable(future_step_dataset_name + '/vz', VELOCITY_Z, 2, must_load_from_database, alpha_old, alpha_future)  

        if self.store_pressure:
            self.UpdateFluidVariable(future_step_dataset_name + '/p', PRESSURE, 3, must_load_from_database, alpha_old, alpha_future)