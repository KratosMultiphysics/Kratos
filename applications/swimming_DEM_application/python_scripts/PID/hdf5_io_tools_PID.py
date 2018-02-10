import hdf5_io_tools
import average_field
import math

BaseLoader = hdf5_io_tools.FluidHDF5Loader

class FluidHDF5LoaderPID(BaseLoader):
    def __init__(self, fluid_model_part, particles_model_part, pp, main_path):
        BaseLoader.__init__(self, fluid_model_part, particles_model_part, pp, main_path)
        self.original_file_name = 'mesh_142752_nodes.hdf5'
        self.dataset_name = 'stationary_field'
        self.averager = average_field.Averager(rotation_axis_initial_point = [0., 0., 0.],
                                               rotation_axis_final_point = [0., 0., 1.],
                                               angular_velocity_module = - 2 * math.pi,
                                               dataset_name = self.dataset_name,
                                               original_file_name = self.original_file_name,
                                               original_file_path = main_path,
                                               initial_time = 0.04,
                                               steps_per_average_step = 1,
                                               calculate_standard_deviations = True,
                                               normalize_standard_deviation = True,
                                               overwrite_previous=True)

    def GetDatasetName(self, time_index_future):
        return self.dataset_name

    def LoadFluid(self, fluid_time):
        self.averager.PerformAverage(reference_time = fluid_time)
        self.file_path = self.averager.GetFilePath()
        BaseLoader.LoadFluid(self, fluid_time)

