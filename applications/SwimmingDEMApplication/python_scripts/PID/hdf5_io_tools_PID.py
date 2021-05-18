import KratosMultiphysics.SwimmingDEMApplication.hdf5_io_tools as hdf5_io_tools
import math

BaseLoader = hdf5_io_tools.FluidHDF5Loader

class FluidHDF5LoaderPID(BaseLoader):
    def __init__(self,
                 parameters,
                 fluid_model_part,
                 particles_model_part,
                 main_path,
                 averager):
        BaseLoader.__init__(self, fluid_model_part, particles_model_part, parameters, main_path)
        self.dataset_name = 'stationary_field'
        self.averager = averager

    def GetFileName(self):
        if not self.project_parameters["averaging_has_already_been_done"].GetBool():
            return BaseLoader.GetFileName(self)
        else:
            original_file_name = self.project_parameters.AddEmptyValue("prerun_fluid_file_name").GetString()
            return original_file_name.replace('.hdf5', '') + '_averaged.hdf5'

    def CheckTimes(self, hdf5_file):
        if not self.project_parameters["averaging_has_already_been_done"].GetBool():
            BaseLoader.CheckTimes(self, hdf5_file)
        else:
            pass

    def GetTimeIndicesAndWeights(self, current_time):
        return 0, 0.0, 1, 1.0

    def GetDatasetName(self, time_index_future):
        return self.dataset_name

    def LoadFluid(self, fluid_time):
        self.file_path = self.averager.GetFilePath()
        BaseLoader.LoadFluid(self, fluid_time)

