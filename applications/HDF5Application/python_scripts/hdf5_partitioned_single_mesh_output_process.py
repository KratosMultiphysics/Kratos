import KratosMultiphysics
import KratosMultiphysics.HDF5Application as KratosHDF5
import hdf5_output

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return HDF5PartitionedSingleMeshOutputProcess(Model, settings["Parameters"])

class HDF5PartitionedSingleMeshOutputProcess(KratosMultiphysics.Process):
    """A process for writing partitioned simulation results for a single mesh to HDF5."""

    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)
        default_settings = KratosMultiphysics.Parameters("""
            {
                "model_part_name" : "please_specify_model_part_name",
                "list_of_variables" : [],
                "output_time_frequency": 1.0,
                "output_step_frequency": 0,
                "echo_level" : 0
            }
            """)
        settings.ValidateAndAssignDefaults(default_settings)
        self.model_part = Model[settings["model_part_name"].GetString()]
        self._list_of_variables = [settings["list_of_variables"][i].GetString() for i in range(settings["list_of_variables"].size())]
        output_time_frequency = settings["output_time_frequency"].GetDouble()
        output_step_frequency = settings["output_step_frequency"].GetInt()
        self._echo_level = settings["echo_level"].GetInt()
        self._hdf5_file_factory = hdf5_output.HDF5ParallelFileFactory("truncate", self._echo_level)
        self._temporal_output = hdf5_output.TemporalOutputProcess(self.model_part, self._hdf5_file_factory, output_time_frequency, output_step_frequency)
        self._temporal_output.AddOutputObject(hdf5_output.PartitionedNodalResultsOutput(self._list_of_variables))

    def ExecuteBeforeSolutionLoop(self):
        file_name = self.model_part.Name + ".h5"
        hdf5_file = self._hdf5_file_factory.Open(file_name)
        hdf5_output.PartitionedModelPartOutput().Execute(self.model_part, hdf5_file)
        hdf5_output.PartitionedNodalResultsOutput(self._list_of_variables).Execute(self.model_part, hdf5_file)
        self._temporal_output.ExecuteBeforeSolutionLoop()

    def ExecuteFinalizeSolutionStep(self):
        self._temporal_output.ExecuteFinalizeSolutionStep()
