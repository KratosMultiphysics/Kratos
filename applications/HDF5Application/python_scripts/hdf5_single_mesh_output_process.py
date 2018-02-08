import KratosMultiphysics
import KratosMultiphysics.HDF5Application as KratosHDF5
import hdf5_output

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return HDF5SingleMeshOutputProcess(Model, settings["Parameters"])

class HDF5SingleMeshOutputProcess(KratosMultiphysics.Process):
    """A process for writing simulation results for a single mesh to HDF5."""

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
        self._temporal_output = hdf5_output.TemporalOutputProcess(self.model_part, output_time_frequency, output_step_frequency)
        self._temporal_output.list_of_writers.append(hdf5_output.NodalResultsWriter(self._list_of_variables, echo_level))

    def ExecuteBeforeSolutionLoop(self):
        file_name = self.model_part.Name + ".h5"
        hdf5_output.ModelPartWriter(self._echo_level).Write(file_name, self.model_part)
        hdf5_output.NodalResultsWriter(self._list_of_variables, self._echo_level).Write(file_name, self.model_part)
        self._temporal_output.ExecuteBeforeSolutionLoop()

    def ExecuteFinalizeSolutionStep(self):
        self._temporal_output.ExecuteBeforeSolutionLoop()
