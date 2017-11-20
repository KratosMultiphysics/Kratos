import KratosMultiphysics
import KratosMultiphysics.HDF5Application as KratosHDF5
import os

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return Hdf5FixedMeshOutputProcess(Model, settings["Parameters"])

class Hdf5FixedMeshOutputProcess(KratosMultiphysics.Process):
    """A process for writing simulation results on a fixed mesh to HDF5."""

    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)
        default_settings = KratosMultiphysics.Parameters("""
            {
                "model_part_name" : "please_specify_model_part_name",
                "list_of_elements" : [],
                "list_of_conditions" : [],
                "list_of_variables" : [],
                "output_time_frequency": 1.0,
                "output_step_frequency": 0,
                "partitioned" : false,
                "time_tag_precision" : 4,
                "echo_level" : 0
            }
            """)
        settings.ValidateAndAssignDefaults(default_settings)
        self.model_part = Model[settings["model_part_name"].GetString()]
        self._list_of_elements = settings["list_of_elements"].Clone()
        self._list_of_conditions = settings["list_of_conditions"].Clone()
        self._list_of_variables = settings["list_of_variables"].Clone()
        self._output_time_frequency = settings["output_time_frequency"].GetDouble()
        self._output_step_frequency = settings["output_step_frequency"].GetInt()
        self._output_time = 0.0
        self._output_step = 0
        self._partitioned = settings["partitioned"].GetBool()
        self._echo_level = settings["echo_level"].GetInt()
        self._time_tag_precision = settings["time_tag_precision"].GetInt()

    def _get_hdf5_file(self, file_name):
        params = KratosMultiphysics.Parameters("""{
            "file_name" : "",
            "file_access_mode" : "",
            "file_driver" : "",
            "echo_level" : 0
        }""")
        params["file_name"].SetString(file_name)
        params["file_access_mode"].SetString("truncate")
        params["echo_level"].SetInt(self._echo_level)
        if self._partitioned:
            params["file_driver"].SetString("mpio")
            return KratosHDF5.HDF5FileParallel(params)
        else:
            params["file_driver"].SetString("sec2")
            return KratosHDF5.HDF5FileSerial(params)

    def _get_model_part_io(self, hdf5_file):
        params = KratosMultiphysics.Parameters("""{
            "prefix" : "/ModelData"
        }""")
        params.AddValue("list_of_elements", self._list_of_elements)
        params.AddValue("list_of_conditions", self._list_of_conditions)
        if self._partitioned:
            return KratosHDF5.HDF5PartitionedModelPartIO(params, hdf5_file)
        else:
            return KratosHDF5.HDF5ModelPartIO(params, hdf5_file)

    def _get_nodal_results_io(self, hdf5_file):
        params = KratosMultiphysics.Parameters("""{
            "prefix" : "/ResultsData",
            "partitioned": ""
        }""")
        params.AddValue("list_of_variables", self._list_of_variables)
        if self._partitioned:
            params["partitioned"].SetBool(True)
        else:
            params["partitioned"].SetBool(False)
        return KratosHDF5.HDF5NodalSolutionStepDataIO(params, hdf5_file)

    def ExecuteBeforeSolutionLoop(self):
        hdf5_file = self._get_hdf5_file(self.model_part.Name + ".h5")
        self._get_model_part_io(hdf5_file).WriteModelPart(self.model_part)
        self._get_nodal_results_io(hdf5_file).WriteNodalResults(self.model_part.Nodes, 0)
        self._output_time = 0.0
        self._output_step = 0

    def ExecuteFinalizeSolutionStep(self):
        delta_time = self.model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
        self._output_time += delta_time
        self._output_step += 1
        if self._output_time >= self._output_time_frequency or self._output_step == self._output_step_frequency:
            fmt = "{:." + str(self._time_tag_precision) + "f}"
            time_tag = "-" + fmt.format(self.model_part.ProcessInfo[KratosMultiphysics.TIME])
            hdf5_file = self._get_hdf5_file(self.model_part.Name + time_tag + ".h5")
            self._get_nodal_results_io(hdf5_file).WriteNodalResults(self.model_part.Nodes, 0)
            self._output_time = 0.0
            self._output_step = 0
