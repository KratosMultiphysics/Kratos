# Core imports
import KratosMultiphysics

# HDF5 imports
import KratosMultiphysics.HDF5Application as HDF5Application
from KratosMultiphysics.HDF5Application.core.operations.aggregated_operations import ControlledOperation
from KratosMultiphysics.HDF5Application.core.operations.aggregated_operations import AggregatedControlledOperations
from KratosMultiphysics.HDF5Application.core.file_io import OpenHDF5File
from KratosMultiphysics.HDF5Application.core.pattern import EvaluatePattern
from KratosMultiphysics.HDF5Application.core.processes import HDF5OutputProcess
from KratosMultiphysics.HDF5Application.core.operations.model_part import VertexValueOutput
from KratosMultiphysics.HDF5Application.core.operations.model_part import VertexCoordinateOutput
from KratosMultiphysics.HDF5Application.core.controllers import SingleTimeController


def Factory(parameters: KratosMultiphysics.Parameters,
            model: KratosMultiphysics.Model):
    if not isinstance(parameters, KratosMultiphysics.Parameters):
        raise Exception("expecting input parameters of type KratosMultiphysics.Parameters, but got {}".format(type(parameters)))
    return PointSetOutputProcess(model, parameters["Parameters"])


class PointSetOutputProcess(KratosMultiphysics.OutputProcess):

    def __init__(self, model: KratosMultiphysics.Model, parameters: KratosMultiphysics.Parameters):
        KratosMultiphysics.OutputProcess.__init__(self)
        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())
        parameters["file_parameters"].ValidateAndAssignDefaults(self.GetDefaultParameters()["file_parameters"])

        self.model_part = model.GetModelPart(parameters["model_part_name"].GetString())
        self.variables = parameters["output_variables"].Clone()

        configuration_name = parameters["search_configuration"].GetString()
        if configuration_name == "initial":
            search_configuration = KratosMultiphysics.Configuration.Initial
        elif configuration_name == "current":
            search_configuration = KratosMultiphysics.Configuration.Current
        else:
            raise RuntimeError("Invalid search configuration: '{}'".format(configuration_name))

        search_tolerance = parameters["search_tolerance"].GetDouble()
        locator = HDF5Application.BruteForcePointLocatorAdaptor(self.model_part,
                                                                search_configuration,
                                                                search_tolerance)

        self.interval_utility = KratosMultiphysics.IntervalUtility(parameters)
        self.output_frequency = parameters["output_frequency"].GetInt()

        self.file_parameters = parameters["file_parameters"].Clone()

        self.is_distributed = self.model_part.IsDistributed()

        self.coordinates_prefix_pattern = parameters["coordinates_prefix"].GetString()
        self.variables_prefix_pattern = parameters["variables_prefix"].GetString()
        self.isHistorical = parameters["historical_value"].GetBool()

        # Create vertices
        self.vertices = HDF5Application.VertexContainer()
        for i_vertex, position in enumerate(parameters["positions"].values()):
            vertex = HDF5Application.Vertex.MakeShared(
                position.GetVector(),
                locator,
                i_vertex)
            if vertex.IsLocated():
                self.vertices.push_back(vertex)


    def IsOutputStep(self):
        time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        step = self.model_part.ProcessInfo[KratosMultiphysics.STEP]
        return self.interval_utility.IsInInterval(time) and (step % self.output_frequency == 0)


    def ExecuteInitialize(self):
        coordinates_path = EvaluatePattern(
            self.coordinates_prefix_pattern,
            self.model_part)

        io_parameters = KratosMultiphysics.Parameters("""{
            "prefix" : "",
            "write_vertex_ids" : true
        }""")
        io_parameters["prefix"].SetString(coordinates_path)
        with OpenHDF5File(self.__GetCurrentFileParameters(), self.model_part) as file:
            # Skip writing if the analysis is restarted and the group already exists.
            if self.model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] and file.HasPath(coordinates_path):
                KratosMultiphysics.Logger.PrintWarning("[PointSetOutputProcess] Path exists", f"Skip writing vertex coordinates to group: {coordinates_path}")
            else:
                HDF5Application.VertexContainerCoordinateIO(io_parameters, file).Write(self.vertices)



    def PrintOutput(self):
        prefix = EvaluatePattern(self.variables_prefix_pattern, self.model_part)

        io_parameters = KratosMultiphysics.Parameters()
        io_parameters.AddString("prefix", prefix)
        io_parameters.AddValue("list_of_variables", self.variables.Clone())

        # Set append mode to avoid overwriting the coordinates
        file_parameters = self.__GetCurrentFileParameters()
        file_parameters["file_access_mode"].SetString("read_write")

        with OpenHDF5File(file_parameters, self.model_part) as file:
            # Skip writing if the analysis is restarted and the group already exists.
            if self.model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] and file.HasPath(prefix):
                KratosMultiphysics.Logger.PrintWarning("[PointSetOutputProcess] Path exists", f"Skip writing vertex data to group: {prefix}")
            else:
                if self.isHistorical:
                    HDF5Application.VertexContainerHistoricalVariableIO(io_parameters, file).Write(self.vertices)
                else:
                    HDF5Application.VertexContainerNonHistoricalVariableIO(io_parameters, file).Write(self.vertices)

    @staticmethod
    def GetDefaultParameters():
        return KratosMultiphysics.Parameters("""{
            "model_part_name"       : "",
            "interval"              : [0.0, "End"],
            "output_frequency"      : 1,
            "positions"             : [[]],
            "output_variables"      : [],
            "historical_value"      : true,
            "search_configuration"  : "initial",
            "search_tolerance"      : 1e-6,
            "coordinates_prefix"    : "/<model_part_name>_point_set_output",
            "variables_prefix"      : "/<model_part_name>_point_set_output/step_<step>",
            "file_parameters"       : {
                "file_name"         : "",
                "file_access_mode"  : "truncate",
                "echo_level"        : 0
            }
        }""")


    def __GetCurrentFileParameters(self):
        parameters = self.file_parameters.Clone()
        file_name = EvaluatePattern(
            parameters["file_name"].GetString(),
            self.model_part)

        parameters["file_name"].SetString(file_name)
        return parameters

class NewPointOutputProcess(HDF5OutputProcess):
    @classmethod
    def GetDefaultParameters(cls) -> KratosMultiphysics.Parameters:
        return KratosMultiphysics.Parameters("""
            {
                "model_part_name": "PLEASE_SPECIFY_MODEL_PART_NAME",
                "file_settings": {
                    "file_name"        : "<model_part_name>-<time>.h5",
                    "time_format"      : "0.4f",
                    "file_access_mode" : "exclusive",
                    "max_files_to_keep": "unlimited",
                    "echo_level"       :  0
                },
                "output_time_settings": {
                    "output_control_type": "time",
                    "output_interval"    : 1.0
                },
                "vertex_output_settings": {
                    "prefix"              : "/VertexData",
                    "time_format"         : "0.4f",
                    "positions"           : [[]],
                    "search_configuration": "initial",
                    "search_tolerance"    : 1e-6,
                },
                "nodal_solution_step_data_settings": {
                    "prefix"           : "/VertexData/VertexSolutionStepData",
                    "list_of_variables": [],
                    "time_format"      : "0.4f"
                },
                "nodal_data_value_settings": {
                    "prefix"           : "/VertexData/VertexDataValues",
                    "list_of_variables": [],
                    "time_format"      : "0.4f"
                }
            }""")

    def __init__(self, model: KratosMultiphysics.Model, parameters: KratosMultiphysics.Parameters) -> None:
        super().__init__()
        parameters.RecursivelyValidateAndAssignDefaults(self.GetDefaultParameters())

        model_part = model[parameters["model_part_name"].GetString()]

        vertex_output_settings = parameters["vertex_output_settings"]

        # create locator
        configuration_name = vertex_output_settings["search_configuration"].GetString()
        if configuration_name == "initial":
            search_configuration = KratosMultiphysics.Configuration.Initial
        elif configuration_name == "current":
            search_configuration = KratosMultiphysics.Configuration.Current
        else:
            raise RuntimeError("Invalid search configuration: '{}'".format(configuration_name))
        search_tolerance = parameters["search_tolerance"].GetDouble()
        locator = HDF5Application.BruteForcePointLocatorAdaptor(model_part,
                                                                search_configuration,
                                                                search_tolerance)

        # create vertices
        self.vertices = HDF5Application.VertexContainer()
        for i_vertex, position in enumerate(vertex_output_settings["positions"].values()):
            vertex = HDF5Application.Vertex.MakeShared(
                position.GetVector(),
                locator,
                i_vertex)
            if vertex.IsLocated():
                self.vertices.push_back(vertex)

        # create temporal controller
        temporal_controller_settings = parameters["output_time_settings"]
        temporal_controller_settings.AddString("model_part_name", model_part.FullName())
        temporal_controller = KratosMultiphysics.TemporalController(model, temporal_controller_settings)

        # create the aggregated operation with hdf5 file settings
        operations = AggregatedControlledOperations(model_part, parameters["file_settings"])

        # adding one time mesh output
        operations.AddControlledOperation(ControlledOperation(VertexCoordinateOutput, parameters["model_part_output_settings"], SingleTimeController(temporal_controller)))

        # now adding temporal outputs.
        operations.AddControlledOperation(ControlledOperation(NodalSolutionStepDataOutput, parameters["nodal_solution_step_data_settings"], temporal_controller))
        operations.AddControlledOperation(ControlledOperation(NodalDataValueOutput, parameters["nodal_data_value_settings"], temporal_controller))
        operations.AddControlledOperation(ControlledOperation(NodalFlagValueOutput, parameters["nodal_flag_value_settings"], temporal_controller))
        operations.AddControlledOperation(ControlledOperation(ElementDataValueOutput, parameters["element_data_value_settings"], temporal_controller))
        operations.AddControlledOperation(ControlledOperation(ElementGaussPointOutput, parameters["element_gauss_point_value_settings"], temporal_controller))
        operations.AddControlledOperation(ControlledOperation(ElementFlagValueOutput, parameters["element_flag_value_settings"], temporal_controller))
        operations.AddControlledOperation(ControlledOperation(ConditionDataValueOutput, parameters["condition_data_value_settings"], temporal_controller))
        operations.AddControlledOperation(ControlledOperation(ConditionGaussPointOutput, parameters["condition_gauss_point_value_settings"], temporal_controller))
        operations.AddControlledOperation(ControlledOperation(ConditionFlagValueOutput, parameters["condition_flag_value_settings"], temporal_controller))

        # now add all operations to PrintOutput method
        self.AddPrintOutput(operations)

    def ExecuteInitialize(self) -> None:
        # generate the vertices


        return super().ExecuteInitialize()