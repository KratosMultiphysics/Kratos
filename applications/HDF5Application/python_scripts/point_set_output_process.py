# Core imports
import KratosMultiphysics

# HDF5 imports
import KratosMultiphysics.HDF5Application as HDF5Application
from KratosMultiphysics.HDF5Application.core.operations.aggregated_operations import ControlledOperation
from KratosMultiphysics.HDF5Application.core.operations.aggregated_operations import AggregatedControlledOperations
from KratosMultiphysics.HDF5Application.core.processes import HDF5OutputProcess
from KratosMultiphysics.HDF5Application.core.operations.model_part import VertexHistoricalValueOutput
from KratosMultiphysics.HDF5Application.core.operations.model_part import VertexNonHistoricalValueOutput
from KratosMultiphysics.HDF5Application.core.operations.model_part import VertexCoordinateOutput
from KratosMultiphysics.HDF5Application.core.controllers import SingleTimeController


def Factory(parameters: KratosMultiphysics.Parameters,
            model: KratosMultiphysics.Model):
    if not isinstance(parameters, KratosMultiphysics.Parameters):
        raise Exception("expecting input parameters of type KratosMultiphysics.Parameters, but got {}".format(type(parameters)))
    if not isinstance(model, KratosMultiphysics.Model):
        raise Exception("expecting input model of type KratosMultiphysics.Model, but got {}".format(type(model)))
    return PointSetOutputProcess(model, parameters["Parameters"])


class PointSetOutputProcess(HDF5OutputProcess):
    @classmethod
    def GetDefaultParameters(cls) -> KratosMultiphysics.Parameters:
        return KratosMultiphysics.Parameters("""
            {
                "model_part_name": "PLEASE_SPECIFY_MODEL_PART_NAME",
                "file_settings": {
                    "file_name"        : "",
                    "time_format"      : "0.4f",
                    "file_access_mode" : "truncate",
                    "max_files_to_keep": "unlimited",
                    "echo_level"       :  0
                },
                "output_time_settings": {
                    "output_control_type": "step",
                    "output_interval"    : 1.0,
                    "interval"           : [0.0, "End"]
                },
                "point_output_settings": {
                    "prefix"              : "/VertexData/Coordinates",
                    "time_format"         : "0.4f",
                    "positions"           : [[]],
                    "search_configuration": "initial",
                    "search_tolerance"    : 1e-6,
                    "custom_attributes"   : {}
                },
                "nodal_solution_step_data_settings": {
                    "prefix"           : "/VertexData/VertexSolutionStepData",
                    "list_of_variables": [],
                    "time_format"      : "0.4f",
                    "custom_attributes": {}
                },
                "nodal_data_value_settings": {
                    "prefix"           : "/VertexData/VertexDataValues",
                    "list_of_variables": [],
                    "time_format"      : "0.4f",
                    "custom_attributes": {}
                }
            }""")

    def __init__(self, model: KratosMultiphysics.Model, parameters: KratosMultiphysics.Parameters) -> None:
        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())
        self.model_part = model[parameters["model_part_name"].GetString()]
        super().__init__()

        point_output_settings = parameters["point_output_settings"]

        # create locator
        configuration_name = point_output_settings["search_configuration"].GetString()
        if configuration_name == "initial":
            search_configuration = KratosMultiphysics.Configuration.Initial
        elif configuration_name == "current":
            search_configuration = KratosMultiphysics.Configuration.Current
        else:
            raise RuntimeError("Invalid search configuration: '{}'".format(configuration_name))
        search_tolerance = point_output_settings["search_tolerance"].GetDouble()
        locator = HDF5Application.BruteForcePointLocatorAdaptor(self.model_part,
                                                                search_configuration,
                                                                search_tolerance)

        # create vertices
        self.vertices = HDF5Application.VertexContainer()
        for i_vertex, position in enumerate(point_output_settings["positions"].values()):
            vertex = HDF5Application.Vertex.MakeShared(
                position.GetVector(),
                locator,
                i_vertex)
            if vertex.IsLocated():
                self.vertices.push_back(vertex)

        # create temporal controller
        temporal_controller_settings = parameters["output_time_settings"]
        temporal_controller_settings.AddString("model_part_name", self.model_part.FullName())
        temporal_controller = KratosMultiphysics.OutputController(model, temporal_controller_settings)
        self.interval_utility = KratosMultiphysics.IntervalUtility(temporal_controller_settings)

        # create the aggregated operation with hdf5 file settings
        operations = AggregatedControlledOperations(self.model_part, self._GetValidatedParameters("file_settings", parameters))

        # adding one time mesh output
        operations.AddControlledOperation(ControlledOperation(VertexCoordinateOutput, self._GetOperationParameters("point_output_settings", parameters), SingleTimeController(temporal_controller), self.vertices))

        # now adding temporal outputs.
        operations.AddControlledOperation(ControlledOperation(VertexHistoricalValueOutput, self._GetOperationParameters("nodal_solution_step_data_settings", parameters), temporal_controller, self.vertices))
        operations.AddControlledOperation(ControlledOperation(VertexNonHistoricalValueOutput, self._GetOperationParameters("nodal_data_value_settings", parameters), temporal_controller, self.vertices))

        # now add all operations to PrintOutput method
        self.AddPrintOutput(operations)

    def IsOutputStep(self) -> bool:
        time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        return self.interval_utility.IsInInterval(time) and super().IsOutputStep()