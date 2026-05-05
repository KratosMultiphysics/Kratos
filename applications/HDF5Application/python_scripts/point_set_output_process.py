# Core imports
import KratosMultiphysics

# HDF5 imports
import KratosMultiphysics.HDF5Application as HDF5Application
from KratosMultiphysics.HDF5Application.core.file_io import OpenHDF5File
from KratosMultiphysics.HDF5Application.core.utils import EvaluatePattern


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

        # TODO: This whole process will be re-organized in coming PRs.
        self.is_historical = parameters["historical_value"].GetBool()

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
                if self.is_historical:
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
