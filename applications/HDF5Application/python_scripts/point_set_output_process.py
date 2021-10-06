import KratosMultiphysics
import KratosMultiphysics.HDF5Application as HDF5Application
from KratosMultiphysics.HDF5Application.core.operations.model_part import Prefix


def Factory(parameters: KratosMultiphysics.Parameters,
            model: KratosMultiphysics.Model):
    if not isinstance(parameters, KratosMultiphysics.Parameters):
        raise Exception("expecting input parameters of type KratosMultiphysics.Parameters, but got {}".format(type(parameters)))
    return PointSetOutputProcess(model, parameters)



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
        self.output_counter = 0

        self.file_parameters = parameters["file_parameters"].Clone()

        communicator = KratosMultiphysics.ParallelEnvironment.GetDefaultDataCommunicator()
        self.isMPIRun = 1 < communicator.Size()

        self.coordinates_prefix_pattern = parameters["coordinates_prefix"].GetString()
        self.variables_prefix_pattern = parameters["variables_prefix"].GetString()
        isHistorical = parameters["historical_value"].GetBool()

        # Create vertices
        self.vertices = HDF5Application.VertexContainer()
        for i_vertex, position in enumerate(parameters["positions"]):
            vertex = HDF5Application.Vertex.MakeShared(
                position.GetVector(),
                locator,
                i_vertex,
                isHistorical)
            if vertex.IsLocated():
                self.vertices.push_back(vertex)


    def IsOutputStep(self):
        time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        return self.interval_utility.IsInInterval(time) and (self.output_counter % self.output_frequency == 0)


    def ExecuteInitialize(self):
        coordinates_path = Prefix(
            self.coordinates_prefix_pattern,
            self.model_part)

        io_parameters = KratosMultiphysics.Parameters("""{
            "prefix" : "",
            "write_vertex_ids" : true
        }""")
        io_parameters["prefix"].SetString(coordinates_path)
        with PointSetOutputProcess.OpenHDF5File(self.__GetCurrentFileParameters(), self.isMPIRun) as file:
            HDF5Application.VertexContainerCoordinateIO(io_parameters, file).Write(self.vertices)


    def PrintOutput(self):
        prefix = Prefix(
            self.variables_prefix_pattern,
            self.model_part)

        io_parameters = KratosMultiphysics.Parameters()
        io_parameters.AddString("prefix", prefix)
        io_parameters.AddValue("list_of_variables", self.variables.Clone())

        with PointSetOutputProcess.OpenHDF5File(self.__GetCurrentFileParameters(), self.isMPIRun) as file:
            HDF5Application.VertexContainerVariableIO(io_parameters, file).Write(self.vertices)


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
                "file_access_mode"  : "read_write",
                "file_driver"       : "sec2",
                "echo_level"        : 0
            }
        }""")


    def __GetCurrentFileParameters(self):
        parameters = self.file_parameters.Clone()
        file_name = Prefix(
            parameters["file_name"].GetString(),
            self.model_part)

        parameters["file_name"].SetString(file_name)
        return parameters


    class OpenHDF5File(object):

        def __init__(self,
                     file_parameters: KratosMultiphysics.Parameters,
                     isMPIRun: bool):
            parameters = file_parameters.Clone()
            if isMPIRun:
                parameters["file_driver"].SetString("mpio")
                self.file =  HDF5Application.HDF5FileParallel(parameters)
            else:
                self.file =  HDF5Application.HDF5FileSerial(parameters)


        def __enter__(self):
            return self.file


        def __exit__(self, exit_type, exit_value, exit_traceback):
            # HDF5::File has RAII, so this is the best we can do to close it
            self.file = None