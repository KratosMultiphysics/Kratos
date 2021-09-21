import KratosMultiphysics
import KratosMultiphysics.HDF5Application as HDF5Application


def Factory(parameters: KratosMultiphysics.Parameters,
            model: KratosMultiphysics.Model):
    if not isinstance(parameters, KratosMultiphysics.Parameters):
        raise Exception("expecting input parameters of type KratosMultiphysics.Parameters, but got {}".format(type(parameters)))
    return PointSetOutputProcess(model, parameters)


class StepPattern:

    def __init__(self,
                 pattern: str,
                 time_format = "",
                 step_format = ""):
        self.pattern = pattern
        self.time_format = time_format
        self.step_format = step_format

        self.time_tag = "<time>"
        self.step_tag = "<step>"

    
    def Substitute(self, model_part: KratosMultiphysics.ModelPart):
        string = self.pattern

        if self.time_tag in string:
            time = model_part.ProcessInfo[KratosMultiphysics.TIME]
            string = string.replace(self.time_tag, format(time, self.time_format))

        if self.step_tag in string:
            step = model_part.ProcessInfo[KratosMultiphysics.STEP]
            string = string.replace(self.step_tag, format(step, self.step_format))

        return string



class PointSetOutputProcess(KratosMultiphysics.Process):

    def __init__(self, model: KratosMultiphysics.Model, parameters: KratosMultiphysics.Parameters):
        KratosMultiphysics.Process.__init__(self)
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

        self.file_parameters = parameters["file_parameters"].Clone()

        communicator = KratosMultiphysics.ParallelEnvironment.GetDefaultDataCommunicator()
        self.isMPIRun = 1 < communicator.Size()

        self.group_prefix = parameters["group_prefix"].GetString()
        self.step_group_pattern = StepPattern(parameters["step_group_pattern"].GetString())
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


    def ExecuteInitialize(self):
        io_parameters = KratosMultiphysics.Parameters()
        io_parameters.AddString("group_prefix", self.group_prefix)
        with PointSetOutputProcess.OpenHDF5File(self.__GetCurrentFileParameters(), self.isMPIRun) as file:
            HDF5Application.VertexContainerIO(io_parameters, file).WriteCoordinatesAndIDs(self.vertices)


    def ExecuteFinalizeSolutionStep(self):
        time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        if self.interval_utility.IsInInterval(time):
            group_name = self.step_group_pattern.Substitute(self.model_part)

            io_parameters = KratosMultiphysics.Parameters()
            io_parameters.AddString("group_prefix", self.group_prefix)
            io_parameters.AddString("variables_path", "/" + group_name)
            io_parameters.AddValue("list_of_variables", self.variables)

            with PointSetOutputProcess.OpenHDF5File(self.__GetCurrentFileParameters(), self.isMPIRun) as file:
                HDF5Application.VertexContainerIO(io_parameters, file).WriteVariables(self.vertices)


    @staticmethod
    def GetDefaultParameters():
        return KratosMultiphysics.Parameters("""{
            "model_part_name"       : "",
            "interval"              : [0.0, "End"],
            "positions"             : [[]],
            "output_variables"      : [],
            "historical_value"      : true,
            "search_configuration"  : "initial",
            "search_tolerance"      : 1e-6,
            "group_prefix"          : "/point_set_output",
            "step_group_pattern"    : "step_<step>",
            "file_parameters"       : {
                "file_name"         : "",
                "file_access_mode"  : "read_write",
                "file_driver"       : "sec2",
                "echo_level"        : 0
            }
        }""")


    def __GetCurrentFileParameters(self):
        parameters = self.file_parameters.Clone()
        file_name = StepPattern(
            parameters["file_name"].GetString(),
            time_format=".4f",
            step_format="").Substitute(self.model_part)

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