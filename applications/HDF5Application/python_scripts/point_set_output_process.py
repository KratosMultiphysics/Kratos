import KratosMultiphysics
import KratosMultiphysics.HDF5Application as HDF5Application


def Factory(parameters: KratosMultiphysics.Parameters,
            model: KratosMultiphysics.Model):
    if not isinstance(parameters, KratosMultiphysics.Parameters):
        raise Exception("expecting input parameters of type KratosMultiphysics.Parameters, but got {}".format(type(parameters)))
    return PointSetOutputProcess(model, parameters)


class PointSetOutputProcess(KratosMultiphysics.Process):

    def __init__(self, model: KratosMultiphysics.Model, parameters: KratosMultiphysics.Parameters):
        KratosMultiphysics.Process.__init__(self)
        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())        

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
        self.locator = HDF5Application.BruteForcePointLocatorAdaptor(self.model_part,
                                                                     search_configuration,
                                                                     search_tolerance)

        self.interval_utility = KratosMultiphysics.IntervalUtility(parameters)

        # Initialize file
        file_parameters = KratosMultiphysics.Parameters()
        file_parameters.AddValue("file_name", parameters["file_path"])

        # TODO: detect mpi run
        is_mpi_run = False

        if is_mpi_run:
            self.file = HDF5Application.HDF5FileParallel(file_parameters)
        else:
            self.file = HDF5Application.HDF5FileSerial(file_parameters)

        self.prefix = parameters["prefix"].GetString()
        self.isHistorical = parameters["historical_value"].GetBool()

        # Create vertices
        self.vertices = HDF5Application.VertexContainer()
        for i_vertex, position in enumerate(parameters["positions"]):
            self.vertices.push_back(HDF5Application.Vertex.MakeShared(
                position.GetVector(),
                self.locator,
                i_vertex,
                self.isHistorical))


    def ExecuteInitialize(self):
        io_parameters = KratosMultiphysics.Parameters()
        io_parameters.AddString("prefix", self.prefix)
        HDF5Application.VertexContainerIO(io_parameters, self.file).WriteCoordinatesAndIDs(self.vertices)


    def ExecuteFinalizeSolutionStep(self):
        time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        if self.interval_utility.IsInInterval(time):
            step_id = self.model_part.ProcessInfo[KratosMultiphysics.STEP]
            group_name = "step_{}".format(step_id)

            io_parameters = KratosMultiphysics.Parameters()
            io_parameters.AddString("prefix", self.prefix)
            io_parameters.AddString("variables_path", "/" + group_name)
            io_parameters.AddValue("list_of_variables", self.variables)

            HDF5Application.VertexContainerIO(io_parameters, self.file).WriteVariables(self.vertices)


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
            "file_path"             : "",
            "prefix"                : "/point_set_output"
        }""")