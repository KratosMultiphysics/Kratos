import KratosMultiphysics
import KratosMultiphysics.ShallowWaterApplication as Shallow

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return InitialPerturbationProcess(Model, settings["Parameters"])

## This process sets the value of a scalar variable using the AssignScalarVariableProcess.
class InitialPerturbationProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "model_part_name"            : "main_model_part",
                "source_type"                : "point or model_part",
                "source_point_coordinates"   : [0.0, 0.0, 0.0],
                "source_model_part_name"     : "main_model_part.sub_model_part",
                "variable_name"              : "FREE_SURFACE_ELEVATION",
                "default_value"              : 0.0,
                "distance_of_influence"      : 1.0,
                "maximum_perturbation_value" : 1.0
            }
            """
            )
        settings.ValidateAndAssignDefaults(default_settings)

        self.variable_name = settings["variable_name"].GetString()
        self.model_part = Model[settings["model_part_name"].GetString()]

        # Creation of the parameters for the c++ process
        cpp_parameters = KratosMultiphysics.Parameters("""{}""")
        cpp_parameters.AddValue("variable_name", settings["variable_name"])
        cpp_parameters.AddValue("default_value", settings["default_value"])
        cpp_parameters.AddValue("distance_of_influence", settings["distance_of_influence"])
        cpp_parameters.AddValue("maximum_perturbation_value", settings["maximum_perturbation_value"])

        if settings["source_type"].GetString() == "point":
            # retrieving the position of the point
            point_position = settings["source_point_coordinates"].GetVector()
            if (point_position.Size() != 3):
                raise Exception('The source_point_coordinates has to be provided with 3 coordinates! It has ', point_position.Size())
            node = KratosMultiphysics.Node(1, point_position[0], point_position[1], point_position[2])
            # Construction of the process with one node
            self.perturbation_process = Shallow.InitialPerturbationProcess(self.model_part, node, cpp_parameters)

        elif settings["source_type"].GetString() == "model_part":
            # Construction of the process with a sub model part
            source_model_part = Model[settings["source_model_part_name"].GetString()]
            self.perturbation_process = Shallow.InitialPerturbationProcess(self.model_part, source_model_part.Nodes, cpp_parameters)

        else:
            raise Exception("InitialPerturbationProcess: unknown source type")

        self.variables_utility = Shallow.ShallowWaterVariablesUtility(self.model_part)

    def ExecuteInitialize(self):
        self.perturbation_process.Execute()
        if self.variable_name == "HEIGHT":
            self.variables_utility.ComputeFreeSurfaceElevation()
        elif self.variable_name == "FREE_SURFACE_ELEVATION":
            self.variables_utility.ComputeHeightFromFreeSurface()
