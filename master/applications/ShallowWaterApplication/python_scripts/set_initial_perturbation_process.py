import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return InitialPerturbationProcess(Model, settings["Parameters"])

## This process sets the value of a scalar variable using the AssignScalarVariableProcess.
class InitialPerturbationProcess(KM.Process):

    def __init__(self, Model, settings):

        KM.Process.__init__(self)

        default_settings = KM.Parameters("""
            {
                "model_part_name"            : "main_model_part",
                "interval"                   : [0.0, 0.0],
                "source_type"                : "coordinates or model_part",
                "source_coordinates"         : [0.0, 0.0, 0.0],
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
        variable = KM.KratosGlobals.GetVariable(self.variable_name)

        # Creation of the parameters for the c++ process
        cpp_parameters = KM.Parameters("""{}""")
        cpp_parameters.AddValue("default_value", settings["default_value"])
        cpp_parameters.AddValue("distance_of_influence", settings["distance_of_influence"])
        cpp_parameters.AddValue("maximum_perturbation_value", settings["maximum_perturbation_value"])

        if settings["source_type"].GetString() == "coordinates":
            # retrieving the position of the source
            source_coordinates = settings["source_coordinates"].GetVector()
            if (source_coordinates.Size() != 3):
                raise Exception('The source_coordinates has to be provided with 3 coordinates! It has ', source_coordinates.Size())
            node = KM.Node(1, source_coordinates[0], source_coordinates[1], source_coordinates[2])
            # Construction of the process with one node
            self.perturbation_process = SW.ApplyPerturbationFunctionToScalar(self.model_part, node, variable, cpp_parameters)

        elif settings["source_type"].GetString() == "model_part":
            # Construction of the process with a sub model part
            source_model_part = Model[settings["source_model_part_name"].GetString()]
            self.perturbation_process = SW.ApplyPerturbationFunctionToScalar(self.model_part, source_model_part.Nodes, variable, cpp_parameters)

        else:
            raise Exception("InitialPerturbationProcess: unknown source type")

    def ExecuteInitialize(self):
        self.perturbation_process.Execute()
        if self.variable_name == "HEIGHT":
            SW.ShallowWaterUtilities().ComputeFreeSurfaceElevation(self.model_part)
        elif self.variable_name == "FREE_SURFACE_ELEVATION":
            SW.ShallowWaterUtilities().ComputeHeightFromFreeSurface(self.model_part)
