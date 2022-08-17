# Importing the Kratos Library
import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

def Factory(settings, model):
    if not isinstance(settings, Kratos.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    if not isinstance(model, Kratos.Model):
        raise Exception("expected input shall be a model object")

    return FluidUtilitiesProcess(model, settings["Parameters"])


class CFLUtility:
    def __init__(self, model_part, params):
        self.model_part = model_part

        default_settings = Kratos.Parameters("""
            {
                "utility_type": "cfl",
                "echo_level"  : 0
            }
            """)
        params.ValidateAndAssignDefaults(default_settings)
        self.echo_level = params["echo_level"].GetInt()

    def Execute(self):
        KratosCFD.FluidCharacteristicNumbersUtilities.CalculateLocalCFL(self.model_part)
        if self.echo_level > 0:
            Kratos.Logger.PrintInfo("CFLUtility", "Calculated CFL numebers of elements in {:s}.".format(self.model_part.FullName()))


class FluidUtilitiesProcess(Kratos.Process):
    """
    A class responsible for calculating values defined in utilities.
    """

    def __init__(self, model, params):
        Kratos.Process.__init__(self)

        default_settings = Kratos.Parameters("""
            {
                "model_part_name" : "PLEASE_PROVIDE_A_MODEL_PART_NAME",
                "interval"        : [0.0, "End"],
                "echo_level"      : 0,
                "utility_settings": {
                    "utility_type": "cfl"
                }
            }
            """)

        self.interval_utility = Kratos.IntervalUtility(params)
        params.ValidateAndAssignDefaults(default_settings)

        # getting the ModelPart from the Model
        self.model_part = model[params["model_part_name"].GetString()]

        if params["utility_settings"].Has("utility_type"):
            utility_type = params["utility_settings"]["utility_type"].GetString()
            if utility_type == "cfl":
                self.utility = CFLUtility(self.model_part, params["utility_settings"])
            else:
                raise RuntimeError("Unsupported utility type requested. [ utility_type = {:s}]. Supported utility types are:\n\t\t \"cfl\"".format(utility_type))
        else:
            raise RuntimeError("Please provide \"utiltiy_type\" string setting under \"utility_settings\"")

    def ExecuteFinalizeSolutionStep(self):
        if self.interval_utility.IsInInterval(self.model_part.ProcessInfo[Kratos.TIME]):
            self.utility.Execute()
