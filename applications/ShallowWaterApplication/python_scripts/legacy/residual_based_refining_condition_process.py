# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW
import KratosMultiphysics.MeshingApplication as MSH

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ResidualBasedRefiningConditionProcess(Model, settings["Parameters"])

class ResidualBasedRefiningConditionProcess(KM.Process):
    def __init__(self, model, settings ):

        ## Settings string in json format
        default_parameters = KM.Parameters("""
        {
            "model_part_name"           : "model_part",
            "error_variable"            : "RESIDUAL_NORM",
            "variable_threshold"        : 1e-3,
            "increase_threshold"        : true,
            "only_refine_wet_domain"    : true
        }
        """)

        # Overwrite the default settings with user-provided parameters
        settings.ValidateAndAssignDefaults(default_parameters)

        model_part = model[settings["model_part_name"].GetString()]

        # Creation of the parameters for the c++ process
        criteria_parameters = KM.Parameters("""{}""")
        criteria_parameters.AddValue("error_variable", settings["error_variable"])
        criteria_parameters.AddValue("only_refine_wet_domain", settings["only_refine_wet_domain"])
        criteria_parameters.AddValue("variable_threshold", settings["variable_threshold"])

        if settings["increase_threshold"].GetBool():
            variable_threshold = settings["variable_threshold"].GetDouble()
            next_subscale_index = model_part.ProcessInfo[MSH.SUBSCALE_INDEX] + 1
            variable_threshold *= next_subscale_index
            criteria_parameters["variable_threshold"].SetDouble(variable_threshold)

        self.process = SW.ElementalRefiningCriteriaProcess(model_part, criteria_parameters)

    def Execute(self):
        self.process.Execute()
