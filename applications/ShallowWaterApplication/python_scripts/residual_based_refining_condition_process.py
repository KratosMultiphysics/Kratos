from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics as KratosMultiphysics
import KratosMultiphysics.ShallowWaterApplication as Shallow

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ResidualBasedRefiningConditionProcess(Model, settings["Parameters"])

class ResidualBasedRefiningConditionProcess(KratosMultiphysics.Process):
    def __init__(self, model, settings ):

        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""
        {
            "model_part_name"           : "model_part",
            "error_variable"            : "RESIDUAL_NORM",
            "variable_threshold"        : 1e-3,
            "only_refine_wet_domain"    : true
        }
        """)

        # Overwrite the default settings with user-provided parameters
        settings.ValidateAndAssignDefaults(default_parameters)

        self.model_part = model[settings["model_part_name"].GetString()]
        self.error_variable = getattr(KratosMultiphysics, settings["error_variable"].GetString())
        self.variable_threshold = settings["variable_threshold"].GetDouble()
        self.only_refine_wet_domain = settings["only_refine_wet_domain"].GetBool()

        self.process = Shallow.ElementalRefiningCriteriaProcess(self.model_part, settings)

    def Execute(self):
        self.process.Execute()
