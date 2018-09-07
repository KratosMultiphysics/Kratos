from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics as KratosMultiphysics
import KratosMultiphysics.MeshingApplication as MeshingApplication

KratosMultiphysics.CheckForPreviousImport()

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return MultiscaleProcess(Model, settings["Parameters"])

class MultiscaleProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):

        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""
        {
            "main_model_part_name"            : "MainModelPart",
            "current_subscale"                : 0,
            "maximum_number_of_subscales"     : 4,    
            "echo_level"                      : 0,
            "advanced_configuration"          : {
                "echo_level"                      : 0,
                "number_of_divisions_at_subscale" : 2,
                "subscale_interface_name"         : "refined_interface",      
                "subscale_boundary_condition"     : "Condition2D2N"
            }
        }
        """)

        # Overwrite the default settings with user-provided parameters
        self.settings = settings
        self.settings.RecursivelyValidateAndAssignDefaults(default_parameters)

        self.model = Model

        self.current_subscale = self.settings['current_subscale'].GetInt()
        self.maximum_number_of_subscales = self.settings['maximum_number_of_subscales'].GetInt()

        self.coarse_model_part_name = self.settings['main_model_part_name'].GetString()
        if (self.current_subscale > 0):
            self.coarse_model_part_name += '_' + str(self.current_subscale)

        if (self.current_subscale < self.maximum_number_of_subscales):
            self._InitializeRefinedModelPart()

    def _InitializeRefinedModelPart(self):
        self.refined_model_part_name = self.settings['main_model_part_name'].GetString() + '_' + str(self.current_subscale + 1)
        coarse_model_part = self.model[self.coarse_model_part_name]
        refined_model_part = KratosMultiphysics.ModelPart(self.refined_model_part_name)
        self.model.AddModelPart(refined_model_part)

        # Create the new subscale process
        self.subscales_utility = MeshingApplication.MultiScaleRefiningProcess(coarse_model_part, refined_model_part, self.settings["advanced_configuration"])
