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
            "model_part_name"                 : "MainModelPart",
            "maximum_number_of_subscales      : 4,    
            "echo_level"                      : 0,
            "number_of_divisions_at_subscale" : 2,
            "refining_interface_model_part"   : "refining_interface",      
            "refining_boundary_condition"     : "Condition2D2N"
        }
        """)

        # Overwrite the default settings with user-provided parameters
        self.settings = settings
        self.settings.RecursivelyValidateAndAssignDefaults(default_parameters)

        self.model = Model

        self.maximum_number_of_subscales = self.settings['maximum_number_of_subscales'].GetInt()
        self.number_of_divisions_at_subscale = self.settings['number_of_divisions_at_subscale'].GetInt()