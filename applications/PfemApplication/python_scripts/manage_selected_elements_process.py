from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.PfemApplication as KratosPfem

def Factory(custom_settings, Model):
    if( not isinstance(custom_settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ManageSelectedElementsProcess(Model, custom_settings["Parameters"])

class ManageSelectedElementsProcess(KratosMultiphysics.Process):
    def __init__(self, Model, custom_settings ):

        KratosMultiphysics.Process.__init__(self)

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
             "model_part_name": "MODEL_PART_NAME"
        }
        """)

        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        self.model = Model

    def ExecuteInitialize(self):

        # set model part
        self.model_part = self.model[self.settings["model_part_name"].GetString()]

        self.ManageElementsProcess = KratosPfem.ManageSelectedElementsProcess(self.model_part)

    def ExecuteBeforeSolutionLoop(self):

        self.ManageElementsProcess.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):

        self.ManageElementsProcess.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):

        self.ManageElementsProcess.ExecuteFinalizeSolutionStep()
