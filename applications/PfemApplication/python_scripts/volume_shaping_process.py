from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.PfemApplication as KratosPfem

def Factory(custom_settings, Model):
    if( not isinstance(custom_settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return VolumeShapingProcess(Model, custom_settings["Parameters"])

## All the processes python should be derived from "Process"
class VolumeShapingProcess(KratosMultiphysics.Process):
    def __init__(self, Model, custom_settings ):
        KratosMultiphysics.Process.__init__(self)

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
             "model_part_name": "main_domain",
             "variable_name" : "VOLUME_WEAR",
             "flags_list" : [],
             "properties" : {}
        }
        """)

        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        self.model = Model

    #
    def GetVariables(self):
        nodal_variables = [self.settings["variable_name"].GetString()]
        return nodal_variables

    #
    def ExecuteInitialize(self):

        # set model part
        self.model_part = self.model[self.settings["model_part_name"].GetString()]

        echo_level = 0;

        params = KratosMultiphysics.Parameters("{}")
        params.AddValue("variable_name", self.settings["variable_name"])
        params.AddValue("flags_list", self.settings["flags_list"])
        params.AddValue("properties", self.settings["properties"])

        self.VolumeShapingProcess = KratosPfem.VolumeShapingProcess(self.model_part, params)

        self.VolumeShapingProcess.ExecuteInitialize()


    def ExecuteFinalizeSolutionStep(self):

        self.VolumeShapingProcess.ExecuteFinalizeSolutionStep()
