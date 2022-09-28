from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library

import KratosMultiphysics
import KratosMultiphysics.ConstitutiveModelsApplication as KratosModels

#from multiprocessing import Pool

def Factory(settings, Model):
    if( not isinstance(settings,KratosMultiphysics.Parameters) ):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return NonLocalPlasticityProcess(Model, settings["Parameters"])


class NonLocalPlasticityProcess(KratosMultiphysics.Process):
    #
    def __init__(self, Model, custom_settings ):

        KratosMultiphysics.Process.__init__(self)

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "echo_level"            : 1,
            "model_part_name"       : "Main_Domain",
            "characteristic_length":  0.0
        }
        """)

        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = Model
        self.model_part_name = self.settings["model_part_name"].GetString()

        self.echo_level        = self.settings["echo_level"].GetInt()

        self.model = Model


    def ExecuteInitialize(self):

        self.main_model_part = self.model[ self.model_part_name]

        self.non_local_cxx_process = KratosModels.NonLocalPlasticityProcess( self.main_model_part, self.settings)

    def ExecuteInitializeSolutionStep(self):

        self.non_local_cxx_process.Execute()


