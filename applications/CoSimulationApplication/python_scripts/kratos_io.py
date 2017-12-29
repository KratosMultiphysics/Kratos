import KratosMultiphysics
# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("CoSimulationApplication")
import KratosMultiphysics.CoSimulationApplication as CoSimApp

# Other imports
import os

def CreateIo(custom_settings):
    return KratosIo(custom_settings)


class KratosIo(CoSimApp.CoSimulationBaseIo):
    def __init__(self, custom_settings):

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
            {
                "type":"kratos_io",
                "settings":{
                }
            }
        """)

        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        super(KratosIo, self).__init__(custom_settings)

        ### Constructing the IO for this solver