try:
    import KratosMultiphysics
    # Check that applications were imported in the main script
    KratosMultiphysics.CheckRegisteredApplications("CoSimulationApplication")
    import KratosMultiphysics.CoSimulationApplication as CoSimApp
except:
except ModuleNotFoundError:
        print(tools.bcolors.FAIL + 'KRATOS is not available ! Please ensure that Kratos is available for usage !'+ tools.bcolors.ENDC)
        exit()

from co_simulation_base_io import CoSimulationBaseIO
import co_simulation_tools as tools

# Other imports
import os

def Create(custom_settings):
    return KratosIo(custom_settings)

class KratosIo(CoSimulationBaseIO):
    def __init__(self, custom_settings):

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
            {
                "type":"kratos_io",
                "settings":{
                }
            }
        """)
        self.settings = tools.ValidateAndAssignInputParameters(default_settings, custom_settings, False)

        super(KratosIo, self).__init__(custom_settings)
        ### Constructing the IO for this solver