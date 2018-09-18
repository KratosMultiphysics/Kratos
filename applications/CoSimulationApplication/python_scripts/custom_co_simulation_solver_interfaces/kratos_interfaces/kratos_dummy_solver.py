from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from co_simulation_base_solver import CoSimulationBaseSolver
import co_simulation_tools as tools

def Create(settings):
    return KratosDummyCoSimulationSolver(settings)

class KratosDummyCoSimulationSolver(CoSimulationBaseSolver):

    def __init__(self, custom_settings):
        default_settings = KratosMultiphysics.Parameters("""
        {
            "name": "dummy",
            "type": "",
            "geometry_list": [],
            "io_type": "kratos_io",
            "settings":{
                "input_file":""
            }
        }""")
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)
        self.name = self.settings['name'].GetString()
        self.geometry = {}

        self.project_parameters_file_name = self.settings['settings']['input_file'].GetString()
        self.geometry = {}


    def Initialize(self):
        pass

    def Finalize(self):
        pass

    def SolveTimeStep(self):
        pass

    def InitializeTimeStep(self):
        # This function is not needed in this solver
        pass
    def FinalizeTimeStep(self):
        # This function is not needed in this solver
        pass

    def ImportData(self, DataName, FromClient):
        pass
    def ImportMesh(self, MeshName, FromClient):
        pass

    def ExportData(self, DataName, ToClient):
        pass
    def ExportMesh(self, MeshName, ToClient):
        pass

    def MakeDataAvailable(self, DataName, ToClient):
        pass
    def MakeMeshAvailable(self, MeshName, ToClient):
        pass

    ###########################################################################