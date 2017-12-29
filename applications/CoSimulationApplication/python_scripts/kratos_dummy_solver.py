import KratosMultiphysics
# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("CoSimulationApplication")
import KratosMultiphysics.CoSimulationApplication as CoSimApp

# Other imports
import os

def CreateSolver(custom_settings):
    return KratosDummySolver(custom_settings)

class KratosDummySolver(CoSimApp.CoSimulationBaseApplication):
    def __init__(self, custom_settings):

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "type": "kratos_dummy_solver",
            "name": "dummy1",
            "io": {
                "type":"kratos_io",
                "settings":{
                }
            },   
            "output_data_field_list": ["DISPLACEMENT"],
            "input_data_field_list": ["POINT_LOAD"],
            "settings":{
            },
            "output_configuration": {
                "file_label"          : "step",
                "output_control_type" : "step",
                "output_frequency"    : 1,
                "nodal_results"       : ["DISPLACEMENT","REACTION"],
                "elemental_results"   : []   
            }
        }""")

        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        super(KratosDummySolver, self).__init__(custom_settings)        

        ### Constructing the IO for this solver
        iOModule = __import__(self.settings['io']['type'].GetString())
        self.io = iOModule.CreateIo(self.settings['io'])

    def InitializeSolutionStep(self):
        pass
    def Predict(self):
        pass
    def Initialize(self):
        print("KratosDummySolver initializing .................. !!")
    def Clear(self):
        pass
    def IsConverged(self):
        pass
    def CalculateOutputData(self):
        pass
    def FinalizeSolutionStep(self):
        pass
    def SolveSolutionStep(self):
        pass
    def SetEchoLevel(self):
        pass
    def GetEchoLevel(self):
        pass
    def GetResidualNorm(self):
        pass
    def Solve(self):
        print("KratosDummySolver solving .................. !!")
    def GetModelPart(self):
        return KratosMultiphysics.ModelPart('app')
    def SynchronizeInputData(self):
        self.io.SynchronizeInputData()
    def SynchronizeOutputData(self):
        self.io.SynchronizeOutputData()
    