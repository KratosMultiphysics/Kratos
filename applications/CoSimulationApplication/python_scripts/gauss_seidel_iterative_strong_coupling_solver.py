import KratosMultiphysics
# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("CoSimulationApplication")
import KratosMultiphysics.CoSimulationApplication as CoSimApp


# Other imports
import os

def CreateSolver(app1, app2):
    return IterativeStrongCouplingSolver(app1, app2)

class IterativeStrongCouplingSolver(CoSimApp.CoSimulationBaseCouplingStrategy):
    def __init__(self, app1, app2):
        super(IterativeStrongCouplingSolver, self).__init__(app1, app2)
        self.AppOne = app1
        self.AppTwo = app2
    def InitializeSolutionStep(self):
        self.AppOne.ImportModelPart()
        self.AppOne.InitializeSolutionStep()
        self.AppTwo.InitializeSolutionStep()
    def Predict(self):
        print("predicting ... !!")
    def Initialize(self):
        pass
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
        self.Predict()
        self.TransferDataField()
        print("solving .................. !!")
        self.AppOne.Solve()
        self.AppTwo.Solve()