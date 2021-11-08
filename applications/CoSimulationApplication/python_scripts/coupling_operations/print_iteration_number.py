# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupling_operation import CoSimulationCouplingOperation

class PrintIterationNumberOperation(CoSimulationCouplingOperation):
    # this is a dummy implementation, this be used with iterative coupled-solvers
    def Initialize(self):
        pass

    def Finalize(self):
        pass


    def InitializeSolutionStep(self):
        pass

    def FinalizeSolutionStep(self):
        pass


    def InitializeCouplingIteration(self):
        pass

    def FinalizeCouplingIteration(self):
        pass


    def Execute(self):
        pass


    def PrintInfo(self):
        pass

    def Check(self):
        pass

