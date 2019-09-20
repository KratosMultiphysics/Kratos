from __future__ import print_function, absolute_import, division

from KratosMultiphysics.python_solver import PythonSolver
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

def CreateTurbulenceModel(model, settings, is_distributed = False):
    if not CheckIfApplicationsAvailable("RANSModellingApplication"):
        msg = "Using a turbulence model requires the RANSModellingApplication. "
        msg += "Please re-install/re-compile with RANSModellingApplication."
        raise Exception(msg)

    from KratosMultiphysics.RANSModellingApplication.turbulence_model_factory import Factory
    return Factory(settings, model, is_distributed)

class TurbulenceModelSolver(PythonSolver):
    def AddVariables(self):
        msg = "Calling the base TurbulenceModelSolver class AddVariables method."
        msg += " Please override it in the derrived class."
        raise Exception(msg)

    def AddDofs(self):
        msg = "Calling the base TurbulenceModelSolver class AddDofs method."
        msg += " Please override it in the derrived class."
        raise Exception(msg)

    def PrepareModelPart(self):
        msg = "Calling the base TurbulenceModelSolver class PrepareModelPart method."
        msg += " Please override it in the derrived class."
        raise Exception(msg)

    def GetTurbulenceSolvingProcess(self):
        msg = "Calling the base TurbulenceModelSolver class GetTurbulenceSolvingProcess method."
        msg += " Please override it in the derrived class to return a KratosMultiphysics.Process."
        raise Exception(msg)

    def Initialize(self):
        msg = "Calling the base TurbulenceModelSolver class Initialize method."
        msg += " Please override it in the derrived class."
        raise Exception(msg)

    def Check(self):
        msg = "Calling the base TurbulenceModelSolver class Check method."
        msg += " Please override it in the derrived class."
        raise Exception(msg)

    def InitializeSolutionStep(self):
        msg = "Calling the base TurbulenceModelSolver class InitializeSolutionStep method."
        msg += " Please override it in the derrived class."
        raise Exception(msg)

    def FinalizeSolutionStep(self):
        msg = "Calling the base TurbulenceModelSolver class FinalizeSolutionStep method."
        msg += " Please override it in the derrived class."
        raise Exception(msg)

    def GetFluidVelocityPressureConditionName(self):
        msg = "Calling the base TurbulenceModelSolver class GetFluidVelocityPressureConditionName method."
        msg += " Please override it in the derrived class to return a condition name."
        raise Exception(msg)

    def SetCommunicator(self, epetra_communicator):
        msg = "Calling the base TurbulenceModelSolver class SetCommunicator method."
        msg += " Please override it in the derrived class to set the epetra_communicator for RANS model parts"
        raise Exception(msg)

    def Finalize(self):
        msg = "Calling the base TurbulenceModelSolver class SetCommunicator method."
        msg += " Please override it in the derrived class"
        raise Exception(msg)