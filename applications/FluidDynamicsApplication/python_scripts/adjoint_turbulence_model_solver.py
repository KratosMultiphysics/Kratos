from __future__ import print_function, absolute_import, division

from KratosMultiphysics.python_solver import PythonSolver
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

def CreateAdjointTurbulenceModel(model, settings):
    if not CheckIfApplicationsAvailable("RANSModellingApplication"):
        msg = "Using a turbulence model requires the RANSModellingApplication. "
        msg += "Please re-install/re-compile with RANSModellingApplication."
        raise Exception(msg)

    from KratosMultiphysics.RANSModellingApplication.adjoint_turbulence_model_factory import Factory
    return Factory(settings, model)

class AdjointTurbulenceModelSolver(PythonSolver):
    def AddVariables(self):
        msg = "Calling the base AdjointTurbulenceModelSolver class AddVariables method."
        msg += " Please override it in the derrived class."
        raise Exception(msg)

    def AddDofs(self):
        msg = "Calling the base AdjointTurbulenceModelSolver class AddDofs method."
        msg += " Please override it in the derrived class."
        raise Exception(msg)

    def GetAdjointElementName(self):
        msg = "Calling the base AdjointTurbulenceModelSolver class GetAdjointElementName method."
        msg += " Please override it in the derrived class to return an element name."
        raise Exception(msg)

    def GetAdjointConditionName(self):
        msg = "Calling the base AdjointTurbulenceModelSolver class GetAdjointConditionName method."
        msg += " Please override it in the derrived class to return an element name."
        raise Exception(msg)

    def GetAdjointResposeFunction(self):
        msg = "Calling the base AdjointTurbulenceModelSolver class GetAdjointResposeFunction method."
        msg += " Please override it in the derrived class to return a response function."
        raise Exception(msg)

    def Initialize(self):
        msg = "Calling the base AdjointTurbulenceModelSolver class Initialize method."
        msg += " Please override it in the derrived class."
        raise Exception(msg)

    def Check(self):
        msg = "Calling the base AdjointTurbulenceModelSolver class Check method."
        msg += " Please override it in the derrived class."
        raise Exception(msg)

    def InitializeSolutionStep(self):
        msg = "Calling the base AdjointTurbulenceModelSolver class InitializeSolutionStep method."
        msg += " Please override it in the derrived class."
        raise Exception(msg)

    def FinalizeSolutionStep(self):
        msg = "Calling the base AdjointTurbulenceModelSolver class FinalizeSolutionStep method."
        msg += " Please override it in the derrived class."
        raise Exception(msg)

    def Finalize(self):
        msg = "Calling the base TurbulenceModelSolver class SetCommunicator method."
        msg += " Please override it in the derrived class"
        raise Exception(msg)