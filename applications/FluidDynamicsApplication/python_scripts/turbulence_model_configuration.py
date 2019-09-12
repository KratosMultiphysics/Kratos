from __future__ import print_function, absolute_import, division

import KratosMultiphysics as Kratos
from python_solver import PythonSolver
from kratos_utilities import CheckIfApplicationsAvailable

if CheckIfApplicationsAvailable("RANSModellingApplication"):
    import KratosMultiphysics.RANSModellingApplication as KratosRANS


def CreateTurbulenceModel(model, settings, parallel_type = "OpenMP"):
    if not CheckIfApplicationsAvailable("RANSModellingApplication"):
        msg = "Using a turbulence model requires the RANSModellingApplication. "
        msg += "Please re-install/re-compile with RANSModellingApplication."
        raise Exception(msg)

    from turbulence_model_factory import Factory
    return Factory(settings, model, parallel_type)

class TurbulenceModelConfiguration(PythonSolver):
    def AddVariables(self):
        msg = "Calling the base TurbulenceModelConfiguration class AddVariables method."
        msg += " Please override it in the derrived class."
        raise Exception(msg)

    def AddDofs(self):
        msg = "Calling the base TurbulenceModelConfiguration class AddDofs method."
        msg += " Please override it in the derrived class."
        raise Exception(msg)

    def PrepareModelPart(self):
        msg = "Calling the base TurbulenceModelConfiguration class PrepareModelPart method."
        msg += " Please override it in the derrived class."
        raise Exception(msg)

    def GetTurbulenceSolvingProcess(self):
        msg = "Calling the base TurbulenceModelConfiguration class GetTurbulenceSolvingProcess method."
        msg += " Please override it in the derrived class to return a KratosMultiphysics.Process."
        raise Exception(msg)

    def Initialize(self):
        msg = "Calling the base TurbulenceModelConfiguration class Initialize method."
        msg += " Please override it in the derrived class."
        raise Exception(msg)

    def Check(self):
        msg = "Calling the base TurbulenceModelConfiguration class Check method."
        msg += " Please override it in the derrived class."
        raise Exception(msg)

    def InitializeSolutionStep(self):
        msg = "Calling the base TurbulenceModelConfiguration class InitializeSolutionStep method."
        msg += " Please override it in the derrived class."
        raise Exception(msg)

    def FinalizeSolutionStep(self):
        msg = "Calling the base TurbulenceModelConfiguration class FinalizeSolutionStep method."
        msg += " Please override it in the derrived class."
        raise Exception(msg)

    def GetFluidVelocityPressureConditionName(self):
        msg = "Calling the base TurbulenceModelConfiguration class GetFluidVelocityPressureConditionName method."
        msg += " Please override it in the derrived class to return a condition name."
        raise Exception(msg)