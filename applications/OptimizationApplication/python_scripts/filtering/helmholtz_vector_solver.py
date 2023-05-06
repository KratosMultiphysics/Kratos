# Importing the Kratos Library
import KratosMultiphysics as KM

# Import applications
import KratosMultiphysics.OptimizationApplication as KOA

# Import baseclass
from KratosMultiphysics.OptimizationApplication.helmholtz_solver_base import HelmholtzSolverBase


def CreateSolver(model, custom_settings):
    return HelmholtzVectorSolver(model, custom_settings)


class HelmholtzVectorSolver(HelmholtzSolverBase):
    def __init__(self, model, custom_settings):
        super().__init__(model, custom_settings)
        KM.Logger.PrintInfo("::[HelmholtzVectorSolver]:: Construction finished")

    def AddVariables(self):
        # Add variables required for the helmholtz filtering
        self.helmholtz_model_part.AddNodalSolutionStepVariable(KOA.HELMHOLTZ_VECTOR)
        KM.Logger.PrintInfo("::[HelmholtzVectorSolver]:: Variables ADDED.")

    def AddDofs(self):
        KM.VariableUtils().AddDof(KOA.HELMHOLTZ_VECTOR_X, self.helmholtz_model_part)
        KM.VariableUtils().AddDof(KOA.HELMHOLTZ_VECTOR_Y, self.helmholtz_model_part)
        KM.VariableUtils().AddDof(KOA.HELMHOLTZ_VECTOR_Z, self.helmholtz_model_part)
        KM.Logger.PrintInfo("::[HelmholtzVectorSolver]:: DOFs ADDED.")