# Importing the Kratos Library
import KratosMultiphysics as KM

# Import applications
import KratosMultiphysics.OptimizationApplication as KOA

# Import baseclass
from KratosMultiphysics.OptimizationApplication.filtering.helmholtz_solver_base import HelmholtzSolverBase


def CreateSolver(model, custom_settings):
    return HelmholtzScalarSolver(model, custom_settings)


class HelmholtzScalarSolver(HelmholtzSolverBase):
    def __init__(self, model, custom_settings):
        super().__init__(model, custom_settings)
        self.filter_radius = self.settings["filter_radius"].GetDouble()
        KM.Logger.PrintInfo("::[HelmholtzScalarSolver]:: Construction finished")

    def AddVariables(self):
        # Add variables required for the helmholtz filtering
        self.original_model_part.AddNodalSolutionStepVariable(KOA.HELMHOLTZ_SCALAR)
        self.helmholtz_model_part.AddNodalSolutionStepVariable(KOA.HELMHOLTZ_SCALAR)
        KM.Logger.PrintInfo("::[HelmholtzScalarSolver]:: Variables ADDED.")

    def AddDofs(self):
        KM.VariableUtils().AddDof(KOA.HELMHOLTZ_SCALAR, self.helmholtz_model_part)
        KM.Logger.PrintInfo("::[HelmholtzScalarSolver]:: DOFs ADDED.")

    def PrepareModelPart(self):
        KM.ConnectivityPreserveModeler().GenerateModelPart(
                self.original_model_part, self.helmholtz_model_part, "HelmholtzSolidElement3D8N")

