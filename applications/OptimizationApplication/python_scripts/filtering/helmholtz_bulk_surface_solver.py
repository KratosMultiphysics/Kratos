# Importing the Kratos Library
import KratosMultiphysics as KM

# Import applications
import KratosMultiphysics.OptimizationApplication as KOA

# Import baseclass
from KratosMultiphysics.OptimizationApplication.filtering.helmholtz_solver_base import HelmholtzSolverBase


def CreateSolver(model, custom_settings):
    return HelmholtzBulkSurfaceSolver(model, custom_settings)


class HelmholtzBulkSurfaceSolver(HelmholtzSolverBase):
    def __init__(self, model, custom_settings):
        super().__init__(model, custom_settings)
        KM.Logger.PrintInfo("::[HelmholtzBulkSurfaceSolver]:: Construction finished")

    def AddVariables(self):
        # Add variables required for the helmholtz filtering
        self.original_model_part.AddNodalSolutionStepVariable(KOA.HELMHOLTZ_VECTOR)
        self.helmholtz_model_part.AddNodalSolutionStepVariable(KOA.HELMHOLTZ_VECTOR)
        KM.Logger.PrintInfo("::[HelmholtzBulkSurfaceSolver]:: Variables ADDED.")

    def AddDofs(self):
        KM.VariableUtils().AddDof(KOA.HELMHOLTZ_VECTOR_X, self.helmholtz_model_part)
        KM.VariableUtils().AddDof(KOA.HELMHOLTZ_VECTOR_Y, self.helmholtz_model_part)
        KM.VariableUtils().AddDof(KOA.HELMHOLTZ_VECTOR_Z, self.helmholtz_model_part)
        KM.Logger.PrintInfo("::[HelmholtzBulkSurfaceSolver]:: DOFs ADDED.")

    def PrepareModelPart(self):
        pass

        # #check elements types
        # is_surface = False
        # num_nodes = None
        # for elem in self.original_model_part.Elements:
        #     geom = elem.GetGeometry()
        #     if geom.WorkingSpaceDimension() != geom.LocalSpaceDimension():
        #         is_surface = True
        #     num_nodes = len(elem.GetNodes())
        #     break

        # if is_surface:
        #     if num_nodes == 3:
        #         KM.ConnectivityPreserveModeler().GenerateModelPart(
        #                 self.original_model_part, self.helmholtz_model_part, "HelmholtzSurfaceElement3D3N")
        #     elif num_nodes == 4:
        #         KM.ConnectivityPreserveModeler().GenerateModelPart(
        #                 self.original_model_part, self.helmholtz_model_part, "HelmholtzSurfaceElement3D4N")
        # else:
        #     if num_nodes == 4:
        #         KM.ConnectivityPreserveModeler().GenerateModelPart(
        #                 self.original_model_part, self.helmholtz_model_part, "HelmholtzSolidElement3D4N")
        #     elif num_nodes == 8:
        #         KM.ConnectivityPreserveModeler().GenerateModelPart(
        #                 self.original_model_part, self.helmholtz_model_part, "HelmholtzSolidElement3D8N")



