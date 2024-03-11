# Importing the Kratos Library
import KratosMultiphysics as KM

# Import applications
import KratosMultiphysics.OptimizationApplication as KOA

# Import baseclass
from KratosMultiphysics.OptimizationApplication.filtering.helmholtz_solver_base import HelmholtzSolverBase

def CreateSolver(model: KM.Model, custom_settings: KM.Parameters):
    return HelmholtzScalarSolver(model, custom_settings)

class HelmholtzScalarSolver(HelmholtzSolverBase):

    def AddVariables(self) -> None:
        # Add variables required for the helmholtz filtering
        self.GetOriginRootModelPart().AddNodalSolutionStepVariable(KOA.HELMHOLTZ_SCALAR)
        KM.Logger.PrintInfo("::[HelmholtzScalarSolver]:: Variables ADDED.")

    def AddDofs(self) -> None:
        KM.VariableUtils().AddDof(KOA.HELMHOLTZ_SCALAR, self.GetOriginRootModelPart())
        KM.Logger.PrintInfo("::[HelmholtzScalarSolver]:: DOFs ADDED.")

    def PrepareModelPart(self) -> None:

        if len(self.GetOriginModelPart().Conditions)>0 and len(self.GetOriginModelPart().Elements)>0:
            KM.Logger.PrintWarning("::[HelmholtzScalarSolver]:: filter model part ", self.GetOriginModelPart().Name, " has both elements and conditions. Giving precedence to elements ")

        if len(self.GetOriginModelPart().Elements)>0:
           filter_container = self.GetOriginModelPart().Elements
        elif len(self.GetOriginModelPart().Conditions)>0:
           filter_container = self.GetOriginModelPart().Conditions

        is_surface_filter = self._IsSurfaceContainer(filter_container)
        num_nodes = self._GetContainerTypeNumNodes(filter_container)

        if is_surface_filter:
            element_name = f"HelmholtzSurfaceElement3D{num_nodes}N"
        else:
            element_name = f"HelmholtzSolidElement3D{num_nodes}N"

        KM.ConnectivityPreserveModeler().GenerateModelPart(self.GetOriginModelPart(), self.helmholtz_model_part, element_name)