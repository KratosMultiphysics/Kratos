# Importing the Kratos Library
import KratosMultiphysics as KM

# Import applications
import KratosMultiphysics.OptimizationApplication as KOA

# Import baseclass
from KratosMultiphysics.OptimizationApplication.filtering.helmholtz_solver_base import HelmholtzSolverBase

def CreateSolver(model: KM.Model, custom_settings: KM.Parameters):
    return HelmholtzScalarSolver(model, custom_settings)

class HelmholtzScalarSolver(HelmholtzSolverBase):
    def GetSolvingVariable(self) -> KM.DoubleVariable:
        return KOA.HELMHOLTZ_SCALAR

    def _GetComputingModelPartName(self) -> str:
        return self.GetOriginModelPart().FullName().replace(".", "_") + "_helmholtz_scalar"

    def SetFilterRadius(self, filter_radius: float) -> None:
        self.GetComputingModelPart().ProcessInfo.SetValue(KOA.HELMHOLTZ_RADIUS, filter_radius)

    def _FillComputingModelPart(self) -> None:
        if self.GetOriginModelPart().NumberOfElements() > 0:
            # here we have to replace the elements, while keeping the
            # geometries the same. Hence using the ConnectivityPreserveModeller
            container = self.GetOriginModelPart().Elements
            num_nodes = self._GetContainerTypeNumNodes(container)
            if self._IsSurfaceContainer(container):
                KM.ConnectivityPreserveModeler().GenerateModelPart(self.GetOriginModelPart(), self.GetComputingModelPart(), f"HelmholtzSurfaceElement3D{num_nodes}N")
            else:
                KM.ConnectivityPreserveModeler().GenerateModelPart(self.GetOriginModelPart(), self.GetComputingModelPart(), f"HelmholtzSolidElement3D{num_nodes}N")
        elif self.GetOriginModelPart().NumberOfConditions() > 0:
            # here we have conditions in the origin model part. Now we have to create elements using
            # the geometries of the conditions. There cannot be volume conditions in Kratos, therefore,
            # there can be only surface conditions, hence only required to create surface elements.
            KOA.OptAppModelPartUtils.GenerateModelPart(self.GetOriginModelPart().Conditions, self.GetComputingModelPart(), f"HelmholtzSurfaceElement3D{self._GetContainerTypeNumNodes(self.GetOriginModelPart().Conditions)}N")
        else:
            raise RuntimeError(f"No elements or conditions found in {self.GetOriginModelPart()}.")