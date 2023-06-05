# Importing the Kratos Library
import KratosMultiphysics as KM

# Import applications
import KratosMultiphysics.OptimizationApplication as KOA

# Import baseclass
from KratosMultiphysics.OptimizationApplication.filtering.helmholtz_solver_base import HelmholtzSolverBase

def CreateSolver(model: KM.Model, custom_settings: KM.Parameters):
    return HelmholtzScalarSolver(model, custom_settings)

class HelmholtzScalarSolver(HelmholtzSolverBase):
    def AddVariables(self):
        # Add variables required for the helmholtz filtering
        self.original_model_part.AddNodalSolutionStepVariable(KOA.HELMHOLTZ_SCALAR)
        self.helmholtz_model_part.AddNodalSolutionStepVariable(KOA.HELMHOLTZ_SCALAR)
        KM.Logger.PrintInfo("::[HelmholtzScalarSolver]:: Variables ADDED.")

    def AddDofs(self):
        KM.VariableUtils().AddDof(KOA.HELMHOLTZ_SCALAR, self.helmholtz_model_part)
        KM.Logger.PrintInfo("::[HelmholtzScalarSolver]:: DOFs ADDED.")

    def PrepareModelPart(self):
        #check elements types
        is_surface = False
        num_nodes = None
        for elem in self.original_model_part.Elements:
            geom = elem.GetGeometry()
            if geom.WorkingSpaceDimension() != geom.LocalSpaceDimension():
                is_surface = True
            num_nodes = len(elem.GetNodes())
            break

        if is_surface:
            element_name = f"HelmholtzSurfaceElement3D{num_nodes}N"
        else:
            element_name = f"HelmholtzSolidElement3D{num_nodes}N"

        KM.ConnectivityPreserveModeler().GenerateModelPart(self.original_model_part, self.helmholtz_model_part, element_name)
