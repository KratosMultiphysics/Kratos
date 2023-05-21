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

        #check original model part has both conditions and elements
        if self.original_model_part.NumberOfConditions()<1:
            raise Exception('::[HelmholtzBulkSurfaceSolver]:: given model part must have surface conditions')

        num_elems_nodes = None
        is_surface = False
        for elem in self.original_model_part.Elements:
            geom = elem.GetGeometry()
            if geom.WorkingSpaceDimension() != geom.LocalSpaceDimension():
                is_surface = True
            num_elems_nodes = len(elem.GetNodes())
            break

        num_conds_nodes = None
        for cond in self.original_model_part.Conditions:
            geom = cond.GetGeometry()
            num_conds_nodes = len(cond.GetNodes())
            break

        if num_elems_nodes != 4:
            raise Exception('::[HelmholtzBulkSurfaceSolver]:: given model part must have only tetrahedral elemenst')
        if num_conds_nodes != 3:
            raise Exception('::[HelmholtzBulkSurfaceSolver]:: given model part must have only triangular conditions')
        if is_surface:
            raise Exception('::[HelmholtzBulkSurfaceSolver]:: given model part must have only solid elemenst')

        KM.ConnectivityPreserveModeler().GenerateModelPart(
            self.original_model_part, self.helmholtz_model_part, "HelmholtzSolidShapeElement3D4N","HelmholtzSurfaceShapeCondition3D3N")

        tmoc = KM.TetrahedralMeshOrientationCheck
        flags = (tmoc.COMPUTE_NODAL_NORMALS).AsFalse() | (tmoc.COMPUTE_CONDITION_NORMALS).AsFalse() | tmoc.ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS
        KM.TetrahedralMeshOrientationCheck(self.helmholtz_model_part, False, flags).Execute()




