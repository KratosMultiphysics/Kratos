# Importing the Kratos Library
import KratosMultiphysics as KM

# Import applications
import KratosMultiphysics.OptimizationApplication as KOA

# Import baseclass
from KratosMultiphysics.OptimizationApplication.filtering.helmholtz_solver_base import HelmholtzSolverBase


def CreateSolver(model, custom_settings):
    return HelmholtzVectorSolver(model, custom_settings)


class HelmholtzVectorSolver(HelmholtzSolverBase):
    def __init__(self, model, custom_settings):
        super().__init__(model, custom_settings)
        if self.settings["filter_type"].GetString() == "bulk_surface_shape":
            self.bulk_surface_shape_filter = True
        else:
            self.bulk_surface_shape_filter = False
        KM.Logger.PrintInfo("::[HelmholtzVectorSolver]:: Construction finished")

    def AddVariables(self):
        # Add variables required for the helmholtz filtering
        self.original_model_part.AddNodalSolutionStepVariable(KOA.HELMHOLTZ_VECTOR)
        self.helmholtz_model_part.AddNodalSolutionStepVariable(KOA.HELMHOLTZ_VECTOR)
        KM.Logger.PrintInfo("::[HelmholtzVectorSolver]:: Variables ADDED.")

    def AddDofs(self):
        KM.VariableUtils().AddDof(KOA.HELMHOLTZ_VECTOR_X, self.helmholtz_model_part)
        KM.VariableUtils().AddDof(KOA.HELMHOLTZ_VECTOR_Y, self.helmholtz_model_part)
        KM.VariableUtils().AddDof(KOA.HELMHOLTZ_VECTOR_Z, self.helmholtz_model_part)
        KM.Logger.PrintInfo("::[HelmholtzVectorSolver]:: DOFs ADDED.")

    def PrepareModelPart(self):

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

        if self.bulk_surface_shape_filter:

            if num_elems_nodes != 4:
                raise Exception('::[HelmholtzVectorSolver]:: given model part must have only tetrahedral elemenst')
            if num_conds_nodes != 3:
                raise Exception('::[HelmholtzVectorSolver]:: given model part must have only triangular conditions')
            if is_surface:
                raise Exception('::[HelmholtzVectorSolver]:: given model part must have only solid elemenst')

            KM.ConnectivityPreserveModeler().GenerateModelPart(
                self.original_model_part, self.helmholtz_model_part, "HelmholtzSolidShapeElement3D4N","HelmholtzSurfaceShapeCondition3D3N")

            # Create property for the solid shape element
            properties = self.helmholtz_model_part.GetProperties(0)
            cl = KOA.HelmholtzJacobianStiffened3DLaw()
            properties.SetValue(KM.CONSTITUTIVE_LAW, cl)

            tmoc = KM.TetrahedralMeshOrientationCheck
            flags = (tmoc.COMPUTE_NODAL_NORMALS).AsFalse() | (tmoc.COMPUTE_CONDITION_NORMALS).AsFalse() | tmoc.ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS
            KM.TetrahedralMeshOrientationCheck(self.helmholtz_model_part, False, flags).Execute()
        else:
            if is_surface:
                if num_elems_nodes == 3:
                    KM.ConnectivityPreserveModeler().GenerateModelPart(
                            self.original_model_part, self.helmholtz_model_part, "HelmholtzVectorSurfaceElement3D3N")
                elif num_elems_nodes == 4:
                    KM.ConnectivityPreserveModeler().GenerateModelPart(
                            self.original_model_part, self.helmholtz_model_part, "HelmholtzVectorSurfaceElement3D4N")
            else:
                if num_elems_nodes == 4:
                    KM.ConnectivityPreserveModeler().GenerateModelPart(
                            self.original_model_part, self.helmholtz_model_part, "HelmholtzVectorSolidElement3D4N")
                elif num_elems_nodes == 8:
                    KM.ConnectivityPreserveModeler().GenerateModelPart(
                            self.original_model_part, self.helmholtz_model_part, "HelmholtzVectorSolidElement3D8N")




