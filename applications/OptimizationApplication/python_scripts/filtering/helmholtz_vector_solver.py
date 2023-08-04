# Importing the Kratos Library
import KratosMultiphysics as KM

# Import applications
import KratosMultiphysics.OptimizationApplication as KOA

# Import baseclass
from KratosMultiphysics.OptimizationApplication.filtering.helmholtz_solver_base import HelmholtzSolverBase

def CreateSolver(model, custom_settings):
    return HelmholtzVectorSolver(model, custom_settings)

class HelmholtzVectorSolver(HelmholtzSolverBase):
    def __init__(self, model: KM.Model, custom_settings: KM.Parameters) -> None:
        super().__init__(model, custom_settings)
        if self.settings["filter_type"].GetString() == "bulk_surface_shape":
            self.bulk_surface_shape_filter = True
        else:
            self.bulk_surface_shape_filter = False
        KM.Logger.PrintInfo("::[HelmholtzVectorSolver]:: Construction finished")

    def AddVariables(self) -> None:
        # Add variables required for the helmholtz filtering
        self.GetOriginRootModelPart().AddNodalSolutionStepVariable(KOA.HELMHOLTZ_VECTOR)
        self.GetComputingModelPart().AddNodalSolutionStepVariable(KOA.HELMHOLTZ_VECTOR)
        KM.Logger.PrintInfo("::[HelmholtzVectorSolver]:: Variables ADDED.")

    def AddDofs(self) -> None:
        KM.VariableUtils().AddDof(KOA.HELMHOLTZ_VECTOR_X, self.GetComputingModelPart())
        KM.VariableUtils().AddDof(KOA.HELMHOLTZ_VECTOR_Y, self.GetComputingModelPart())
        KM.VariableUtils().AddDof(KOA.HELMHOLTZ_VECTOR_Z, self.GetComputingModelPart())
        KM.Logger.PrintInfo("::[HelmholtzVectorSolver]:: DOFs ADDED.")

    def PrepareModelPart(self) -> None:
        num_elems_nodes = None
        is_surface = False
        for elem in self.GetOriginModelPart().Elements:
            geom = elem.GetGeometry()
            if geom.WorkingSpaceDimension() != geom.LocalSpaceDimension():
                is_surface = True
            num_elems_nodes = len(elem.GetNodes())
            break

        num_conds_nodes = None
        for cond in self.GetOriginModelPart().Conditions:
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
                self.GetOriginModelPart(), self.GetComputingModelPart(), "HelmholtzSolidShapeElement3D4N","HelmholtzSurfaceShapeCondition3D3N")

            material_properties = self.settings["material_properties"]
            defaults = KM.Parameters("""{
                "properties_id": 1,
                "Material": {
                    "constitutive_law": {
                        "name": "HelmholtzJacobianStiffened3D"
                    },
                    "Variables": {
                        "POISSON_RATIO": 0.3
                    }
                }
            }""")
            material_properties.RecursivelyAddMissingParameters(defaults)
            self._AssignProperties(material_properties)

            tmoc = KM.TetrahedralMeshOrientationCheck
            flags = (tmoc.COMPUTE_NODAL_NORMALS).AsFalse() | (tmoc.COMPUTE_CONDITION_NORMALS).AsFalse() | tmoc.ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS
            KM.TetrahedralMeshOrientationCheck(self.GetComputingModelPart(), False, flags).Execute()
        else:
            if is_surface:
                element_name = f"HelmholtzVectorSurfaceElement3D{num_elems_nodes}N"
            else:
                element_name = f"HelmholtzVectorSolidElement3D{num_elems_nodes}N"

            KM.ConnectivityPreserveModeler().GenerateModelPart(self.GetOriginModelPart(), self.GetComputingModelPart(), element_name)
