# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.OptimizationApplication as KOA

# Import baseclass
from KratosMultiphysics.OptimizationApplication.filtering.helmholtz_solver_base import HelmholtzSolverBase

def CreateSolver(model: KM.Model, custom_settings: KM.Parameters):
    return HelmholtzShapeSolver(model, custom_settings)

class HelmholtzShapeSolver(HelmholtzSolverBase):
    def GetSolvingVariable(self) -> KM.Array1DVariable3:
        return KOA.HELMHOLTZ_VECTOR

    def _GetComputingModelPartName(self) -> str:
        return self.GetOriginModelPart().FullName().replace(".", "_") + "_helmholtz_shape"

    def SetFilterRadius(self, filter_radius: float) -> None:
        self.GetComputingModelPart().ProcessInfo.SetValue(KOA.HELMHOLTZ_RADIUS, filter_radius)
        KOA.ImplicitFilterUtils.SetBulkRadiusForShapeFiltering(self.GetComputingModelPart())

    def _FillComputingModelPart(self) -> None:
        if self.GetOriginModelPart().NumberOfElements() > 0:
            # here we have to replace the elements, while keeping the
            # geometries the same. Hence using the ConnectivityPreserveModeller
            container = self.GetOriginModelPart().Elements
            num_nodes = self._GetContainerTypeNumNodes(container)
            if self._IsSurfaceContainer(container):
                KM.ConnectivityPreserveModeler().GenerateModelPart(self.GetOriginModelPart(), self.GetComputingModelPart(), f"HelmholtzVectorSurfaceElement3D{num_nodes}N")
            else:
                raise RuntimeError(f"Helmholtz shape solver found volume elements in {self.GetOriginModelPart().FullName()}. It can work with only either surface elements or surface conditions.")
        elif self.GetOriginModelPart().NumberOfConditions() > 0:
            # we have only conditions in the origin model part. Since this solver is used to
            # do the mesh motion as well, then we need to work on the origin model part's root model part.
            # here also we can use the connectivity preserve modeller.
            element_num_nodes = self._GetContainerTypeNumNodes(self.GetOriginRootModelPart().Elements)
            KM.ConnectivityPreserveModeler().GenerateModelPart(self.GetOriginRootModelPart(), self.GetComputingModelPart(), f"HelmholtzSolidShapeElement3D{element_num_nodes}N", f"HelmholtzSurfaceShapeCondition3D{self._GetContainerTypeNumNodes(self.GetOriginRootModelPart().Conditions)}N")

            # now we have to add the material properties
            properties = KM.Parameters("""{
                "properties" : [
                    {
                        "model_part_name": \"""" + self.GetComputingModelPart().FullName() + """\",
                        "properties_id"  : """ + str(self.GetComputingModelPart().NumberOfProperties() + 1) + """,
                        "Material": {
                            "constitutive_law": {
                                "name": "HelmholtzJacobianStiffened3D"
                            },
                            "Variables": {
                                "POISSON_RATIO": 0.3
                            }
                        }
                    }
                ]
            }""")
            KM.ReadMaterialsUtility(self.model).ReadMaterials(properties)

            # now we do the tetrahedral mesh orientation check and parent element
            # assignment.
            if element_num_nodes == 4:
                tmoc = KM.TetrahedralMeshOrientationCheck
                flags = (tmoc.COMPUTE_NODAL_NORMALS).AsFalse() | (tmoc.COMPUTE_CONDITION_NORMALS).AsFalse() | tmoc.ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS
                KM.TetrahedralMeshOrientationCheck(self.GetComputingModelPart(), False, flags).Execute()
            else:
                raise RuntimeError("Helmholtz shape solver only supports tetrahedral elements.")
        else:
            raise RuntimeError(f"No elements or conditions found in {self.GetOriginModelPart()}.")
