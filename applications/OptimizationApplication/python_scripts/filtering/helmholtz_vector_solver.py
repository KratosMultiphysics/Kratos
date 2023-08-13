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
        KM.Logger.PrintInfo("::[HelmholtzVectorSolver]:: Variables ADDED.")

    def AddDofs(self) -> None:
        KM.VariableUtils().AddDof(KOA.HELMHOLTZ_VECTOR_X, self.GetOriginRootModelPart())
        KM.VariableUtils().AddDof(KOA.HELMHOLTZ_VECTOR_Y, self.GetOriginRootModelPart())
        KM.VariableUtils().AddDof(KOA.HELMHOLTZ_VECTOR_Z, self.GetOriginRootModelPart())
        KM.Logger.PrintInfo("::[HelmholtzVectorSolver]:: DOFs ADDED.")

    def PrepareModelPart(self) -> None:

        if self.bulk_surface_shape_filter:

            num_root_elems_nodes = self._GetContainerTypeNumNodes(self.GetOriginModelPart().GetRootModelPart().Elements)
            is_root_surface = self._IsSurfaceContainer(self.GetOriginModelPart().GetRootModelPart().Elements)
            num_cond_nodes = self._GetContainerTypeNumNodes(self.GetOriginModelPart().Conditions)
            is_surface_condition = self._IsSurfaceContainer(self.GetOriginModelPart().Conditions)

            if num_root_elems_nodes != 4:
                raise Exception('::[HelmholtzVectorSolver]:: given model part must have only tetrahedral elemenst')
            if num_cond_nodes != 3:
                raise Exception('::[HelmholtzVectorSolver]:: given model part must have only triangular conditions')
            if not is_surface_condition or is_root_surface:
                raise Exception('::[HelmholtzVectorSolver]:: bulk surface should have solid tetrahedral elemenst and triangular conditions')

            # add nodes
            for node in self.GetOriginModelPart().GetRootModelPart().Nodes:
                self.helmholtz_model_part.AddNode(node)

            # add elems
            elem_properties = self.helmholtz_model_part.CreateNewProperties(self.helmholtz_model_part.NumberOfProperties()+1)
            elem_index = len(self.helmholtz_model_part.Elements) + 1
            for elem in self.GetOriginModelPart().GetRootModelPart().Elements:
                element_nodes_ids = []
                for node in elem.GetNodes():
                    element_nodes_ids.append(node.Id)
                self.helmholtz_model_part.CreateNewElement("HelmholtzSolidShapeElement3D4N", elem_index, element_nodes_ids, elem_properties)
                elem_index += 1

            # add conds
            cond_properties = self.helmholtz_model_part.CreateNewProperties(self.helmholtz_model_part.NumberOfProperties()+1)
            cond_index = len(self.helmholtz_model_part.GetRootModelPart().Conditions) + 1
            for cond in self.GetOriginModelPart().Conditions:
                cond_nodes_ids = []
                for node in cond.GetNodes():
                    cond_nodes_ids.append(node.Id)
                self.helmholtz_model_part.CreateNewCondition("HelmholtzSurfaceShapeCondition3D3N", cond_index, cond_nodes_ids, cond_properties)
                cond_index += 1

            material_properties = self.settings["material_properties"]
            defaults = KM.Parameters("""{
                "properties_id": 10000000000000000,
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
            KM.TetrahedralMeshOrientationCheck(self.helmholtz_model_part, False, flags).Execute()

        else:

            if len(self.GetOriginModelPart().Conditions)>0 and len(self.GetOriginModelPart().Elements)>0:
                KM.Logger.PrintWarning("::[HelmholtzVectorSolver]:: filter model part ", self.GetOriginModelPart().Name, " has both elements and conditions. Giving precedence to conditions ")

            if len(self.GetOriginModelPart().Conditions)>0:
               filter_container = self.GetOriginModelPart().Conditions
            elif len(self.GetOriginModelPart().Elements)>0:
               filter_container = self.GetOriginModelPart().Elements

            is_surface_filter = self._IsSurfaceContainer(filter_container)
            num_nodes = self._GetContainerTypeNumNodes(filter_container)

            if is_surface_filter:
                element_name = f"HelmholtzVectorSurfaceElement3D{num_nodes}N"
            else:
                element_name = f"HelmholtzVectorSolidElement3D{num_nodes}N"

            filter_properties = self.helmholtz_model_part.GetRootModelPart().CreateNewProperties(self.helmholtz_model_part.GetRootModelPart().NumberOfProperties()+1)
            for node in self.GetOriginModelPart().Nodes:
                self.helmholtz_model_part.AddNode(node)

            elem_index = len(self.helmholtz_model_part.GetRootModelPart().Elements) + 1
            for cond in filter_container:
                element_nodes_ids = []
                for node in cond.GetNodes():
                    element_nodes_ids.append(node.Id)
                self.helmholtz_model_part.CreateNewElement(element_name, elem_index, element_nodes_ids, filter_properties)
                elem_index += 1