import typing

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
        self.__containers: 'list[typing.Union[KM.ConditionsArray, KM.ElementsArray]]' = []
        self.__element_name = ""
        self.__condition_name = ""
        self.__materials = KM.Parameters("""{}""")
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

    def __InitializeBulkSurfaceShapeModelPartConfiguration(self) -> None:
        self.__containers = [self.GetOriginModelPart().Conditions, self.GetOriginRootModelPart().Elements]
        self.__condition_name = f"HelmholtzSurfaceShapeCondition3D{self._GetContainerTypeNumNodes(self.__containers[0])}N"
        self.__element_name = f"HelmholtzSolidShapeElement3D{self._GetContainerTypeNumNodes(self.__containers[1])}N"
        self.__materials = KM.Parameters("""{
                "constitutive_law": {
                    "name": "HelmholtzJacobianStiffened3D"
                },
                "Variables": {
                    "POISSON_RATIO": 0.3
                }
        }""")

    def InitializeModelPartConfiguration(self) -> None:
        filter_type = self.GetFilterType()
        number_of_conditions = self.GetOriginModelPart().NumberOfConditions()
        number_of_elements = self.GetOriginModelPart().NumberOfElements()
        self.__materials = KM.Parameters("""{}""")

        if filter_type == "shape":
            # this is the filter type which is used to switch automatically between
            # bulk_surface_shape and general_vector filters based on the type
            # of entities present in the filtering model part.

            if number_of_elements > 0:
                # if there are elements in the filtering model part, then those
                # elements should be surface elements (such as shell), otherwise
                # we cannot do shape filtering. Shape filtering must only be done
                # on surfaces. If someone wants to do filtering on conditions where
                # there are elements in the model part, then
                # the proper sub-model part should be chosen for filtering model part.
                self.__containers = [self.GetOriginModelPart().Elements]
                self.__element_name = f"HelmholtzVectorSurfaceElement3D{self._GetContainerTypeNumNodes(self.__containers[0])}N"
            elif number_of_conditions > 0:
                # now this model part only has conditions. So, we check whether they are
                # surface conditions and add both conditions and root model part elements
                self.__InitializeBulkSurfaceShapeModelPartConfiguration()
                raise RuntimeError(1)
            else:
                raise RuntimeError(f"The filtering model part \"{self.GetOriginModelPart().FullName()} does not have elements or conditions.")
        elif filter_type == "bulk_surface_shape":
            self.__InitializeBulkSurfaceShapeModelPartConfiguration()
        elif filter_type == "general_vector":
            # this is to use helmholtz for general vector filtering. Again here we have to
            # check whether we have elements and conditions. This needs to support
            # surface elements, bulk elements, surface conditions.
            if number_of_elements > 0:
                # filter surface has elements.
                self.__containers = [self.GetOriginModelPart().Elements]
                if self._IsSurfaceContainer(self.__containers[0]):
                    # filter is used on a model part with surface elements
                    self.__element_name = f"HelmholtzVectorSurfaceElement3D{self._GetContainerTypeNumNodes(self.__containers[0])}N"
                else:
                    # filter is used on a model part with solid elements
                    self.__element_name = f"HelmholtzVectorSolidElement3D{self._GetContainerTypeNumNodes(self.__containers[0])}N"
            elif number_of_conditions > 0:
                self.__containers = [self.GetOriginModelPart().Conditions]
                if self._IsSurfaceContainer(self.__containers[0]):
                    # filter is used on a model part with surface conditions
                    self.__element_name = f"HelmholtzVectorSurfaceElement3D{self._GetContainerTypeNumNodes(self.__containers[0])}N"
                else:
                    raise RuntimeError("Unsupported model part.")
            else:
                raise RuntimeError(f"The filtering model part \"{self.GetOriginModelPart().FullName()} does not have elements or conditions.")

    def GetElementName(self) -> str:
        return self.__element_name

    def GetConditionName(self) -> str:
        return self.__condition_name

    def GetContainers(self) -> 'list[typing.Union[KM.ConditionsArray | KM.ElementsArray]]':
        return self.__containers

    def GetMaterialProperties(self) -> KM.Parameters:
        return self.__materials
