import typing

# Importing the Kratos Library
import KratosMultiphysics as KM

# Import applications
import KratosMultiphysics.OptimizationApplication as KOA

# Import baseclass
from KratosMultiphysics.OptimizationApplication.filtering.helmholtz_solver_base import HelmholtzSolverBase

def CreateSolver(model: KM.Model, custom_settings: KM.Parameters):
    return HelmholtzScalarSolver(model, custom_settings)

class HelmholtzScalarSolver(HelmholtzSolverBase):
    def __init__(self, model: KM.Model, custom_settings: KM.Parameters):
        super().__init__(model, custom_settings)
        self.__containers: 'typing.Union[KM.ConditionsArray, KM.ElementsArray]' = []
        self.__element_name = ""
        KM.Logger.PrintInfo("::[HelmholtzScalarSolver]:: Construction finished")

    def AddVariables(self) -> None:
        # Add variables required for the helmholtz filtering
        self.GetOriginRootModelPart().AddNodalSolutionStepVariable(KOA.HELMHOLTZ_SCALAR)
        KM.Logger.PrintInfo("::[HelmholtzScalarSolver]:: Variables ADDED.")

    def AddDofs(self) -> None:
        KM.VariableUtils().AddDof(KOA.HELMHOLTZ_SCALAR, self.GetOriginRootModelPart())
        KM.Logger.PrintInfo("::[HelmholtzScalarSolver]:: DOFs ADDED.")

    def InitializeModelPartConfiguration(self) -> None:
        number_of_conditions = self.GetOriginModelPart().NumberOfConditions()
        number_of_elements = self.GetOriginModelPart().NumberOfElements()

        if number_of_conditions > 0 and number_of_elements > 0:
            KM.Logger.PrintWarning("::[HelmholtzScalarSolver]::", f"Filter model part {self.GetOriginModelPart().FullName()} has both elements and conditions. Giving precedence to elements.")

        if number_of_elements > 0:
            self.__containers = [self.GetOriginModelPart().Elements]
        elif number_of_conditions > 0:
           self.__containers = [self.GetOriginModelPart().Conditions]
        else:
           raise RuntimeError(f"No elements or conditions found in the filtering model part \"{self.GetOriginModelPart().FullName()}\"")

        num_nodes = self._GetContainerTypeNumNodes(self.__containers[0])

        if self._IsSurfaceContainer(self.__containers[0]):
            self.__element_name = f"HelmholtzSurfaceElement3D{num_nodes}N"
        else:
            self.__element_name = f"HelmholtzSolidElement3D{num_nodes}N"

    def GetElementName(self) -> str:
        return self.__element_name

    def GetConditionName(self) -> str:
        return ""

    def GetContainers(self) -> 'list[typing.Union[KM.ConditionsArray | KM.ElementsArray]]':
        return self.__containers

    def GetMaterialProperties(self) -> KM.Parameters:
        return KM.Parameters("""{}""")