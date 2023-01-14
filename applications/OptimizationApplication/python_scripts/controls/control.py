from abc import ABC
from abc import abstractmethod

import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.utilities.helper_utils import ContainerEnum

class Control(ABC):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        self.model = model
        self.parameters = parameters
        self.optimization_info = optimization_info
        self.__is_control_update_computed = False
        self.__control_update = Kratos.Vector()

    def Initialize(self):
        pass

    def InitializeSolutionStep(self):
        self.__is_control_update_computed = False

    def FinalizeSolutionStep(self):
        pass

    def Finalize(self):
        pass

    def SetControlUpdatesVector(self, control_update: Kratos.Vector):
        self.__is_control_update_computed = True
        self.__control_update = control_update

    def GetControlUpdateVector(self) -> Kratos.Vector:
        if not self.IsControlUpdateComputed():
            raise RuntimeError("Control update is not computed for current solution step.")
        return self.__control_update

    def IsControlUpdateComputed(self) -> bool:
        return self.__is_control_update_computed

    @abstractmethod
    def GetContainerType(self) -> ContainerEnum:
        """Returns the contaienr type on which the control is acted upon

        This method returns the container type, where the sensitivites/updates
        are carried on.

        Returns:
            ContainerEnum: Controls contaienr type
        """
        pass

    @abstractmethod
    def GetControlUpdateVariable(self) -> any:
        """Returns the control update variable

        This method returns the update value holding variable for the control.
        Eg: For DENSITY control variable, one can use DENSITY_UPDATE

        Returns:
            any: Control update variable
        """
        pass

    @abstractmethod
    def GetControlSensitivityVariable(self) -> any:
        """Returns the sensitivity computation variable

        This method returns the variable on which the computed sensitivities
        for a respective response is stored.
        Eg: For DENSITY control variable, DENSITY_SENSITIVITY can be used.

        Returns:
            any: Returns the sensitivity storage variable
        """
        pass

    @abstractmethod
    def GetModelPart(self) -> Kratos.ModelPart:
        """Returns the model part on which this controller is used.

        Returns:
            Kratos.ModelPart: Controlled model part
        """
        pass

