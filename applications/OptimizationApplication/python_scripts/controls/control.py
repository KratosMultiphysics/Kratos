from abc import ABC, abstractmethod
import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes

class Control(ABC):
    """Base abstract control class.

    The base abstract control class has the following responsibilities.
        1. Initializes and finalizes the control.
        2. Maps physical space gradients to control space gradients, if needed. Otherwise, it passes.
        2. Updating the controlled parts of the model part with new given design.
        3. Retrieve control field.
        4. Retrieve physical space kratos variables (If more than one physical space kratos variables are controlled by the given control field.)

    This control should only work on one model part and one kratos control variable. Hence, if multiple model parts required then,
    a single model part should be created using Kratos.ModelPartOperationUtilities.

    """
    def __init__(self, control_name: str) -> None:
        self.__name = control_name

    def GetName(self) -> str:
        return self.__name

    @abstractmethod
    def Initialize(self) -> None:
        """Initializes the control.

        This method initializes the control. This is only called once in the whole optimization process.

        """
        pass

    @abstractmethod
    def Check(self) -> None:
        """Checks the control.

        This method checks the control. This is only called once in the whole optimization process.

        """
        pass

    @abstractmethod
    def Finalize(self) -> None:
        """Finalizes the control.

        This method finalizes the control. This is only called once in the whole optimization process.

        """
        pass

    @abstractmethod
    def GetPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        """Returns list of physical control variables controlled by this control.

        This method returns list of physical variable which are controlled by the control. In most of the cases, there is only one
        control variable per control. But in some rare cases, there may be the requirement to control two physical variables. Anyhow,
        the control space should only have one variable per control in the control space, but they can have multiple physical variables
        in their design space.

        Returns:
            list[SupportedSensitivityFieldVariableTypes]: List of physical control variables.
        """
        pass

    @abstractmethod
    def GetEmptyField(self) -> Kratos.TensorAdaptors.DoubleTensorAdaptor:
        """Returns a new empty data field holder with correct dimensionality information.

        This returns a new empty data field holder to give information about on which model part's container
        this model part is acting on. It does not read any data from respective model part's container.

        Returns:
            Kratos.TensorAdaptors.DoubleTensorAdaptor: Returns a new empty @ref Kratos.TensorAdaptor corresponding to control's model part's respective container.
        """
        pass

    @abstractmethod
    def GetControlField(self) -> Kratos.TensorAdaptors.DoubleTensorAdaptor:
        """Returns the current control field of the control.

        This method returns the control field of the current design.

        Returns:
            Kratos.TensorAdaptors.DoubleTensorAdaptor: Current designs control field.
        """
        pass

    @abstractmethod
    def MapGradient(self, physical_gradient_variable_tensor_adaptor_map: 'dict[SupportedSensitivityFieldVariableTypes, Kratos.TensorAdaptors.DoubleTensorAdaptor]') -> Kratos.TensorAdaptors.DoubleTensorAdaptor:
        """Maps physical space gradients to the control space.

        This method is used to map the given physical space gradients to the control space. The input should be as in the following example:
            physical_gradient_variable_tensor_adaptor_map = {
                Kratos.YOUNG_MODULUS: Kratos.DoubleTensorAdaptor,
                Kratos.DENSITY      : Kratos.DoubleTensorAdaptor,
                Kratos.SHAPE        : Kratos.DoubleTensorAdaptor
            }

        All the gradients w.r.t. @see GetPhysicalKratosVariables() variables will be given in @ref physical_gradient_variable_tensor_adaptor_map.
        If the response does not depend on some of them or all, then @ref Kratos.TensorAdaptor with correctly sized zero tensors will be passed.

        Args:
            physical_gradient_variable_tensor_adaptor_map (dict[SupportedSensitivityFieldVariableTypes, Kratos.TensorAdaptors.DoubleTensorAdaptor]): Map of physical space variable and @ref Kratos::TensorAdaptor with sensitivities.

        Returns:
            Kratos.TensorAdaptors.DoubleTensorAdaptor: Gradients mapped in to control space.
        """
        pass

    @abstractmethod
    def Update(self, control_field: Kratos.TensorAdaptors.DoubleTensorAdaptor) -> bool:
        """Modifies the current control with the given control field.

        Args:
            control_field (Kratos.TensorAdaptors.DoubleTensorAdaptor): The control field in control space.

        Returns:
            bool: True if the control field was applied to obtain a new design, otherwise False
        """
        pass


