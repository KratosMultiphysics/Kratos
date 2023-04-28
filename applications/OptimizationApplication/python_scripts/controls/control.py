from abc import ABC, abstractmethod
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes

class Control(ABC):
    """Base abstract control class.

    The base abstract control class has the following responsibilities.
        1. Initalizes and finalizes the control.
        2. Mapping physical space gradients to control space gradients.
        2. Updating the controlled parts of the model part with new design given in control space.
        3. Retrieve control field.
        4. Retrieve physical space kratos variables (If more than one physical space kratos variables are controlled by the given control field.)

    This control should only work on one model part and one kratos control variable. Hence, if multiple model parts required then,
    a single model part should be created using Kratos.ModelPartBinaryOperators.

    """
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

        If the case is that this control has one physical variable, then it should be same as the control variable.

        Returns:
            list[SupportedSensitivityFieldVariableTypes]: List of physical control variables.
        """
        pass

    @abstractmethod
    def GetEmptyControlField(self) -> ContainerExpressionTypes:
        """Returns a new empty control field data holder.

        This returns a new empty control field data holder to give information about on which model part's container
        this model part is acting on. This has O(1) complexity, hence has the least cost because it does not read
        any data from respective model part's container.

        Returns:
            ContainerExpressionTypes: Returns a new empty ContainerExpression corresponding to control's model part's respective container.
        """
        pass

    @abstractmethod
    def MapGradient(self, physical_gradient_variable_container_expression_map: 'dict[SupportedSensitivityFieldVariableTypes, ContainerExpressionTypes]') -> ContainerExpressionTypes:
        """Maps physical space gradients to the control space.

        This method is used to map the given physical space gradients to the control space. The input should be as in the following example:
            sensitivity_variable_collective_expression_info = {
                Kratos.YOUNG_MODULUS: Kratos.ContainerExpressions.NodalNonHistoricalContainer,
                Kratos.DENSITY      : Kratos.ContainerExpressions.ElementNonHistoricalContainer,
                Kratos.SHAPE        : Kratos.ContainerExpressions.ElementNonHistoricalContainer
            }

        Args:
            physical_gradient_variable_container_expression_map (dict[SupportedSensitivityFieldVariableTypes, ContainerExpressionTypes]): Map of physical space variable and ContainerExpression with sensitivities.

        Returns:
            ContainerExpressionTypes: Gradients mapped in to control space.
        """
        pass

    @abstractmethod
    def Update(self, control_field: ContainerExpressionTypes) -> bool:
        """Modifies the current control with the given control field.

        This method modifies the control with the given control field from @ref control_field. This
        control_field should be in the control space, and it should not represent a change of design.

        Args:
            control_field (ContainerExpressionTypes): The control field in control space.

        Returns:
            bool: True if the control field was applied to obtain a new design, otherwise False
        """
        pass


