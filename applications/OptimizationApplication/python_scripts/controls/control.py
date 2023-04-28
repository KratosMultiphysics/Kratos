from abc import ABC, abstractmethod
from typing import Any
import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes

class Control(ABC):
    """Base abstract control class.

    The base abstract control class has the following responsibilities.
        1. Initalizes and finalizes the control.
        2. Mapping physical space gradients to control space gradients.
        2. Updating the controlled parts of the model part with final design given in control space.
        3. Retrieve control variable.
        4. Retrieve physical space variables (If more than one physical space variables are controlled by the given control variable.)

    This control should only work on one model part and one control variable. Hence, if required multiple model parts can be used
    by using the union operator in Kratos.ModelPartBinaryOperators to create one model part within the control.

    """
    def Initialize(self) -> None:
        """Initializes the control.

        This method initializes the control. This is only called once in the whole optimization process.

        """
        pass

    def Finalize(self) -> None:
        """Finalizes the control.

        This method finalizes the control. This is only called once in the whole optimization process.

        """
        pass

    def GetPhysicalVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        """Returns list of physical control variables controlled by this control.

        This method returns list of physical variable which are controlled by the control. In most of the cases, there is only one
        control variable per control. But in some rare cases, there may be the requirement to control two physical variables. Anyhow,
        the control space should only have one variable per control in the control space, but they can have multiple physical variables
        in their design space.

        If the case is that this control has one physical variable, then it should be same as the control variable which is the
        default behaviour.

        Returns:
            list[SupportedSensitivityFieldVariableTypes]: List of physical control variables.
        """
        return [self.GetControlVariable()]

    @staticmethod
    @abstractmethod
    def Create(model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo) -> Any:
        """Returns a new control constructed with standard input parameters.

        Args:
            model (Kratos.Model): Kratos model the control belongs to.
            parameters (Kratos.Parameters): Parameters to define behaviour of the control.
            optimization_info (OptimizationInfo): Optimization problem data.

        Returns:
            Any: The constructed control.
        """
        pass

    @abstractmethod
    def GetControlVariable(self) -> SupportedSensitivityFieldVariableTypes:
        """Returns the variable used in this control.

        Returns:
            SupportedSensitivityFieldVariableTypes: Variable used in this control.
        """
        pass

    @abstractmethod
    def GetEmptyControlField(self) -> ContainerExpressionTypes:
        """Returns an empty control field data holder.

        This returns an empty control field data holder to give information about on which model part's container
        this model part is acting on. This has O(1) complexity, hence has the least cost because it does not read
        any data from respective model part's container.

        Returns:
            ContainerExpressionTypes: Returns an empty ContainerExpression corresponding to control's model part's respective container.
        """
        pass

    @abstractmethod
    def MapGradient(self, sensitivity_variable_container_expression_map: 'dict[SupportedSensitivityFieldVariableTypes, ContainerExpressionTypes]') -> ContainerExpressionTypes:
        """Maps physical space sensitivity gradients to the control space

        This method is used to map the given physical space gradients to the control space. The input should be as in the following example:
            sensitivity_variable_collective_expression_info = {
                Kratos.YOUNG_MODULUS: Kratos.ContainerExpressions.NodalNonHistoricalContainer,
                Kratos.DENSITY      : Kratos.ContainerExpressions.ElementNonHistoricalContainer,
                Kratos.SHAPE        : Kratos.ContainerExpressions.ElementNonHistoricalContainer
            }

        Args:
            sensitivity_variable_container_expression_map (dict[SupportedSensitivityFieldVariableTypes, ContainerExpressionTypes]): Map of design variable and ContainerExpression with sensitivities

        Returns:
            ContainerExpressionTypes: Gradients mapped in to control space.
        """
        pass

    @abstractmethod
    def Update(self, update_container_expression: ContainerExpressionTypes) -> None:
        """Updates the current control.

        This method updates the control with the given update from @ref update_container_expression. This
        update should be in the control space, and it should be the final design, not a change of design.

        Args:
            update_container_expression (ContainerExpressionTypes): Final design in control space as an update.
        """
        pass


