import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import CallOnAll

class MasterControl:
    def __init__(self) -> None:
        self.__list_of_controls: 'list[Control]' = []

    def AddControl(self, control: Control) -> None:
        self.__list_of_controls.append(control)

    def GetControls(self) -> 'list[Control]':
        return self.__list_of_controls

    def Initialize(self) -> None:
        CallOnAll(self.__list_of_controls, Control.Initialize)

    def Finalize(self) -> None:
        CallOnAll(self.__list_of_controls, Control.Finalize)

    def GetPhysicalVariables(self) -> 'list[list[SupportedSensitivityFieldVariableTypes]]':
        physical_variables: 'list[list[SupportedSensitivityFieldVariableTypes]]' = []

        for control in self.__list_of_controls:
            physical_variables.append(control.GetPhysicalVariables())

        return physical_variables

    def GetControlVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        control_variables: 'list[SupportedSensitivityFieldVariableTypes]' = []

        for control in self.__list_of_controls:
            control_variables.append(control.GetControlVariable())

        return control_variables

    def GetEmptyControlFields(self) -> KratosOA.ContainerExpression.CollectiveExpressions:
        empty_control_fields = KratosOA.ContainerExpression.CollectiveExpressions()

        for control in self.__list_of_controls:
            empty_control_fields.Add(control.GetEmptyControlField())

        return empty_control_fields

    def MapGradients(self, physical_space_sensitivity_variable_and_collective_expressions_map: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.ContainerExpression.CollectiveExpressions]') -> KratosOA.ContainerExpression.CollectiveExpressions:
        mapped_gradients = KratosOA.ContainerExpression.CollectiveExpressions()

        for control in self.__list_of_controls:
            control_physical_sensitivities_container_expression_map = {}
            for physical_control_variable in control.GetPhysicalVariables():
                if physical_control_variable in physical_space_sensitivity_variable_and_collective_expressions_map.keys():
                    control_expression = control.GetEmptyControlField()
                    sensitivity_collective_expression = physical_space_sensitivity_variable_and_collective_expressions_map[physical_control_variable]
                    for container_expression in sensitivity_collective_expression.GetContainerExpressions():
                        if control_expression.GetContainer() == container_expression.GetContainer():
                            control_expression.CopyFrom(container_expression)
                            control_physical_sensitivities_container_expression_map[physical_control_variable] = control_expression
                            break

            mapped_gradients.Add(control.MapGradient(control_physical_sensitivities_container_expression_map))

        return mapped_gradients

    def Update(self, update_collective_expressions: KratosOA.ContainerExpression.CollectiveExpressions) -> None:
        if len(self.__list_of_controls) != len(update_collective_expressions.GetContainerExpressions()):
            raise RuntimeError(f"Controls size and update size mismatch [ number of controls: {len(self.__list_of_controls)}, number of container expressions: {len(update_collective_expressions.GetContainerExpressions())} ].")

        for control, container_expression in zip(self.__list_of_controls, update_collective_expressions.GetContainerExpressions()):
            control.Update(container_expression)

