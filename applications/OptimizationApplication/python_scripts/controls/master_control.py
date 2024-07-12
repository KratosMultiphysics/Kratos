import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import IsSameContainerExpression
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import HasContainerExpression
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import time_decorator
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import CallOnAll

class MasterControl:
    """Master control class.

    This class is used to simplify working with many controls at once. Following responsibilities are assumed:
        1. Maps physical gradients from different domains to one CollectiveExpression (using MapGradient).
        2. Updates each respective domain from updates given by one CollectiveExpression (using Update).

    There should be only one master control class per optimization problem.

    """
    def __init__(self) -> None:
        self.__list_of_controls: 'list[Control]' = []

    def AddControl(self, control: Control) -> None:
        """Adds a given control to the master control.

        Args:
            control (Control): Control to be added
        """
        self.__list_of_controls.append(control)

    def GetListOfControls(self) -> 'list[Control]':
        """Returns the list of controls in the master control.

        Returns:
            list[Control]: List of controls.
        """
        return self.__list_of_controls

    def GetPhysicalKratosVariableCollectiveExpressionsMap(self) -> 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]':
        """Returns map of physical variables and collective expressions from each control.

        This returns a map of physical control variables and a collective expressions. The collective expressions will contain
        all the container expressions for respective control's control domains for each physical control variable. The repeated container
        expressions for the same physical control variable is omitted to avoid double calculations.

        Returns:
            dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]: Physical control variable and collective expressions map.
        """
        physical_variable_collective_expressions: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]' = {}

        for control in self.__list_of_controls:
            for physical_variable in control.GetPhysicalKratosVariables():
                # check whether the physical variable is already there.
                if not physical_variable in physical_variable_collective_expressions.keys():
                    physical_variable_collective_expressions[physical_variable] = KratosOA.CollectiveExpression()

                current_variable_collective_expression = physical_variable_collective_expressions[physical_variable]

                # check whether the container for that physical variable is already there.
                control_container_expression = control.GetEmptyField()
                if not HasContainerExpression(control_container_expression, current_variable_collective_expression.GetContainerExpressions()):
                    current_variable_collective_expression.Add(control_container_expression)

        return physical_variable_collective_expressions

    def GetEmptyField(self) -> KratosOA.CollectiveExpression:
        """Returns empty CollectiveExpression containing empty ContainerExpressions for each control.

        Returns:
            KratosOA.CollectiveExpression: Empty CollectiveExpression
        """
        empty_control_fields = KratosOA.CollectiveExpression()

        for control in self.__list_of_controls:
            empty_control_fields.Add(control.GetEmptyField())

        return empty_control_fields

    def GetControlField(self) -> KratosOA.CollectiveExpression:
        """Returns CollectiveExpression containing control field ContainerExpressions for each control.

        Returns:
            KratosOA.CollectiveExpression: Control field CollectiveExpression
        """
        control_fields = KratosOA.CollectiveExpression()

        for control in self.__list_of_controls:
            control_fields.Add(control.GetControlField())

        return control_fields

    @time_decorator()
    def MapGradient(self, physical_space_gradient_variable_and_collective_expressions_map: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]') -> KratosOA.CollectiveExpression:
        """Maps physical space gradients to a collective expression.

        This method maps sensitivities w.r.t. physical space variables to control space by using each control. It is done by converting input
        map with CollectiveExpression to a map with ContainerExpressions. The following is a pseudo example:

        input:

        physical_space_gradient_variable_and_collective_expressions_map = {
            "YOUND_MODULUS": [ControlDomain1, ControlDomain2],
            "DENSITY"      : [ControlDomain1],
            "VISCOSITY"    : [ControlDomain2],
        }

        control info:
        control1: domain = ControlDomain1
                  physical_vars = YOUND_MODULUS, DENSITY

        control2: domain = ControlDomain2
                  physical_vars = VISCOSITY, DENSITY


        from above information, following two maps are created and passed to each control to obtain one ContainerExpression from each control.
        control_specific_maps:
        for control1: {YOUND_MODULUS: ControlDomain1, DENSITY: ControlDomain1}
        for control2: {VISCOSITY: ControlDomain2, DENSITY: A zero valued ControlDomain}

        then returned mapped gradients are added to one CollectiveExpression.

        In here, the missing gradients for required physical variables will be assumed to be zero.

        This converts given map of sensitivities w.r.t. different physical space variables to one sensitivities in control space and aggregated
        to one CollectiveExpression.

        Args:
            physical_space_gradient_variable_and_collective_expressions_map (dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]): Map of physical space sensitivities.

        Returns:
            KratosOA.CollectiveExpression: Control space sensitivities.
        """


        mapped_gradients = KratosOA.CollectiveExpression()

        for control in self.__list_of_controls:
            # iterate through each control to create its own container expression map from the collective expressions map given as input
            control_physical_sensitivities_container_expression_map = {}
            for physical_control_variable in control.GetPhysicalKratosVariables():
                # first assume the gradients for this physical_control_variable is zero, hence get the zero valued expression.
                control_expression = control.GetEmptyField()

                # get the required physical variables from control.
                if physical_control_variable in physical_space_gradient_variable_and_collective_expressions_map.keys():
                    # if the sensitivities for the given physical control variable exists, then try to find whether
                    # for this specific control the physical control variable sensitivities exists.

                    sensitivity_collective_expression = physical_space_gradient_variable_and_collective_expressions_map[physical_control_variable]
                    for container_expression in sensitivity_collective_expression.GetContainerExpressions():
                        if IsSameContainerExpression(control_expression, container_expression):
                            # there exists for this control's physical variables sensitivities. then copy it to the expression
                            # this copy moves the underlying vectors, hence cheap.
                            control_expression.CopyFrom(container_expression)
                            break

                # now add it to the map. If it is found from input gradients, the control_expression will have those values,
                # otherwise it will have representative zero control_expression.
                control_physical_sensitivities_container_expression_map[physical_control_variable] = control_expression

            # map the physical control variable sensitivities to one control space
            mapped_gradients.Add(control.MapGradient(control_physical_sensitivities_container_expression_map))

        return mapped_gradients

    @time_decorator()
    def Update(self, update_collective_expressions: KratosOA.CollectiveExpression) -> 'dict[Control, bool]':
        """Update each control with given collective expression's respective container expression.

        Args:
            update_collective_expressions (KratosOA.CollectiveExpression): Update

        Raises:
            RuntimeError: If number of controls and number of container expressions mismatch.

        Returns:
            dict[Control, bool]: A map with control and a boolean whether the update changed anything in that control.
        """
        if len(self.__list_of_controls) != len(update_collective_expressions.GetContainerExpressions()):
            raise RuntimeError(f"Controls size and update size mismatch [ number of controls: {len(self.__list_of_controls)}, number of container expressions: {len(update_collective_expressions.GetContainerExpressions())} ].")

        update_map: 'dict[Control, bool]' = {}
        for control, container_expression in zip(self.__list_of_controls, update_collective_expressions.GetContainerExpressions()):
            update_map[control] = control.Update(container_expression)

        return update_map

    def Check(self) -> None:
        CallOnAll(self.__list_of_controls, Control.Check)

    def Initialize(self) -> None:
        CallOnAll(self.__list_of_controls, Control.Initialize)

    def Finalize(self) -> None:
        CallOnAll(self.__list_of_controls, Control.Finalize)
