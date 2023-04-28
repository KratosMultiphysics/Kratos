import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl

class ResponseRoutine:
    def __init__(self, master_control: MasterControl, response_name: str, optimization_info: OptimizationInfo) -> None:
        self.__master_control = master_control
        self.__response_name = response_name
        self.__response = optimization_info.GetResponse(response_name)
        self.__optimization_info = optimization_info
        self.__reponse_value = None
        self.__contributin_controls_list: 'list[Control]' = []

    def Initialize(self):
        """

        this should create a map of map with keys (control) and values bool stating, whether each control has a contribution
        for this response.

        """

        self.__contributin_controls_list: 'list[Control]' = []

        # create the required physical control fields to compute gradients
        self.__required_physical_gradients = self.__master_control.GetPhysicalVariableCollectiveExpressionsMap()

        # now check which are the dependent physical space variables for the response, if not then remove
        # that variable
        for required_physical_variable in self.__required_physical_gradients.keys():
            if required_physical_variable not in self.__response.GetDependentPhysicalKratosVariables():
                del self.__required_physical_gradients[required_physical_variable]        

    def GetReponseName(self) -> str:
        return self.__response_name

    def CalculateValue(self, control_field: KratosOA.ContainerExpression.CollectiveExpressions) -> float:
        # update using the master control and get updated states.
        update_states = self.__master_control.Update(control_field)

        dependent_kratos_variables = self.__response.GetDependentPhysicalKratosVariables()

        # now get modified model parts
        compute_value = False
        for control, is_updated in update_states.items():
            if is_updated and control in self.__contributin_controls_list:
                compute_value = True

        compute_value = compute_value or self.__reponse_value is None

        # now set the response value in response data container
        if len(modified_model_parts) > 0:
            # if there are any modified model parts, then this changes state of is_executed to false
            # in the execution policies which has some intersection with the modified model parts.
            ChangeExecutionPolicyStates(modified_model_parts, False, self.__optimization_info)


        if compute_reponse_value:
            self.__response_value = self.__response.CalculateValue()

        return self.__response_value

    def CalculateGradient(self) -> KratosOA.ContainerExpression.CollectiveExpressions:
        """Returns Collective expression containing all the control space gradients for all control variable types (fields).

        Notes:
            1. It expects that the CalcaulteValue is called.
            2. The gradients are computed with respect to updates from master control.

        Returns:
            KratosOA.ContainerExpression.CollectiveExpressions: _description_
        """
        # fills the proper physical gradients from the response
        self.__response.CalculateGradient(self.__required_physical_gradients)

        # calculate and return the control space gradients from respective controls
        return self.__master_control.MapGradients(self.__required_physical_gradients)

