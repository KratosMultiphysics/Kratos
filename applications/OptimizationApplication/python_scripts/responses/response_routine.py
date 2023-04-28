import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl

class ResponseRoutine:
    def __init__(self, master_control: MasterControl, response_name: str, optimization_info: OptimizationInfo) -> None:
        self.__master_control = master_control
        self.__response_name = response_name
        self.__response = optimization_info.GetResponse(response_name)
        self.__optimization_info = optimization_info

    def CalculateValue(self, control_space_updates: KratosOA.ContainerExpression.CollectiveExpressions) -> float:
        # update using the master control and get updated states.
        update_states = self.__master_control.Update(control_space_updates)

        # now get modified model parts
        modified_model_parts: 'list[Kratos.ModelPart]' = []
        for control, is_updated in update_states.items():
            if is_updated:
                modified_model_parts.append(control.GetEmptyControlField().GetModelPart())

        response_data = self.__optimization_info.GetReponseData(self.__response_name)

        # now set the response value in response data container
        if len(modified_model_parts) > 0:
            # if there are any modified model parts, then this changes state of is_executed to false
            # in the execution policies which has some intersection with the modified model parts.
            ChangeExecutionPolicyStates(modified_model_parts, False, self.__optimization_info)
            response_data.SetValue("value", self.__response.CalculateValue())

        if not response_data.HasValue("value"):
            # this can happen, if there is no change in model part and we need the response value
            # for a new design iteration, but the response is already converged.
            # then we try to get the response value from the previous iteration.
            if response_data.HasValue("value", 1):
                response_data.SetValue("value", response_data.GetValue("value", 1))
            else:
                # this block is reached usually in the initial iteration where no model parts are updated
                # as well as no value is found for this response in current or previous iteration.
                # then we need to calculate the response value.
                response_data.SetValue("value", self.__response.CalculateValue())

        # return response value from the container
        return response_data.GetValue("value")

    def CalculateGradient(self, control_space_updates: KratosOA.ContainerExpression.CollectiveExpressions) -> KratosOA.ContainerExpression.CollectiveExpressions:
        #
        # update using the master control and get updated states.
        update_states = self.__master_control.Update(control_space_updates)

        # now get modified model parts
        modified_model_parts: 'list[Kratos.ModelPart]' = []
        for control, is_updated in update_states.items():
            if is_updated:
                modified_model_parts.append(control.GetEmptyControlField().GetModelPart())

        # create the required physical control fields to compute gradients
        required_physical_gradients = self.__master_control.GetPhysicalVariableCollectiveExpressionsMap()

        # now check which are the dependent physical space variables for the response, if not then remove
        # that variable
        for required_physical_variable in required_physical_gradients.keys():
            if required_physical_variable not in self.__response.GetDependentPhysicalVariable():
                del required_physical_gradients[required_physical_variable]

        response_data = self.__optimization_info.GetReponseData(self.__response_name)

        # now set the response value in response data container
        if len(modified_model_parts) > 0:
            # if there are any modified model parts, then this changes state of is_executed to false
            # in the execution policies which has some intersection with the modified model parts.
            ChangeExecutionPolicyStates(modified_model_parts, False, self.__optimization_info)
            response_data.SetValue("value", self.__response.CalculateValue())