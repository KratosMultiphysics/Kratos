import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
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

        if len(modified_model_parts) > 0:
            # if there are any modified model parts, then this changes state of is_executed in
            # the execution policies which has some intersection with the modified
            #  model parts to false.
            ChangeExecutionPolicyStates(modified_model_parts, False, self.__optimization_info)
            return self.__response.CalculateValue()
        else:
            # nothing has changed in the model parts. Hence
            if response_data.HasValue("value"):
                return response_data.GetValue("value")
            else:



