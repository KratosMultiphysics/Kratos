import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes

class ResponseRoutine:
    def __init__(self, master_control: MasterControl, response_name: str, optimization_problem: OptimizationProblem) -> None:
        # set the master control
        self.__master_control = master_control

        # set the response
        self.__response = optimization_problem.GetResponse(response_name)
        self.__response_name = response_name
        self.__response_value = None

        self.__contributing_controls_list: 'list[Control]' = []
        self.__required_physical_gradients: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.ContainerExpression.CollectiveExpressions]' = {}

    def Initialize(self):
        """

        this should create a map of map with keys (control) and values bool stating, whether each control has a contribution
        for this response.

        """
        # create the required physical control fields to compute gradients
        self.__required_physical_gradients = self.__master_control.GetPhysicalKratosVariableCollectiveExpressionsMap()

        # now check which are the dependent physical space variables for the response, if not then remove
        # that variable
        for required_physical_variable in self.__required_physical_gradients.keys():
            if required_physical_variable not in self.__response.GetDependentPhysicalKratosVariables():
                del self.__required_physical_gradients[required_physical_variable]

        for control in self.__master_control.GetListOfControls():
            if set(control.GetPhysicalKratosVariables()).intersection(self.__required_physical_gradients.keys()):
                if Kratos.ModelPartOperationUtilities.HasIntersection(self.__response.GetModelPart(), [self.__response.GetModelPart(), control.GetEmptyControlField().GetModelPart()]):
                    self.__contributing_controls_list.append(control)

        if not self.__contributing_controls_list:
            raise RuntimeError(f"The controls does not have any influence over the response {self.GetReponseName()}.")

    def GetReponseName(self) -> str:
        return self.__response_name

    def CalculateValue(self, control_field: KratosOA.ContainerExpression.CollectiveExpressions) -> float:
        # update using the master control and get updated states.
        update_states = self.__master_control.Update(control_field)

        compute_response_value_flag = False
        for control, is_updated in update_states.items():
            if is_updated and control in self.__contributing_controls_list:
                compute_response_value_flag = True

        compute_response_value_flag = compute_response_value_flag or self.__response_value is None

        # TODO: In the case of having two analysis with the same mesh (model parts) for two different
        #       responses, we need to flag all the anayses which are affected by the control update_state
        #       from the first call, otherwise the second call will not update anything, hence no execution
        #       policies will be executed.

        #       Example: Drag 1 with 10m/s, Drag 2 with 20 m/s using the same mesh and model parts with two
        #                different analyses.
        # # now set the response value in response data container
        # if len(modified_model_parts) > 0:
        #     # if there are any modified model parts, then this changes state of is_executed to false
        #     # in the execution policies which has some intersection with the modified model parts.
        #     ChangeExecutionPolicyStates(modified_model_parts, False, self.__optimization_problem)

        if compute_response_value_flag:
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
        return self.__master_control.MapGradient(self.__required_physical_gradients)

