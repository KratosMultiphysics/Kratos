import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.responses.response_routine import ResponseRoutine
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView

class StandardizedObjective(ResponseRoutine):
    """Standardized objective response function

    This class creates instances to standardize any response function (make it a minimization problem).
    Supported objective types:
        "minimization",
        "maximization"

    """
    def __init__(self, parameters: Kratos.Parameters, master_control: MasterControl, optimization_problem: OptimizationProblem, required_buffer_size: int = 2):
        default_parameters = Kratos.Parameters("""{
            "response_name": "",
            "type"         : "",
            "scaling"      : 1.0
        }""")
        parameters.ValidateAndAssignDefaults(default_parameters)

        response = optimization_problem.GetResponse(parameters["response_name"].GetString())

        super().__init__(master_control, response)

        if required_buffer_size < 2:
            raise RuntimeError(f"Standardized objective requires 2 as minimum buffer size. [ response name = {self.GetReponse().GetName()} ]")

        self.__optimization_problem = optimization_problem
        self.__component_data_view = ComponentDataView(response, optimization_problem)
        self.__component_data_view.SetDataBuffer(required_buffer_size)

        scaling = parameters["scaling"].GetDouble()
        if scaling < 0.0:
            raise RuntimeError(f"Scaling should be always positive [ given scale = {scaling}]")

        self.__objective_type = parameters["type"].GetString()
        if self.__objective_type == "minimization":
            self.__scaling = scaling
        elif self.__objective_type == "maximization":
            self.__scaling = -scaling
        else:
            raise RuntimeError(f"Requesting unsupported type {self.__objective_type} for objective response function. Supported types are: \n\tminimization\n\tmaximization")

        self.__initial_response_value = None

    def GetInitialValue(self) -> float:
        if self.__initial_response_value is None:
            return self.__initial_response_value
        else:
            raise RuntimeError(f"Response value for {self.GetReponse().GetName()} is not calculated yet.")

    def CalculateStandardizedValue(self, control_field: KratosOA.ContainerExpression.CollectiveExpressions, save_value: bool = True) -> float:
        response_value = self.CalculateValue(control_field)
        standardized_response_value = response_value * self.__scaling

        if not self.__initial_response_value:
            self.__initial_response_value = standardized_response_value

        if save_value:
            data = self.__component_data_view.GetBufferedData()
            data["value"] = response_value
            data["standardized_value"] = standardized_response_value

        return standardized_response_value

    def GetValue(self, step_index: int) -> float:
        return self.__component_data_view.GetBufferedData().GetValue("value", step_index)

    def GetStandardizedValue(self, step_index: int) -> float:
        return self.__component_data_view.GetBufferedData().GetValue("standardized_value", step_index)

    def CalculateStandardizedGradient(self, save_value: bool = True) -> KratosOA.ContainerExpression.CollectiveExpressions:
        gradient_collective_expression = self.CalculateGradient()

        if save_value:
            for gradient_container_expression, control in zip(gradient_collective_expression.GetContainerExpressions(), self.GetMasterControl().GetListOfControls()):
                self.__component_data_view.GetUnBufferedData()[f"d{self.GetReponse().GetName()}_d{control.GetName()}"] = gradient_container_expression.Clone()

        return gradient_collective_expression * self.__scaling

    def GetRelativeChange(self) -> float:
        if self.__optimization_problem.GetStep() > 1:
            return self.GetStandardizedValue() / self.GetStandardizedValue(1) - 1.0 if abs(self.GetStandardizedValue(1)) > 1e-12 else self.GetStandardizedValue()
        else:
            return 0.0

    def GetAbsoluteChange(self) -> float:
        return self.GetStandardizedValue() / self.GetInitialValue() - 1.0 if abs(self.GetInitialValue()) > 1e-12 else self.GetStandardizedValue()

    def UpdateOptimizationProblemData(self) -> None:
        response_problem_data = self.__component_data_view.GetBufferedData()
        response_problem_data["rel_change"] = self.GetRelativeChange()
        response_problem_data["abs_change"] = self.GetAbsoluteChange()

    def GetInfo(self) -> str:
        msg = "\tObjective info:"
        msg += f"\n\t\t name          : {self.GetReponse().GetName()}"
        msg += f"\n\t\t type          : {self.__objective_type}"
        msg += f"\n\t\t value         : {self.GetValue():0.6e}"
        msg += f"\n\t\t abs_change    : {self.GetAbsoluteChange():0.6e}"
        msg += f"\n\t\t rel_change [%]: {self.GetRelativeChange() * 100.0:0.6e}"
        return msg