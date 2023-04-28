import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.responses.response_routine import ResponseRoutine
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl

class StandardizedObjective(ResponseRoutine):
    """Standardized objective response function

    This class creates instances to standardize any response function (make it a minimization problem).
    Supported objective types:
        "minimization",
        "maximization"

    """
    def __init__(self, parameters: Kratos.Parameters, master_control: MasterControl, optimization_info: OptimizationProblem, required_buffer_size: int = 2):
        default_parameters = Kratos.Parameters("""{
            "response_name": "",
            "type"         : "",
            "scaling"      : 1.0
        }""")
        parameters.ValidateAndAssignDefaults(default_parameters)

        super().__init__(master_control, parameters["response_name"].GetString(), optimization_info, required_buffer_size)

        # TODO: setting the buffer sizes correclty

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
            raise RuntimeError(f"Response value for {self.GetReponseName()} is not calculated yet.")
    
    def CalculateStandardizedValue(self, control_space_updates: KratosOA.ContainerExpression.CollectiveExpressions, save_value: bool = True) -> float:
        standardized_response_value = self.CalculateValue(control_space_updates) * self.__scaling
        if not self.__initial_response_value:
            self.__initial_response_value = standardized_response_value

        # put the saving mechanism

        return standardized_response_value

    def GetReponseValue(self, step_index):
        # return the value
        pass

    def CalculateStandardizedGradient(self) -> KratosOA.ContainerExpression.CollectiveExpressions:
        self.CalculateGradient() * self.__scaling

    def UpdateObjectiveData(self) -> None:
        response_problem_data = self.__utility.GetBufferedData()
        response_problem_data["type"] = self.GetType()
        response_problem_data["rel_change"] = self.__utility.GetRelativeChange()
        response_problem_data["abs_change"] = self.__utility.GetAbsoluteChange(self.GetInitialValue())

    def GetObjectiveInfo(self) -> str:
        msg = "\tObjective info:"
        msg += f"\n\t\t name          : {self.GetName()}"
        msg += f"\n\t\t type          : {self.GetType()}"
        msg += f"\n\t\t value         : {self.__utility.GetScaledValue():0.6e}"
        msg += f"\n\t\t rel_change [%]: {self.__utility.GetRelativeChange() * 100.0:0.6e}"
        msg += f"\n\t\t abs_change [%]: {self.__utility.GetAbsoluteChange(self.GetInitialValue()) * 100.0:0.6e}"
        return msg