import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.utilities.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.utilities.communicators.optimization_component_communicator import OptimizationComponentCommunicator
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes

class ResponseFunctionCommunicator:
    def __init__(self, name: str, optimization_info: OptimizationInfo, required_buffer_size: int = 2):
        # minimum buffer size of 2 is required for the communicator
        # algorithm can request higher buffer sizes if required.
        required_buffer_size = max(required_buffer_size, 2)

        self.__optimization_component_communicator = OptimizationComponentCommunicator(optimization_info)
        problem_data = self.__optimization_component_communicator.GetProblemData()

        if not problem_data.HasValue(f"responses/{name}"):
            # first create the buffered data containers to store light objects
            # such as primitive variables
            self.__response_problem_data = OptimizationInfo(required_buffer_size)
            problem_data[f"responses/{name}"] = self.__response_problem_data

            # now create the sensitivities data container without buffer because
            # they store heavy objects containing sensitivities, and buffer
            # is not required.
            self.__response_sensitivity_data = OptimizationInfo(1)
            problem_data[f"responses/{name}/sensitivities"] = self.__response_sensitivity_data
        else:
            self.__response_problem_data: OptimizationInfo = problem_data[f"responses/{name}"]
            if self.__response_problem_data.GetBufferSize() < required_buffer_size:
                raise RuntimeError(f"The required buffer size is not satisfied with the existing problem data container. [ response data container buffer size = {self.__response_problem_data.GetBufferSize()}, required buffer size = {required_buffer_size}")

            self.__response_sensitivity_data: OptimizationInfo = problem_data[f"responses/{name}/sensitivities"]

        self.__response_function = self.__optimization_component_communicator.GetResponseFunction(name)
        self.__name = name

    def GetName(self) -> str:
        return self.__name

    def GetResponseFunction(self) -> ResponseFunction:
        return self.__response_function

    def GetResponseProblemData(self) -> OptimizationInfo:
        return self.__response_problem_data

    def GetResponseSensitivityData(self) -> OptimizationInfo:
        return self.__response_sensitivity_data

    def GetScaledValue(self, step_index: int = 0, scaling: float = 1.0) -> float:
        response_problem_data = self.GetResponseProblemData()
        if step_index == 0 and not response_problem_data.HasValue("value"):
            response_problem_data["value"] = self.__response_function.CalculateValue()

        return response_problem_data["value"] * scaling

    def GetRelativeChange(self) -> float:
        return self.GetScaledValue() / self.GetScaledValue(1) - 1.0 if self.__optimization_component_communicator.GetStep() > 1 else 0.0

    def GetAbsoluteChange(self, reference_value: float) -> float:
        return self.GetScaledValue() / reference_value - 1.0 if abs(reference_value) > 1e-12 else self.GetScaledValue()

    def CalculateScaledSensitivity(self, sensitivity_variable_collective_expression_info: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.ContainerExpression.CollectiveExpressions]', scaling: float = 1.0) -> None:
        response_sensitivity_data = self.GetResponseSensitivityData()

        # check whether the same sensitivity is computed for each sensitivity model part w.r.t. same sensitivity variable for the same response function.
        required_sensitivity_model_part_variable_info: 'dict[SupportedSensitivityFieldVariableTypes, list[Kratos.ModelPart]]' = {}
        for sensitivity_variable, collective_expression in sensitivity_variable_collective_expression_info.items():
            is_already_computed = True
            for container_expression in collective_expression.GetContainerExpressions():
                key = f"{container_expression.GetModelPart().FullName()}/raw_{self.GetName()}_{sensitivity_variable.Name()}"
                is_already_computed = is_already_computed and response_sensitivity_data.HasValue(key)

            if not is_already_computed:
                # translates container expressions to model parts which is known by responses.
                required_sensitivity_model_part_variable_info[sensitivity_variable] = [container_expression.GetModelPart() for container_expression in collective_expression.GetContainerExpressions()]

        self.__response_function.CalculateSensitivity(required_sensitivity_model_part_variable_info)

        # now store the computed sensitivities from resepective containers to container expressions
        for sensitivity_variable in required_sensitivity_model_part_variable_info.keys():
            for container_expression in sensitivity_variable_collective_expression_info[sensitivity_variable].GetContainerExpressions():
                key = f"{container_expression.GetModelPart().FullName()}/raw_{self.GetName()}_{sensitivity_variable.Name()}"
                container_expression.Read(sensitivity_variable)
                # this cloning is cloning of an expression pointer, hence it is not expensive.
                # no cloning of underlying sensitivity data is involved.
                response_sensitivity_data[key] = container_expression.Clone()

        # now assign the sensitivities for all the collective expressions correctly
        for sensitivity_variable, collective_expression in sensitivity_variable_collective_expression_info.items():
            for container_expression in collective_expression.GetContainerExpressions():
                key = f"{container_expression.GetModelPart().FullName()}/raw_{self.GetName()}_{sensitivity_variable.Name()}"
                # this copying is copying of an expression pointer, hence it is not expensive.
                # no copying of underlying sensitivity data is involved.
                container_expression.CopyFrom(response_sensitivity_data[key])
                container_expression *= scaling
