from typing import Optional

import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ResponseFunction:
    if not parameters.Has("name"):
        raise RuntimeError(f"CombinedResponseFunction instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"CombinedResponseFunction instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")

    return CombinedResponseFunction(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)

class CombinedResponseFunction(ResponseFunction):
    """Response function to combine multiple response functions.

    This response function can combine multiple response functions to a single response function.
    Following methods are used to combine the response values:
              1. weighted sum - Weighted summation of responses list is used as the combined response function value.

    """
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        super().__init__(name)

        default_settings = Kratos.Parameters("""{
            "combining_method"   : "sum",
            "combining_responses": [
                {
                    "response_name"  : "",
                    "response_weight": 0.0
                }
            ]
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)
        self.combining_method = parameters["combining_method"].GetString()

        self.model = model
        self.model_part: Optional[Kratos.ModelPart] = None
        self.optimization_problem = optimization_problem

        self.list_of_responses: 'list[tuple[ResponseFunction, float]]' = []
        for response_params in parameters["combining_responses"].values():
            response_params.ValidateAndAssignDefaults(default_settings["combining_responses"].values()[0])
            response_name = response_params["response_name"].GetString()

            if response_name not in [response.GetName() for response in optimization_problem.GetListOfResponses()]:
                raise RuntimeError(f"\"{response_name}\" not found in the optimization problem. Please check whether this response is defined before the \"{self.GetName()}\".")

            self.list_of_responses.append(
                (
                    self.optimization_problem.GetResponse(response_name),
                    response_params["response_weight"].GetDouble()
                ))

        if len(self.list_of_responses) == 0:
            raise RuntimeError(f"Combined response \"{self.GetName()}\" does not have any responses to combine.")

    def GetImplementedPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        variables_list: 'list[SupportedSensitivityFieldVariableTypes]' = []
        for response, _ in self.list_of_responses:
            variables_list.extend(response.GetImplementedPhysicalKratosVariables())
        return list(set(variables_list))

    def GetInfluencingModelPart(self) -> Kratos.ModelPart:
        influencing_model_part_names = [response.GetInfluencingModelPart().FullName() for response, _ in self.list_of_responses]
        model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"response_{self.GetName()}", influencing_model_part_names, False)
        self.model_part = model_part_operation.GetModelPart()
        return self.model_part

    def Initialize(self) -> None:
        for response, _ in self.list_of_responses:
            response.Initialize()

        component_data_view = ComponentDataView(self, self.optimization_problem)
        component_data_view.SetDataBuffer(1)
        self.buffered_data = component_data_view.GetBufferedData()
        self.unbuffered_data = component_data_view.GetUnBufferedData()

    def Check(self) -> None:
        for response, _ in self.list_of_responses:
            response.Check()

    def Finalize(self) -> None:
        for response, _ in self.list_of_responses:
            response.Finalize()

    def CalculateValue(self) -> float:
        if self.combining_method == "sum":
            value = 0.0
            for response, weight in self.list_of_responses:
                    response_value = response.CalculateValue()
                    value += response_value * weight
                    self.buffered_data.SetValue(f"{response.GetName()}", response_value, overwrite=True)
        else:
            raise RuntimeError(f"Unsupported combining_method = \"{self.combining_method}\". Followings are supported:\n\tsum")
        return value

    def CalculateGradient(self, physical_variable_combined_tensor_adaptor: 'dict[SupportedSensitivityFieldVariableTypes, Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor]') -> None:
        # make everything zeros
        for _, combined_tensor_adaptor in physical_variable_combined_tensor_adaptor.items():
            combined_tensor_adaptor.data[:] = 0.0
            Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(combined_tensor_adaptor, perform_store_data_recursively=False, copy=False).StoreData()

        for response_function, weight in self.list_of_responses:
            # get the map for the response function
            sub_physical_variable_combined_tensor_adaptor: 'dict[SupportedSensitivityFieldVariableTypes, Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor]' = {}
            for physical_variable, combined_tensor_adaptor in physical_variable_combined_tensor_adaptor.items():
                if physical_variable in response_function.GetImplementedPhysicalKratosVariables():
                    sub_physical_variable_combined_tensor_adaptor[physical_variable] = combined_tensor_adaptor.Clone()

            # calculate the gradients
            response_function.CalculateGradient(sub_physical_variable_combined_tensor_adaptor)

            # now aggregate the gradients
            for physical_variable, local_combined_tensor_adaptor in sub_physical_variable_combined_tensor_adaptor.items():
                self.unbuffered_data.SetValue(f"d{response_function.GetName()}_d{physical_variable.Name()}", local_combined_tensor_adaptor.Clone(), overwrite=True)
                if self.combining_method == "sum":
                    physical_variable_combined_tensor_adaptor[physical_variable].data[:] += local_combined_tensor_adaptor.data[:] * weight
                else:
                    raise RuntimeError(f"Unsupported combining_method = \"{self.combining_method}\". Followings are supported:\n\tsum")

        for _, combined_tensor_adaptor in physical_variable_combined_tensor_adaptor.items():
            Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(combined_tensor_adaptor, perform_store_data_recursively=False, copy=False).StoreData()