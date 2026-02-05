import numpy
from typing import Optional

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, _) -> ResponseFunction:
    if not parameters.Has("name"):
        raise RuntimeError(f"DiscreteValueResidualResponseFunction instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"DiscreteValueResidualResponseFunction instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")

    return DiscreteValueResidualResponseFunction(parameters["name"].GetString(), model, parameters["settings"])

class DiscreteValueResidualResponseFunction(ResponseFunction):
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters):
        super().__init__(name)

        default_settings = Kratos.Parameters("""{
            "evaluated_model_part_names"     : [
                "PLEASE_PROVIDE_A_MODEL_PART_NAME"
            ],
            "container_type"         : "node_historical",
            "variable_name"          : "",
            "residual_type"          : "exact",
            "list_of_discrete_values": [0.0]
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.model = model
        self.variable: SupportedSensitivityFieldVariableTypes = Kratos.KratosGlobals.GetVariable(parameters["variable_name"].GetString())
        if not isinstance(self.variable, Kratos.DoubleVariable):
            raise RuntimeError(f"Only supports double variables. Provided variable is \"{self.variable.Name()}\".")

        self.list_of_discrete_values = parameters["list_of_discrete_values"].GetVector()

        self.residual_type = parameters["residual_type"].GetString()
        if self.residual_type  not in ["exact", "logarithm"]:
            raise RuntimeError(f"Unsupported residual_type = \"{self.residual_type}\" requested. Followings are supported:\n\texact\n\tlogarithm")

        container_type = parameters["container_type"].GetString()
        if container_type == "node_historical":
            self.ta_getter = lambda model_part: Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor(model_part.Nodes, self.variable)
        elif container_type == "node_non_historical":
            self.ta_getter = lambda model_part : Kratos.TensorAdaptors.VariableTensorAdaptor(model_part.Nodes, self.variable)
        elif container_type == "condition":
            self.ta_getter = lambda model_part : Kratos.TensorAdaptors.VariableTensorAdaptor(model_part.Conditions, self.variable)
        elif container_type == "condition_properties":
            self.ta_getter = lambda model_part : KratosOA.TensorAdaptors.PropertiesVariableTensorAdaptor(model_part.Conditions, self.variable)
        elif container_type == "element":
            self.ta_getter = lambda model_part : Kratos.TensorAdaptors.VariableTensorAdaptor(model_part.Elements, self.variable)
        elif container_type == "element_properties":
            self.ta_getter = lambda model_part : KratosOA.TensorAdaptors.PropertiesVariableTensorAdaptor(model_part.Elements, self.variable)
        else:
            raise RuntimeError(f"Unsupported container_type = \"{container_type}\" requested. Followings are supported:\n\tnode_historical\n\tnode_non_historical\n\tcondition\n\tcondition_properties\n\telement\n\telement_properties")

        evaluated_model_part_names = parameters["evaluated_model_part_names"].GetStringArray()
        if len(evaluated_model_part_names) == 0:
            raise RuntimeError(f"No model parts were provided for DiscreteValueResidualResponseFunction. [ response name = \"{self.GetName()}\"]")

        self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"response_{self.GetName()}", evaluated_model_part_names, False)
        self.model_part: Optional[Kratos.ModelPart] = None

    def GetImplementedPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [self.variable]

    def GetInfluencingModelPart(self) -> Kratos.ModelPart:
        if self.model_part is None:
            raise RuntimeError("Please call DiscreteValueResidualResponseFunction::Initialize first.")
        return self.model_part

    def Initialize(self) -> None:
        self.model_part = self.model_part_operation.GetModelPart()

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    def CalculateValue(self) -> float:
        ta = self.ta_getter(self.model_part)
        ta.CollectData()

        resultant = Kratos.TensorAdaptors.DoubleTensorAdaptor(ta)

        if self.residual_type == "exact":
            resultant.data[:] = 1.0
        elif self.residual_type == "logarithm":
            resultant.data[:] = 0.0

        for value in self.list_of_discrete_values:
            if self.residual_type == "exact":
                resultant.data[:] *= (ta.data[:] - value) ** 2
            elif self.residual_type == "logarithm":
                resultant.data[:] += numpy.log((ta.data[:] - value) ** 2)

        return numpy.sum(resultant.data)

    def CalculateGradient(self, physical_variable_combined_tensor_adaptor: 'dict[SupportedSensitivityFieldVariableTypes, Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor]') -> None:
        values = self.ta_getter(self.model_part)
        values.CollectData()

        # calculate the gradients
        for physical_variable, cta in physical_variable_combined_tensor_adaptor.items():
            if physical_variable == self.variable:

                # initialize the current tensor adaptor
                for ta in cta.GetTensorAdaptors():
                    ta.data[:] = 0.0

                for ta in cta.GetTensorAdaptors():
                    for i, value_i in enumerate(self.list_of_discrete_values):
                        if self.residual_type == "exact":
                            partial_gradient_ta = Kratos.TensorAdaptors.DoubleTensorAdaptor(ta)
                            partial_gradient_ta.data[:] = 1.0
                            for j, value_j in enumerate(self.list_of_discrete_values):
                                if i == j:
                                    partial_gradient_ta.data *= (values.data[:] - value_j) * 2.0
                                else:
                                    partial_gradient_ta.data *= (values.data[:] - value_j) ** 2
                            ta.data[:] += partial_gradient_ta.data
                        elif self.residual_type == "logarithm":
                            ta.data[:] += (((values.data - value_i) ** (-2)) * (values.data - value_i) * 2.0)

                Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(cta, perform_collect_data_recursively=False, copy=False).CollectData()
            else:
                raise RuntimeError(f"Unsupported sensitivity w.r.t. {physical_variable.Name()} requested. Followings are supported sensitivity variables:\n\t{self.variable.Name()}")
