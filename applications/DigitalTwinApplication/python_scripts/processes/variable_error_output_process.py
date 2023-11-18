import typing
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.HDF5Application as KratosHDF5
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.HDF5Application.core.operations.model_part import *
from KratosMultiphysics.HDF5Application.core.file_io import OpenHDF5File
from KratosMultiphysics.DigitalTwinApplication.utilities.data_io import SupportedVariableUnionType
from KratosMultiphysics.DigitalTwinApplication.utilities.expression_utils import ExpressionUnionType

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
    if not isinstance(parameters, Kratos.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    if not isinstance(model, Kratos.Model):
        raise Exception("expected input shall be a model object")

    return VariableErrorOutputProcess(model, parameters["settings"], optimization_problem)

class ErrorComputation:
    def __init__(self, model_part: Kratos.ModelPart, parameters: Kratos.Parameters) -> None:
        self.model_part = model_part
        self.parameters = parameters

        default_parameters = Kratos.Parameters("""{
            "name"          : "",
            "container_type": "element_properties_value",
            "prefix"        : "/TestData/",
            "variable_names": ["YOUNG_MODULUS"]
        }""")

        self.parameters.ValidateAndAssignDefaults(default_parameters)
        self.name = self.parameters["name"].GetString()
        self.h5_params = Kratos.Parameters("""{
            "prefix": ""
        }""")
        self.h5_params["prefix"].SetString(self.parameters["prefix"].GetString())
        self.list_of_vars: 'list[SupportedVariableUnionType]' = [Kratos.KratosGlobals.GetVariable(var_name) for var_name in self.parameters["variable_names"].GetStringArray()]

        container_type = self.parameters["container_type"].GetString()
        self.expression_type: 'typing.Type[ExpressionUnionType]'
        self.read_function: 'typing.Callable[[ExpressionUnionType, SupportedVariableUnionType], None]'
        if container_type == "nodal_historical_value":
            self.expression_type = Kratos.Expression.NodalExpression(self.model_part)
            self.read_function = lambda exp, var: Kratos.Expression.VariableExpressionIO.Read(exp, var, True)
        elif container_type == "nodal_non_historical_value":
            self.expression_type = Kratos.Expression.NodalExpression
            self.read_function = lambda exp, var: Kratos.Expression.VariableExpressionIO.Read(exp, var, False)
        elif container_type == "condition_data_value":
            self.expression_type = Kratos.Expression.ConditionExpression
            self.read_function = lambda exp, var: Kratos.Expression.VariableExpressionIO.Read(exp, var)
        elif container_type == "condition_properties_value":
            self.expression_type = Kratos.Expression.ConditionExpression
            self.read_function = lambda exp, var: KratosOA.PropertiesVariableExpressionIO.Read(exp, var)
        elif container_type == "element_data_value":
            self.expression_type = Kratos.Expression.ElementExpression
            self.read_function = lambda exp, var: Kratos.Expression.VariableExpressionIO.Read(exp, var)
        elif container_type == "element_properties_value":
            self.expression_type = Kratos.Expression.ElementExpression
            self.read_function = lambda exp, var: KratosOA.PropertiesVariableExpressionIO.Read(exp, var)
        else:
            raise RuntimeError(f"Unsupported container_type = \"{container_type}\".")

        self.reference_values_list: 'list[ExpressionUnionType]' = []

    def ReadReferenceValues(self, h5_file: KratosHDF5.HDF5File) -> float:
        exp_io = KratosHDF5.ExpressionIO(self.h5_params, h5_file)
        for var in self.list_of_vars:
            exp = self.expression_type(self.model_part)
            exp_io.Read(var.Name(), exp)
            self.reference_values_list.append(exp)

    def ComputeError(self) -> float:
        resulting_error = 0.0
        for var, ref_exp in zip(self.list_of_vars, self.reference_values_list):
            current_value_expression = self.expression_type(self.model_part)
            self.read_function(current_value_expression, var)
            resulting_error += KratosOA.ExpressionUtils.NormL2(ref_exp - current_value_expression) ** 2
        return (resulting_error / len(self.list_of_vars)) ** 0.5

    def GetName(self) -> str:
        return self.name

class VariableErrorOutputProcess(Kratos.OutputProcess):
    def GetDefaultParameters(self) -> Kratos.Parameters:
        return Kratos.Parameters("""
        {
            "model_part_name"        : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "reference_file_settings": {
                "file_name": "Output/<model_part_name>.h5",
                "file_access_mode": "read_only"
            },
            "list_of_variable_settings": [
                {
                    "name"          : "",
                    "container_type": "element_properties_value",
                    "prefix"        : "/TestData/",
                    "variable_names": ["YOUNG_MODULUS"]
                }
            ]
        }""")

    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> None:
        Kratos.OutputProcess.__init__(self)
        self.parameters = parameters
        self.parameters.RecursivelyValidateAndAssignDefaults(self.GetDefaultParameters())
        self.model_part = model[self.parameters["model_part_name"].GetString()]
        self.optimization_problem = optimization_problem
        self.list_of_error_computations = [ErrorComputation(self.model_part, params) for params in self.parameters["list_of_variable_settings"].values()]

    def ExecuteInitialize(self) -> None:
        with OpenHDF5File(self.parameters["reference_file_settings"], self.model_part) as h5_file:
            for error_computation in self.list_of_error_computations:
                error_computation.ReadReferenceValues(h5_file)

    def PrintOutput(self) -> None:
        unbuffered_data = ComponentDataView("algorithm", self.optimization_problem).GetUnBufferedData()
        buffered_data = ComponentDataView("algorithm", self.optimization_problem).GetBufferedData()
        for error_computation in self.list_of_error_computations:
            current_error = error_computation.ComputeError()
            if not unbuffered_data.HasValue(error_computation.GetName()):
                unbuffered_data[error_computation.GetName()] = current_error

            buffered_data[f"{error_computation.GetName()}_abs"] = current_error
            buffered_data[f"{error_computation.GetName()}_rel"] = current_error / unbuffered_data[error_computation.GetName()]

    def IsOutputStep(self) -> bool:
        return True

