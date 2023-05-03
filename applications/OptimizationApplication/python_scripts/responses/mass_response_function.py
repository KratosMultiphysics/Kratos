import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.response_function import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import ConvertCollectiveExpressionValueMapToModelPartValueMap

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ResponseFunction:
    return MassResponseFunction(model, parameters, optimization_problem)

class MassResponseFunction(ResponseFunction):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        default_settings = Kratos.Parameters("""{
            "combined_output_model_part_name": "<RESPONSE_NAME>_combined_no_neighbours_response",
            "evaluated_model_part_names"     : [
                "PLEASE_PROVIDE_A_MODEL_PART_NAME"
            ]
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.output_model_part_name = parameters["combined_output_model_part_name"].GetString()
        self.model_part_names = parameters["evaluated_model_part_names"].GetStringArray()
        self.optimization_problem = optimization_problem
        self.model = model
        self.model_part: Kratos.ModelPart = None

        if len(self.model_part_names) == 0:
            raise RuntimeError("No model parts were provided for MassResponseFunction.")

    def GetDependentPhysicalKratosVariables(self) -> list[SupportedSensitivityFieldVariableTypes]:
        return [KratosOA.SHAPE, Kratos.DENSITY, Kratos.THICKNESS, KratosOA.CROSS_AREA]

    def Initialize(self) -> None:
        # get the model part name
        response_name = self.optimization_problem.GetComponentName(self)
        output_model_part_name = self.output_model_part_name.replace("<RESPONSE_NAME>", response_name)

        # get root model part
        model_parts_list = [self.model[model_part_name] for model_part_name in self.model_part_names]
        root_model_part = model_parts_list[0].GetRootModelPart()

        # create the combined model part
        if not root_model_part.HasSubModelPart(output_model_part_name):
            self.model_part = Kratos.ModelPartOperationUtilities.Merge(output_model_part_name, root_model_part, model_parts_list, False)
        else:
            self.model_part = root_model_part.GetSubModelPart(output_model_part_name)

    def Check(self) -> None:
        KratosOA.ResponseUtils.MassResponseUtils.Check(self.model_part)

    def Finalize(self) -> None:
        pass

    def GetEvaluatedModelPart(self) -> Kratos.ModelPart:
        if self.model_part is None:
            raise RuntimeError("Please call MassResponseFunction::Initialize first.")
        return self.model_part

    def GetAnalysisModelPart(self) -> None:
        return None

    def CalculateValue(self) -> float:
        return KratosOA.ResponseUtils.MassResponseUtils.CalculateValue(self.model_part)

    def CalculateGradient(self, physical_variable_collective_expressions: dict[SupportedSensitivityFieldVariableTypes, KratosOA.ContainerExpression.CollectiveExpressions]) -> None:
        # first calculate the gradients
        KratosOA.ResponseUtils.MassResponseUtils.CalculateGradient(self.model_part, ConvertCollectiveExpressionValueMapToModelPartValueMap(physical_variable_collective_expressions))

        # now fill the collective expressions
        for variable, collective_expression in physical_variable_collective_expressions.items():
            collective_expression.Read(Kratos.KratosGlobals.GetVariable(variable.Name() + "_SENSITIVITY"))

    def __str__(self) -> str:
        return f"Response [type = {self.__class__.__name__}, name = {self.optimization_problem.GetComponentName(self)}, model part name = {self.model_part.FullName()}]"