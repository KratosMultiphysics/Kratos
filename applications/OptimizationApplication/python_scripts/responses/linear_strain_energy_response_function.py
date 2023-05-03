import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_decorator import ExecutionPolicyDecorator
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.response_function import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import ConvertCollectiveExpressionValueMapToModelPartValueMap

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ResponseFunction:
    return LinearStrainEnergyResponseFunction(model, parameters, optimization_problem)

class LinearStrainEnergyResponseFunction(ResponseFunction):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        default_settings = Kratos.Parameters("""{
            "combined_output_model_part_name": "<RESPONSE_NAME>_combined_no_neighbours_response",
            "primal_analysis_name"           : "",
            "perturbation_size"              : 1e-8,
            "evaluated_model_part_names"     : [
                "PLEASE_PROVIDE_A_MODEL_PART_NAME"
            ]
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.output_model_part_name = parameters["combined_output_model_part_name"].GetString()
        self.model_part_names = parameters["evaluated_model_part_names"].GetStringArray()
        self.perturbation_size = parameters["perturbation_size"].GetDouble()

        self.model = model
        self.optimization_problem = optimization_problem
        self.model_part: Kratos.ModelPart = None
        self.primal_analysis_execution_policy_decorator: ExecutionPolicyDecorator = self.optimization_problem.GetExecutionPolicy(parameters["primal_analysis_name"].GetString())

        if len(self.model_part_names) == 0:
            raise RuntimeError("No model parts were provided for LinearStrainEnergyResponseFunction.")

    def GetDependentPhysicalKratosVariables(self) -> list[SupportedSensitivityFieldVariableTypes]:
        return [KratosOA.SHAPE, Kratos.YOUNG_MODULUS, Kratos.THICKNESS, Kratos.POISSON_RATIO]

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
        pass

    def Finalize(self) -> None:
        pass

    def GetEvaluatedModelPart(self) -> Kratos.ModelPart:
        if self.model_part is None:
            raise RuntimeError("Please call LinearStrainEnergyResponseFunction::Initialize first.")
        return self.model_part

    def GetAnalysisModelPart(self) -> Kratos.ModelPart:
        return self.primal_analysis_execution_policy_decorator.GetExecutionPolicy().GetAnalysisModelPart()

    def CalculateValue(self) -> float:
        return KratosOA.ResponseUtils.LinearStrainEnergyResponseUtils.CalculateValue(self.model_part)

    def CalculateGradient(self, physical_variable_collective_expressions: dict[SupportedSensitivityFieldVariableTypes, KratosOA.ContainerExpression.CollectiveExpressions]) -> None:
        # first calculate the gradients
        KratosOA.ResponseUtils.LinearStrainEnergyResponseUtils.CalculateGradient(self.GetAnalysisModelPart(), self.GetEvaluatedModelPart(), ConvertCollectiveExpressionValueMapToModelPartValueMap(physical_variable_collective_expressions), self.perturbation_size)

        # now fill the collective expressions
        for variable, collective_expression in physical_variable_collective_expressions.items():
            collective_expression.Read(Kratos.KratosGlobals.GetVariable(variable.Name() + "_SENSITIVITY"))

    def __str__(self) -> str:
        return f"Response [type = {self.__class__.__name__}, name = {self.optimization_problem.GetComponentName(self)}, model part name = {self.model_part.FullName()}]"