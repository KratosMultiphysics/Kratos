import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_wrapper import ExecutionPolicyWrapper
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction

from KratosMultiphysics.StructuralMechanicsApplication import structural_mechanics_analysis
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication


class DisplacementResponseFunction(ResponseFunction):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        default_settings = Kratos.Parameters("""{
            "evaluated_model_part_name": "PLEASE_PROVIDE_A_MODEL_PART_NAME",
            "primal_analysis_name"     : "",
            "perturbation_size"        : 1e-8
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)
        self.model_part = model[parameters["evaluated_model_part_name"].GetString()]
        self.primal_analysis_execution_policy_wrapper: ExecutionPolicyWrapper = optimization_info.GetOptimizationRoutine(ExecutionPolicyWrapper, parameters["primal_analysis_name"].GetString())
        self.perturbation_size = parameters["perturbation_size"].GetDouble()

        self.adjoint_model = Kratos.Model()
        adjoint_file_name = "sensitivity_parameters.json"
        with open(adjoint_file_name, 'r') as parameter_file:
            self.adjoint_parameters = Kratos.Parameters(parameter_file.read())

        self.adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(self.adjoint_model, self.adjoint_parameters)
        print(self.adjoint_model)
        print(self.adjoint_analysis)

        names = self.adjoint_model.GetModelPartNames()
        print(names)

    def Check(self):
        pass

    def CalculateValue(self) -> float:
        # execute the primal analysis
        self.primal_analysis_execution_policy_wrapper.Execute()

        nodes = self.model_part.GetNodes()
        response_value = 0
        for i, node in enumerate(nodes):
            response_value += abs(node.GetSolutionStepValue(Kratos.DISPLACEMENT_X))
            response_value += abs(node.GetSolutionStepValue(Kratos.DISPLACEMENT_Y))
            response_value += abs(node.GetSolutionStepValue(Kratos.DISPLACEMENT_Z))

        return response_value

    def CalculateSensitivity(self, sensitivity_variable: any, sensitivity_model_part: Kratos.ModelPart):

        if sensitivity_variable == KratosOA.YOUNG_MODULUS_SENSITIVITY:

            self.adjoint_analysis.Run()
            adjoint_model_part = self.adjoint_analysis.model.GetModelPart("Structure_sensitivity")

            counter = 0
            for element in adjoint_model_part.GetElements():
                print(element.GetValue(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT))
                print(element.GetValue(StructuralMechanicsApplication.YOUNG_MODULUS_SENSITIVITY))
                counter += 1

            print(self.adjoint_model)
            print(self.adjoint_analysis)
            print(f"Number of total elements: {counter}")

            # print(sensitivity_model_part)
            # print(self.adjoint_parameters["solver_settings"]["response_function_settings"])

            # adjoint_displacement_response_function = StructuralMechanicsApplication.AdjointNodalDisplacementResponseFunction(
            #     self.model_part, self.adjoint_parameters["solver_settings"]["response_function_settings"])
            # adjoint_displacement_response_function.

            # print(adjoint_displacement_response_function)

            return 0
            #     KratosOA.LinearStrainEnergyResponseUtilities.CalculateStrainEnergyShapeSensitivity(sensitivity_model_part, self.perturbation_size, sensitivity_variable)
            # elif sensitivity_variable in [KratosOA.DENSITY_SENSITIVITY, KratosOA.YOUNG_MODULUS_SENSITIVITY, KratosOA.POISSON_RATIO_SENSITIVITY]:
            #     primal_variable = Kratos.KratosGlobals.GetVariable(sensitivity_variable.Name()[:-12])
            #     KratosOA.LinearStrainEnergyResponseUtilities.CalculateStrainEnergyElementPropertiesSensitivity(sensitivity_model_part, self.perturbation_size, primal_variable, sensitivity_variable)
        else:
            msg = f"Unsupported sensitivity w.r.t. {sensitivity_variable.Name()} requested for {sensitivity_model_part.FullName()}."
            msg += "Followings are supported options:"
            msg += "\n\tYOUNG_MODULUS_SENSITIVITY"
            raise RuntimeError(msg)
