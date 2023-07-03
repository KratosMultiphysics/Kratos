import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_decorator import ExecutionPolicyDecorator
from KratosMultiphysics.OptimizationApplication.responses.measurement_residual_response_function import MeasurementResidualResponseFunction

model = Kratos.Model()
optimization_problem = OptimizationProblem()
model_part = model.CreateModelPart("Structure")
model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3

# create the primal analysis execution policy wrapper
with kratos_unittest.WorkFolderScope("measurement_residual_test", __file__):
    # creating the execution policy wrapper
    execution_policy_wrapper_settings = Kratos.Parameters("""{
            "name"    : "primal_analysis",
            "module": "KratosMultiphysics.OptimizationApplication.execution_policies",
            "type": "stepping_analysis_execution_policy",
            "settings": {
                "model_part_names" : ["Structure"],
                "analysis_module"  : "KratosMultiphysics.StructuralMechanicsApplication",
                "analysis_type"    : "StructuralMechanicsAnalysis",
                "analysis_settings": {
                    "@include_json": "primal_parameters.json"
                }
            },
            "pre_operations"           : [],
            "post_operations"          : [],
            "log_in_file"              : false,
            "log_file_name"            : "structure.log"
        }""")
    execution_policy_decorator = ExecutionPolicyDecorator(model, execution_policy_wrapper_settings, optimization_problem)
    optimization_problem.AddComponent(execution_policy_decorator)

    Kratos.ModelPartIO("model_file", Kratos.ModelPartIO.READ | Kratos.ModelPartIO.MESH_ONLY).ReadModelPart(model_part)

    # creating the response function wrapper
    response_function_settings = Kratos.Parameters("""{
            "evaluated_model_part_names": ["Structure"],
            "primal_analysis_name"      : "primal_analysis",
            "perturbation_size"         : 1e-8,
            "measurement_data_file": "MeasurementData_right_half.json"
        }""")
    response_function: MeasurementResidualResponseFunction = MeasurementResidualResponseFunction("measurement_residual", model, response_function_settings, optimization_problem)
    optimization_problem.AddComponent(response_function)

    execution_policy_decorator.Initialize()
    response_function.Initialize()

    # now replace the properties
    KratosOA.OptimizationUtils.CreateEntitySpecificPropertiesForContainer(model["Structure.all_nodes_elements_model_part"], model_part.Elements)

    for element in model_part.Elements:
        element.Properties[Kratos.YOUNG_MODULUS] += 30000000000.0 * (element.Id + 1)

    execution_policy_decorator.Execute()
    ref_value = response_function.CalculateValue()
    print(f"Objective value {ref_value}")

    sensitivity_expressions = KratosOA.CollectiveExpression([Kratos.Expression.ElementExpression(model_part), Kratos.Expression.ElementExpression(model_part)])
    response_function.CalculateGradient({Kratos.YOUNG_MODULUS: sensitivity_expressions})
    for expression in sensitivity_expressions.GetContainerExpressions():
        data = expression.Evaluate()
        # print(data)

    results = response_function.CalculateGradientWithFiniteDifferencing({Kratos.YOUNG_MODULUS: sensitivity_expressions})

    import numpy as np
    np_fd_data = np.array(results)
    # print(results)

    import matplotlib.pyplot as plt
    x = np.arange(0, len(model_part.Elements))
    plt.plot(x, -np_fd_data, ".", label="fd")
    plt.plot(x, data, "--", label="ad")
    plt.legend(loc="upper right")
    plt.show()
    print(np.linalg.norm(-np_fd_data - data))
