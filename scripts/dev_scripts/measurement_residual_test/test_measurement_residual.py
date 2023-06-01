import os
import sys

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_decorator import ExecutionPolicyDecorator
from KratosMultiphysics.OptimizationApplication.responses.measurement_likelihood_response_function import MeasurementLikelihoodResponseFunction

from KratosMultiphysics.StructuralMechanicsApplication import structural_mechanics_analysis
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# model = Kratos.Model()
# optimization_problem = OptimizationProblem()
# model_part = model.CreateModelPart("Structure")
# model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3

# # create the primal analysis execution policy wrapper
# with kratos_unittest.WorkFolderScope("linear_strain_energy_test", __file__):
#     # creating the execution policy wrapper
#     execution_policy_wrapper_settings = Kratos.Parameters("""{
#             "name"    : "primal",
#             "module": "KratosMultiphysics.OptimizationApplication.execution_policies",
#             "type": "stepping_analysis_execution_policy",
#             "settings": {
#                 "model_part_names" : ["Structure"],
#                 "analysis_module"  : "KratosMultiphysics.StructuralMechanicsApplication",
#                 "analysis_type"    : "StructuralMechanicsAnalysis",
#                 "analysis_settings": {
#                     "@include_json": "primal_parameters.json"
#                 }
#             },
#             "pre_operations"           : [],
#             "post_operations"          : [],
#             "log_in_file"              : false,
#             "log_file_name"            : "structure.log"
#         }""")
#     execution_policy_decorator = ExecutionPolicyDecorator(model, execution_policy_wrapper_settings, optimization_problem)
#     optimization_problem.AddComponent(execution_policy_decorator)

#     Kratos.ModelPartIO("Structure", Kratos.ModelPartIO.READ | Kratos.ModelPartIO.MESH_ONLY).ReadModelPart(model_part)

#     # creating the response function wrapper
#     response_function_settings = Kratos.Parameters("""{
#             "evaluated_model_part_names": ["Structure"],
#             "primal_analysis_name"      : "primal",
#             "perturbation_size"         : 1e-8
#         }""")
#     response_function: MeasurementLikelihoodResponseFunction = MeasurementLikelihoodResponseFunction("measurement_residual", model, response_function_settings, optimization_problem)
#     optimization_problem.AddComponent(response_function)

#     execution_policy_decorator.Initialize()
#     response_function.Initialize()

#     # now replace the properties
#     KratosOA.OptimizationUtils.CreateEntitySpecificPropertiesForContainer(model["Structure.whole_structure"], model_part.Elements)

#     execution_policy_decorator.Execute()
#     ref_value = response_function.CalculateValue()
#     print(f"Objective value {ref_value}")

#     sensitivity_expression = KratosOA.ContainerExpression.CollectiveExpressions([KratosOA.ContainerExpression.ElementPropertiesExpression(model_part)])
#     response_function.CalculateGradient({Kratos.YOUNG_MODULUS: sensitivity_expression})
#     for element in model_part.Elements:
#         print(f"Sensitivity: {element.Properties[KratosOA.YOUNG_MODULUS_SENSITIVITY]}")


primal_file_name = "/media/meister/localdata/GitRepos/Kratos_In_Progress_SysId/scripts/dev_scripts/measurement_residual_test/linear_shell_test_parameters.json"
adjoint_file_name = "/media/meister/localdata/GitRepos/Kratos_In_Progress_SysId/scripts/dev_scripts/measurement_residual_test/linear_shell_test_nodal_disp_adjoint_parameters.json"

with open(primal_file_name, 'r') as parameter_file:
    primal_parameters = Kratos.Parameters(parameter_file.read())
with open(adjoint_file_name, 'r') as parameter_file:
    adjoint_parameters = Kratos.Parameters(parameter_file.read())
problem_name = primal_parameters["problem_data"]["problem_name"].GetString()
model_part_name = primal_parameters["solver_settings"]["model_part_name"].GetString()

# To avoid many prints
if (primal_parameters["problem_data"]["echo_level"].GetInt() == 0 or adjoint_parameters["problem_data"]["echo_level"].GetInt() == 0):
    Kratos.Logger.GetDefaultOutput().SetSeverity(Kratos.Logger.Severity.WARNING)

# SelectAndVerifyLinearSolver(primal_parameters, skipTest)
# SelectAndVerifyLinearSolver(adjoint_parameters, skipTest)

# solve primal problem
model_primal = Kratos.Model()
primal_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model_primal, primal_parameters)
primal_analysis.Run()
# create adjoint analysis
model_adjoint = Kratos.Model()
adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model_adjoint, adjoint_parameters)
adjoint_analysis.Initialize()
adjoint_analysis.RunSolutionLoop()
adjoint_analysis.Finalize()

# for node in model_adjoint.GetModelPart("Structure").GetNodes():
#     print(f"Displacement: {node.GetSolutionStepValue(Kratos.DISPLACEMENT)}")
#     print(f"Displacement X: {node.GetSolutionStepValue(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT_X)}")
#     print(f"Displacement Y: {node.GetSolutionStepValue(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT_Y)}")
#     print(f"Displacement Z: {node.GetSolutionStepValue(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT_Z)}")

for element in model_adjoint.GetModelPart("Structure").GetElements():
    # print(f"Sensi. Thickness: {element.GetValue(StructuralMechanicsApplication.THICKNESS_SENSITIVITY)}")
    print(f"Sensi. Youngs Modulus: {element.GetValue(StructuralMechanicsApplication.YOUNG_MODULUS_SENSITIVITY)}")
