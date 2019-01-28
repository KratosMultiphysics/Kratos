# Import Kratos core and apps
import KratosMultiphysics as Kratos
import KratosMultiphysics.StructuralMechanicsApplication as KratosCSM
import KratosMultiphysics.ExternalSolversApplication as KratosExternalSolvers
import structural_response_function_factory

# Additional imports
import os
from KratosMultiphysics.KratosUnittest import TestCase
from KratosMultiphysics import kratos_utilities

# ==============================================================================
# Define parameters
kratos_response_settings_1 = Kratos.Parameters("""
{
            "response_type"           : "strain_energy",
            "primal_settings"         : "primal_parameters.json",
            "gradient_mode"           : "semi_analytic",
            "step_size"               : 1e-10,
            "consider_discretization" : false
}""")

kratos_response_settings_2 = Kratos.Parameters("""
{
    "response_type"    : "adjoint_linear_strain_energy",
    "gradient_mode"    : "semi_analytic",
    "step_size"        : 1e-10,
    "primal_settings"  : "primal_parameters.json",
    "adjoint_settings" : "adjoint_parameters.json"
}""")

kratos_response_settings_3 = Kratos.Parameters("""
{
    "response_type"             : "adjoint_external_function",
    "gradient_mode"             : "semi_analytic",
    "step_size"                 : 1e-10,
    "primal_settings"           : "primal_parameters.json",
    "adjoint_settings"          : "adjoint_parameters.json"
}""")

# ==============================================================================
# Test original strain enery

model_1 = Kratos.Model()

response_1 = structural_response_function_factory.CreateResponseFunction(kratos_response_settings_1["response_type"].GetString(), kratos_response_settings_1, model_1)

primal_model_part = response_1.primal_model_part

response_1.Initialize()
response_1.InitializeSolutionStep()
response_1.CalculateValue()
response_1.CalculateGradient()

response_value_1 = response_1.GetValue()
response_gradient_1 = response_1.GetShapeGradient()

response_1.FinalizeSolutionStep()
response_1.Finalize()

# ==============================================================================
# Test adjoint strain enery

model_2 = Kratos.Model()

response_2 = structural_response_function_factory.CreateResponseFunction(kratos_response_settings_2["response_type"].GetString(), kratos_response_settings_2, model_2)

adjoint_model_part = response_2.adjoint_model_part

response_2.Initialize()
response_2.InitializeSolutionStep()
response_2.CalculateValue()
response_2.CalculateGradient()

response_value_2 = response_2.GetValue()
response_gradient_2 = response_2.GetShapeGradient()

response_2.FinalizeSolutionStep()
response_2.Finalize()

# ==============================================================================
# Test adjoint strain enery with external input

model_3 = Kratos.Model()

response_3 = structural_response_function_factory.CreateResponseFunction(kratos_response_settings_3["response_type"].GetString(), kratos_response_settings_3, model_3)

primal_model_part = response_3.primal_model_part
adjoint_model_part = response_3.adjoint_model_part

response_3.Initialize()

# External input
primal_model_part.ProcessInfo[KratosCSM.RESPONSE_VALUE] = 0.01098096046370471
for node in adjoint_model_part.Nodes:

    local_load_vector = primal_model_part.Nodes[node.Id].GetSolutionStepValue(KratosCSM.POINT_LOAD)
    local_load_vector = [0.5*value for value in local_load_vector]

    node.SetSolutionStepValue(KratosCSM.DFDU,local_load_vector)
    node.SetSolutionStepValue(KratosCSM.DFDX,[0,0,0])

response_3.InitializeSolutionStep()
response_3.CalculateValue()
response_3.CalculateGradient()

response_value_3 = response_3.GetValue()
response_gradient_3 = response_3.GetShapeGradient()

response_3.FinalizeSolutionStep()
response_3.Finalize()

# ==============================================================================
# Compare results

# Some output
print("-------------------------------------------------")
print("valueresponse_value_1 = ", response_value_1)
print("valueresponse_value_2 = ", response_value_2)
print("valueresponse_value_3 = ", response_value_3)
print("-------------------------------------------------")

# Test values
decimal_accuray_value = 10

print("> Starting to test value_2...")
TestCase().assertAlmostEqual(response_value_1, response_value_2, decimal_accuray_value)
print("> OK")

print("> Starting to test value_3...")
TestCase().assertAlmostEqual(response_value_1, response_value_3, decimal_accuray_value)
print("> OK")

# Test gradients
decimal_accuray_gradient = 10

print("> Starting to test gradient_2...")
for node_id, gradient in response_gradient_2.items():

    reference = response_gradient_1[node_id][0]
    TestCase().assertAlmostEqual(reference, gradient[0], decimal_accuray_gradient)

    reference = response_gradient_1[node_id][1]
    TestCase().assertAlmostEqual(reference, gradient[1], decimal_accuray_gradient)

    reference = response_gradient_1[node_id][2]
    TestCase().assertAlmostEqual(reference, gradient[2], decimal_accuray_gradient)

print("> OK")

print("> Starting to test gradient_3...")
for node_id, gradient in response_gradient_3.items():

    reference = response_gradient_1[node_id][0]
    TestCase().assertAlmostEqual(reference, gradient[0], decimal_accuray_gradient)

    reference = response_gradient_1[node_id][1]
    TestCase().assertAlmostEqual(reference, gradient[1], decimal_accuray_gradient)

    reference = response_gradient_1[node_id][2]
    TestCase().assertAlmostEqual(reference, gradient[2], decimal_accuray_gradient)

print("> OK")
print("-------------------------------------------------")

# Cleaning
original_directory = os.getcwd()
kratos_utilities.DeleteDirectoryIfExisting("__pycache__")
kratos_utilities.DeleteFileIfExisting(os.path.basename(original_directory)+".post.lst")
kratos_utilities.DeleteFileIfExisting("primal_results.post.bin")
kratos_utilities.DeleteFileIfExisting("adjoint_results.post.bin")
kratos_utilities.DeleteFileIfExisting("ProjectParametersOutput.json")

# ==============================================================================