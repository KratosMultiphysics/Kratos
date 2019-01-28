# Import Kratos core and apps
from KratosMultiphysics import *
from KratosMultiphysics.StructuralMechanicsApplication import *
from KratosMultiphysics.ShapeOptimizationApplication import *
import structural_response_function_factory

# ==============================================================================
# Define parameters
kratos_response_settings = Parameters("""
{
    "response_type"     : "adjoint_local_stress",
    "gradient_mode"     : "semi_analytic",
    "step_size"         : 1.0e-8,
    "traced_element_id" : 22147,
    "stress_type"       : "VON_MISES_STRESS",
    "stress_treatment"  : "mean",
    "primal_settings"   : "primal_parameters.json",
    "adjoint_settings"  : "adjoint_parameters.json"
}""")

# Note that depending on move direction the initial pertubations might be too big leading to a crashing simulation
move_node_id = 15691
levels_of_finite_differencing = 12

# ==============================================================================
# Some preprocessing
with open(kratos_response_settings["primal_settings"].GetString(),'r') as primal_parameters:
    primal_params = Parameters(primal_parameters.read())
model_part_name = primal_params["solver_settings"]["model_part_name"].GetString()

# Compute finite difference values
deltas = []
response_values = []
response_gradients = []

for itr in range(2,levels_of_finite_differencing):

    # Compute pertubation
    current_delta = 1*10**(-itr)
    deltas.append(current_delta)
    kratos_response_settings["step_size"].SetDouble(current_delta)

    # Compute reference values
    model = Model()
    response = structural_response_function_factory.CreateResponseFunction(kratos_response_settings["response_type"].GetString(), kratos_response_settings, model)
    response.RunCalculation(True)

    # Store values
    response_values.append(response.GetValue())
    response_gradients.append(response.GetShapeGradient()[move_node_id])

    print(response_values)

    err


# ==============================================================================

# Print results
print("----------------------------------------")
print("values:")
for itr in range(0,len(response_values)):
    print(response_values[itr])
print("----------------------------------------")
print("\nresponse_gradients:")
for itr in range(0,len(response_gradients)):
    print(deltas[itr],response_gradients[itr])
print("----------------------------------------")

# ==============================================================================