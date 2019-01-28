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
    "traced_element_id" : 1887,
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

# Compute reference values
model = Model()
response = structural_response_function_factory.CreateResponseFunction(kratos_response_settings["response_type"].GetString(), kratos_response_settings, model)
response.RunCalculation(True)
response_value_reference = response.GetValue()
response_gradient_reference = response.GetShapeGradient()

# Compute finite difference values
deltas = []
response_values_x = []
response_values_y = []
response_values_z = []

for itr in range(0,levels_of_finite_differencing):

    # Compute pertubation
    current_delta = 1*10**(-itr)
    deltas.append(current_delta)

    # Delta in x ######################################

    model = Model()
    response = structural_response_function_factory.CreateResponseFunction(kratos_response_settings["response_type"].GetString(), kratos_response_settings, model)

    # Read mdpa
    response.Initialize()

    # Perturb
    model_part = model.GetModelPart(model_part_name)
    model_part.Nodes[move_node_id].X0 += current_delta
    model_part.Nodes[move_node_id].X += current_delta

    # Perform analysis
    response.InitializeSolutionStep()
    response.CalculateValue()
    response.FinalizeSolutionStep()
    response.Finalize()

    # Store value
    response_values_x.append(response.GetValue())

    # Delta in y ######################################

    model = Model()
    response = structural_response_function_factory.CreateResponseFunction(kratos_response_settings["response_type"].GetString(), kratos_response_settings, model)

    # Read mdpa
    response.Initialize()

    # Perturb
    model_part = model.GetModelPart(model_part_name)
    model_part.Nodes[move_node_id].Y0 += current_delta
    model_part.Nodes[move_node_id].Y += current_delta

    # Perform analysis
    response.InitializeSolutionStep()
    response.CalculateValue()
    response.FinalizeSolutionStep()
    response.Finalize()

    # Store value
    response_values_y.append(response.GetValue())

    # Delta in z ######################################

    model = Model()
    response = structural_response_function_factory.CreateResponseFunction(kratos_response_settings["response_type"].GetString(), kratos_response_settings, model)

    # Read mdpa
    response.Initialize()

    # Perturb
    model_part = model.GetModelPart(model_part_name)
    model_part.Nodes[move_node_id].Z0 += current_delta
    model_part.Nodes[move_node_id].Z += current_delta

    # Perform analysis
    response.InitializeSolutionStep()
    response.CalculateValue()
    response.FinalizeSolutionStep()
    response.Finalize()

    # Store value
    response_values_z.append(response.GetValue())

# ==============================================================================
# Evaluation of FD gradients
fd_gradients = []
for itr in range(0,len(response_values_x)):
    grad_x = (response_values_x[itr] - response_value_reference) / deltas[itr]
    grad_y = (response_values_y[itr] - response_value_reference) / deltas[itr]
    grad_z = (response_values_z[itr] - response_value_reference) / deltas[itr]
    fd_gradients.append([grad_x,grad_y,grad_z])

# Print results
print("\n----------------------------------------")
print("Reference value = ", response_value_reference)
print("\n----------------------------------------")
print("values after local pertubation: x,y,z")
for itr in range(0,len(response_values_x)):
    print(response_values_x[itr], response_values_x[itr], response_values_x[itr])
print("\n----------------------------------------")
print("fd_gradients:")
for itr in range(0,len(fd_gradients)):
    print(deltas[itr],fd_gradients[itr])
print("\n----------------------------------------")
print("semi_analytic_gradient at move_node:")
print(response_gradient_reference[move_node_id])
print("----------------------------------------")

# ==============================================================================