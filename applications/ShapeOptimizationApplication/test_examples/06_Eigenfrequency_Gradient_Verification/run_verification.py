# Import Kratos core and apps
from KratosMultiphysics import *
from KratosMultiphysics.StructuralMechanicsApplication import *
import structural_response_function_factory

# ==============================================================================
# Define parameters
kratos_response_settings = Parameters("""
{
            "response_type"          : "eigenfrequency",
            "primal_settings"        : "parameters.json",
            "gradient_mode"          : "semi_analytic",
            "step_size"              : 1e-8,
            "traced_eigenfrequencies": [1]
}""")

# Note that depending on move direction the initial pertubations might be too big leading to a crashing simulation
move_node_id = 28
move_direction = 'z'
levels_of_finite_differencing = 12

# ==============================================================================
# Some preprocessing
with open(kratos_response_settings["primal_settings"].GetString(),'r') as parameters:
    parameters = Parameters(parameters.read())
model_part_name = parameters["solver_settings"]["model_part_name"].GetString()

# Compute reference values
model = Model()
response = structural_response_function_factory.CreateResponseFunction(kratos_response_settings["response_type"].GetString(), kratos_response_settings, model)
response.RunCalculation(True)
response_values = [response.GetValue()]
response_gradient_reference = response.GetShapeGradient()

# Compute finite difference values
deltas = ['-']
for itr in range(0,levels_of_finite_differencing):

    # Compute pertubation
    current_delta = 1*10**(-itr)
    deltas.append(current_delta)

    model = Model()
    response = structural_response_function_factory.CreateResponseFunction(kratos_response_settings["response_type"].GetString(), kratos_response_settings, model)

    # Read mdpa
    response.Initialize()

    # Perturb
    model_part = model.GetModelPart(model_part_name)
    if move_direction == 'x':
        model_part.Nodes[move_node_id].X0 += current_delta
        model_part.Nodes[move_node_id].X += current_delta
    elif move_direction == 'y':
        model_part.Nodes[move_node_id].Y0 += current_delta
        model_part.Nodes[move_node_id].Y += current_delta
    elif move_direction == 'z':
        model_part.Nodes[move_node_id].Z0 += current_delta
        model_part.Nodes[move_node_id].Z += current_delta

    # Perform analysis
    response.InitializeSolutionStep()
    response.CalculateValue()
    response.CalculateGradient()
    response.FinalizeSolutionStep()
    response.Finalize()

    # Store value
    response_values.append(response.GetValue())

# ==============================================================================
# Evaluation of FD gradients
fd_gradients = ['-']
reference_value = response_values[0]
for itr in range(1,len(response_values)):
    fd_gradients.append((response_values[itr] - reference_value) / deltas[itr])

# Print results
print("----------------------------------------")
print("values:")
for itr in range(0,len(response_values)):
    print(response_values[itr])
print("----------------------------------------")
print("\nfd_gradients:")
for itr in range(0,len(fd_gradients)):
    print(deltas[itr],fd_gradients[itr])
print("----------------------------------------")
print("semi_analytic_gradient at move_node:")
print(response_gradient_reference[move_node_id])
print("----------------------------------------")

# ==============================================================================