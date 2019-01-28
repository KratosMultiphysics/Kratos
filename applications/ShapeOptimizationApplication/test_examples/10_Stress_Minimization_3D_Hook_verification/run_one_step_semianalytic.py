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

# ==============================================================================
# Some preprocessing
with open(kratos_response_settings["primal_settings"].GetString(),'r') as primal_parameters:
    primal_params = Parameters(primal_parameters.read())
model_part_name = primal_params["solver_settings"]["model_part_name"].GetString()

# Compute reference values
model = Model()
response = structural_response_function_factory.CreateResponseFunction(kratos_response_settings["response_type"].GetString(), kratos_response_settings, model)
response.RunCalculation(True)

print("response_value = ", response.GetValue())