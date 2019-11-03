# Import Kratos core and apps
import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as KCSM

# Additional imports
import KratosMultiphysics.kratos_utilities as kratos_utilities
from KratosMultiphysics.StructuralMechanicsApplication import structural_response_function_factory
import time, os
from decimal import Decimal

# ==============================================================================
# User defined parameters
# ==============================================================================
kratos_response_settings = KM.Parameters("""
{
    "response_type"     : "adjoint_local_stress",
    "gradient_mode"     : "semi_analytic",
    "step_size"         : 1.0e-8,
    "traced_element_id" : 1,
    "stress_type"       : "VON_MISES_STRESS",
    "stress_treatment"  : "mean",
    "primal_settings"   : "primal_parameters.json",
    "adjoint_settings"  : "auto",
    "sensitivity_settings" : {
        "sensitivity_model_part_name"     : "Parts_structure",
        "nodal_solution_step_sensitivity_variables"     : ["SHAPE_SENSITIVITY"],
        "element_data_value_sensitivity_variables"   : [],
        "condition_data_value_sensitivity_variables" : [],
        "build_mode": "static"
    }
}""")

list_of_move_nodes = [5, 65]
finite_difference_levels = [2,3,4,5,6,7,8,9,10,11,12]

results_filename = "fd_semianalytic_results.txt"
delete_results_afterwards = True

# ==============================================================================
# Preprocessing
# ==============================================================================
with open(kratos_response_settings["primal_settings"].GetString(),'r') as primal_parameters:
    primal_parameters = KM.Parameters(primal_parameters.read())
model_part_name = primal_parameters["solver_settings"]["model_part_name"].GetString()
domain_size = primal_parameters["solver_settings"]["domain_size"].GetInt()

open_file = open(results_filename,"w")
open_file.close()

# ==============================================================================
# Analysis
# ==============================================================================
for move_node_id in list_of_move_nodes:

    # Initialize results file
    with open(results_filename,'a') as open_file:
        open_file.write("-----------------------------------------------------------------------\n")
        open_file.write("time\t\t\t= "+str(time.ctime())+"\n")
        open_file.write("move_node_id\t= "+str(move_node_id)+"\n")
        open_file.write("---\n")
        open_file.write("step,\tgradient_x,\tgradient_y,\tgradient_z,\ttime\n")

    for itr in finite_difference_levels:

        print("\n\n##########################################################################")
        print(">> Starting finite difference level ",itr, "for node ", str(move_node_id))
        print("##########################################################################\n\n")

        # Perturb
        current_delta = 1*10**(-itr)
        kratos_response_settings["step_size"].SetDouble(current_delta)

        # Run simulation
        model = KM.Model()
        response = structural_response_function_factory.CreateResponseFunction(kratos_response_settings["response_type"].GetString(), kratos_response_settings, model)
        response.RunCalculation(True)

        # Write results
        gradient_of_interest = response.GetShapeGradient()[move_node_id]
        with open(results_filename,'a') as open_file:
            line_to_write = '%.0E' % Decimal(str(current_delta)) + ",\t"
            line_to_write += '%.6E' % Decimal(str(gradient_of_interest[0])) + ",\t"
            line_to_write += '%.6E' % Decimal(str(gradient_of_interest[1])) + ",\t"
            line_to_write += '%.6E' % Decimal(str(gradient_of_interest[2])) + ",\t"
            line_to_write += str(time.ctime()) + "\n"
            open_file.write(line_to_write)

# ==============================================================================
# Postprocessing
# ==============================================================================
if delete_results_afterwards:
    kratos_utilities.DeleteDirectoryIfExisting("__pycache__")
    kratos_utilities.DeleteFileIfExisting("sensitivity_verification_process_test.post.lst")
    kratos_utilities.DeleteFileIfExisting(results_filename)
    kratos_utilities.DeleteFileIfExisting("structure.post.bin")