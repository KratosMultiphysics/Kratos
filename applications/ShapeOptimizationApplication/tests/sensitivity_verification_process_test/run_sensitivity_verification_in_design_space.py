# Import Kratos core and apps
import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as KCSM
import KratosMultiphysics.ShapeOptimizationApplication as KSO

# Additional imports
import KratosMultiphysics.kratos_utilities as kratos_utilities
from KratosMultiphysics.StructuralMechanicsApplication import structural_response_function_factory
from KratosMultiphysics.ShapeOptimizationApplication.custom_variable_utilities import WriteDictionaryDataOnNodalVariable
from KratosMultiphysics.gid_output_process import GiDOutputProcess
import time
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
        "nodal_sensitivity_variables"     : ["SHAPE_SENSITIVITY"],
        "element_sensitivity_variables"   : [],
        "condition_sensitivity_variables" : [],
        "build_mode": "static"
    }
}""")

mapper_settings = KM.Parameters("""
{
    "filter_function_type"       : "linear",
    "filter_radius"              : 0.61,
    "max_nodes_in_filter_radius" : 1000,
    "matrix_free_filtering"      : false
}""")

list_of_move_nodes = [5, 65]
finite_difference_levels = [2,3,4,5,6,7,8,9,10,11,12]

results_filename = "fd_design_space_results.txt"
delete_results_afterwards = True

# ==============================================================================
# Helper function
# ==============================================================================
def OutputModelPart( output_mdpa, output_filename, nodal_results ):
    output_parameters = KM.Parameters(""" { "result_file_configuration" : { "nodal_results" : [] } }""")
    for entry in nodal_results:
        output_parameters["result_file_configuration"]["nodal_results"].Append(entry)
    gid_output_original = GiDOutputProcess( output_mdpa, output_filename, output_parameters )
    gid_output_original.ExecuteInitialize()
    gid_output_original.ExecuteBeforeSolutionLoop()
    gid_output_original.ExecuteInitializeSolutionStep()
    gid_output_original.PrintOutput()
    gid_output_original.ExecuteFinalizeSolutionStep()
    gid_output_original.ExecuteFinalize()

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

    print("\n\n##########################################################################")
    print(">> Starting to compute reference of node "+str(move_node_id))
    print("##########################################################################\n\n")

    # Run analysis
    model = KM.Model()
    response = structural_response_function_factory.CreateResponseFunction(kratos_response_settings["response_type"].GetString(), kratos_response_settings, model)
    model_part = model.GetModelPart(model_part_name)
    model_part.AddNodalSolutionStepVariable(KSO.DF1DX)
    model_part.AddNodalSolutionStepVariable(KSO.DF1DX_MAPPED)
    response.Initialize()
    response.InitializeSolutionStep()
    response.CalculateValue()
    response.CalculateGradient()
    response.FinalizeSolutionStep()
    response.Finalize()

    # Perform mapping
    WriteDictionaryDataOnNodalVariable(response.GetShapeGradient(), model_part, KSO.DF1DX)
    vm_mapper = KSO.MapperVertexMorphingMatrixFree(model_part, model_part, mapper_settings)
    vm_mapper.InverseMap(KSO.DF1DX, KSO.DF1DX_MAPPED)

    # Output to verify mapping
    OutputModelPart( model_part, "reference_part", ["DF1DX","DF1DX_MAPPED"])

    # Get results
    reference_value = response.GetValue()
    reference_gradient = model_part.Nodes[move_node_id].GetSolutionStepValue(KSO.DF1DX_MAPPED)

    # Initialize results file
    with open(results_filename,'a') as open_file:
        open_file.write("-----------------------------------------------------------------------\n")
        open_file.write("time\t\t\t\t= "+str(time.ctime())+"\n")
        open_file.write("move_node_id\t\t= "+str(move_node_id)+"\n")
        open_file.write("reference_value\t\t= "+str(reference_value)+"\n")
        open_file.write("reference_gradient\t= "+str(reference_gradient)+"\n")
        open_file.write("---\n")
        open_file.write("step,\tfd_gradient_x,\tfd_gradient_y,\tfd_gradient_z,\ttime\n")

    for itr in finite_difference_levels:
        current_delta = 1*10**(-itr)

        for dim_itr in range(domain_size):

            print("\n\n##########################################################################")
            print(">> Starting finite difference level ",itr, " dimension ",dim_itr+1, "for node ", str(move_node_id))
            print("##########################################################################\n\n")

            # Initialize new analysis
            model = KM.Model()
            response = structural_response_function_factory.CreateResponseFunction(kratos_response_settings["response_type"].GetString(), kratos_response_settings, model)
            model_part = model.GetModelPart(model_part_name)
            model_part.AddNodalSolutionStepVariable(KSO.CONTROL_POINT_UPDATE)
            model_part.AddNodalSolutionStepVariable(KSO.SHAPE_UPDATE)
            response.Initialize()

            # Perturb
            if dim_itr == 0:
                model_part.Nodes[move_node_id].SetSolutionStepValue(KSO.CONTROL_POINT_UPDATE, [current_delta,0,0])
            elif dim_itr == 1:
                model_part.Nodes[move_node_id].SetSolutionStepValue(KSO.CONTROL_POINT_UPDATE, [0,current_delta,0])
            elif dim_itr == 2:
                model_part.Nodes[move_node_id].SetSolutionStepValue(KSO.CONTROL_POINT_UPDATE, [0,0,current_delta])

            vm_mapper = KSO.MapperVertexMorphingMatrixFree(model_part, model_part, mapper_settings)
            vm_mapper.Map(KSO.CONTROL_POINT_UPDATE, KSO.SHAPE_UPDATE)

            KSO.MeshControllerUtilities(model_part).UpdateMeshAccordingInputVariable(KSO.SHAPE_UPDATE)
            KSO.MeshControllerUtilities(model_part).SetReferenceMeshToMesh()

            # Output to verify mapping
            OutputModelPart( model_part, "perturbed_part", ["CONTROL_POINT_UPDATE","SHAPE_UPDATE"])

            # Perform analysis
            response.InitializeSolutionStep()
            response.CalculateValue()
            response.FinalizeSolutionStep()
            response.Finalize()
            current_value = response.GetValue()

            # Evaluate fd gradients
            if dim_itr == 0:
                fd_gradient_x = (current_value-reference_value)/current_delta
            elif dim_itr == 1:
                fd_gradient_y = (current_value-reference_value)/current_delta
            elif dim_itr == 2:
                fd_gradient_z = (current_value-reference_value)/current_delta

        # Write results
        with open(results_filename,'a') as open_file:
            line_to_write = '%.0E' % Decimal(str(current_delta)) + ",\t"
            line_to_write += '%.6E' % Decimal(str(fd_gradient_x)) + ",\t"
            line_to_write += '%.6E' % Decimal(str(fd_gradient_y)) + ",\t"
            line_to_write += '%.6E' % Decimal(str(fd_gradient_z)) + ",\t"
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
    kratos_utilities.DeleteFileIfExisting("reference_part.post.bin")
    kratos_utilities.DeleteFileIfExisting("perturbed_part.post.bin")
