# Import Kratos core and apps
import os, shutil
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication
import structural_response_function_factory
import KratosMultiphysics.kratos_utilities as kratos_utils
import KratosMultiphysics.ExternalSolversApplication as ExternalSolversApplication
import structural_mechanics_analysis
import time as timer


with open("adjoint_strain_energy_response_parameters_nonlinear_structure.json",'r') as parameter_file:
    parameters = KratosMultiphysics.Parameters( parameter_file.read())

model = KratosMultiphysics.Model()
response_function = structural_response_function_factory.CreateResponseFunction("dummy", parameters["kratos_response_settings"], model)

model_part_primal = response_function.primal_model_part
model_part_adjoint = response_function.adjoint_model_part

response_function.RunCalculation(calculate_gradient=True)

## run primal analysis only without adjoint analysis
# with open("rectangular_plate_parameters_nonlinear.json",'r') as parameter_file:
#     parameters = KratosMultiphysics.Parameters( parameter_file.read())

# model = KratosMultiphysics.Model()

# model_part_primal = model.CreateModelPart("Structure" , 2)
# primal_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model, parameters)
# primal_analysis.Initialize()

# startTime = timer.time()
# if not primal_analysis.time < primal_analysis.end_time:
#     primal_analysis.end_time += 1

# primal_analysis.RunSolutionLoop()


##  trial to load serialized data using restart utilities
# import restart_utility
# import save_restart_process as save_rest_proc
# current_model = KratosMultiphysics.Model()
# ## hard coded
# load_model = current_model.CreateModelPart("rectangular_plate_structure")
# restart_parameters = KratosMultiphysics.Parameters("""
# {
#     "input_filename"                 : "test_restart_file",
#     "restart_load_file_label"        : "",
#     "serializer_trace"               : "trace_error",
#     "load_restart_files_from_folder" : false
# }
# """)
# restart_parameters["restart_load_file_label"].SetString( str(round(0.1, 3) ))
# rest_utility = restart_utility.RestartUtility(load_model, restart_parameters)
# rest_utility.LoadRestart()


model_part_name = "rectangular_plate_structure_primal"
for name in os.listdir():
    if name.find(model_part_name) == 0:
        kratos_utils.DeleteFileIfExisting(name)