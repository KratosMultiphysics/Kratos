from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#import kratos core and applications
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# For time measures
import time as timer

# ======================================================================================================================================
# Solver preparation
# ======================================================================================================================================

#### PARSING THE PARAMETERS ####

#import define_output
parameter_file = open("ProjectParameters.json",'r')
ProjectParameters = Parameters( parameter_file.read())

#set echo level
echo_level = ProjectParameters["problem_data"]["echo_level"].GetInt()

#### model_part settings start ####

#defining the model_part
main_model_part = ModelPart(ProjectParameters["problem_data"]["model_part_name"].GetString())
main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, ProjectParameters["problem_data"]["domain_size"].GetInt())

###TODO replace this "model" for real one once available in kratos core
Model = {ProjectParameters["problem_data"]["model_part_name"].GetString() : main_model_part}

#construct the solver (main setting methods are located in the solver_module)
solver_module = __import__(ProjectParameters["solver_settings"]["solver_type"].GetString())
CSM_solver = solver_module.CreateSolver(main_model_part, ProjectParameters["solver_settings"])

#add variables (always before importing the model part) (it must be integrated in the ImportModelPart)
# if we integrate it in the model part we cannot use combined solvers
CSM_solver.AddVariables()

# --------------------------------------------------------------------------
import optimization_settings

# # Initalize model_part here to have it available for further use in this main script
# main_model_part = ModelPart(optimization_settings.model_input_filename)

# Create an optimizer 
# Note that internally variables related to the optimizer are added to the model part
import optimizer_factory
optimizer = optimizer_factory.CreateOptimizer( main_model_part, optimization_settings )

# Create solver for all response functions specified in the optimization settings 
# Note that internally variables related to the individual functions are added to the model part
import response_function_factory
response_function_solver = response_function_factory.CreateSolver( main_model_part, optimization_settings )

# --------------------------------------------------------------------------
#read model_part (note: the buffer_size is set here) (restart can be read here)
CSM_solver.ImportModelPart()

#add dofs (always after importing the model part) (it must be integrated in the ImportModelPart)
# if we integrate it in the model part we cannot use combined solvers
CSM_solver.AddDofs()

# --------------------------------------------------------------------------

# Build sub_model_parts or submeshes (rearrange parts for the application of custom processes)
## Get the list of the submodel part in the object Model
for i in range(ProjectParameters["solver_settings"]["processes_sub_model_part_list"].size()):
    part_name = ProjectParameters["solver_settings"]["processes_sub_model_part_list"][i].GetString()
    if( main_model_part.HasSubModelPart(part_name) ):
        Model.update({part_name: main_model_part.GetSubModelPart(part_name)})

#print model_part and properties
if(echo_level>1):
    print("")
    print(main_model_part)
    for properties in main_model_part.Properties:
        print(properties)

#### model_part settings end ####


#### processes settings start ####

import process_factory
#the process order of execution is important
list_of_processes  = process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["constraints_process_list"] )
list_of_processes += process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["loads_process_list"] )
if(ProjectParameters.Has("problem_process_list")):
    list_of_processes += process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["problem_process_list"] )
if(ProjectParameters.Has("output_process_list")):
    list_of_processes += process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["output_process_list"] )
            
#print list of constructed processes
if(echo_level>1):
    for process in list_of_processes:
        print(process)

for process in list_of_processes:
    process.ExecuteInitialize()

#### processes settings end ####

#### START SOLUTION ####

computing_model_part = CSM_solver.GetComputingModelPart()

#### output settings start ####

problem_path = os.getcwd()
problem_name = ProjectParameters["problem_data"]["problem_name"].GetString()

# initialize GiD  I/O (gid outputs, file_lists)
from gid_output_process import GiDOutputProcess
output_settings = ProjectParameters["output_configuration"]
gid_output = GiDOutputProcess(computing_model_part,
                              problem_name,
                              output_settings)

gid_output.ExecuteInitialize()

#### output settings end ####

## Sets strategies, builders, linear solvers, schemes and solving info, and fills the buffer
CSM_solver.Initialize()
CSM_solver.SetEchoLevel(echo_level)


# --------------------------------------------------------------------------
# Initialize response function solvers
for response_id in response_function_solver:
    response_function_solver[response_id].initialize()

# --------------------------------------------------------------------------

# Start process
for process in list_of_processes:
    process.ExecuteBeforeSolutionLoop()
## Set results when are written in a single file
gid_output.ExecuteBeforeSolutionLoop()

# --------------------------------------------------------------------------
# Call function to solve structure for the given state
def solve_structure(opt_itr): 

    print("Solving structure for step, ",opt_itr)

    # processes to be executed at the begining of the solution step
    for process in list_of_processes:
        process.ExecuteInitializeSolutionStep()

    gid_output.ExecuteInitializeSolutionStep()
        
    # Actual solution
    CSM_solver.Solve()
       
    for process in list_of_processes:
        process.ExecuteFinalizeSolutionStep()
    
    gid_output.ExecuteFinalizeSolutionStep()

    # processes to be executed at the end of the solution step
    for process in list_of_processes:
        process.ExecuteFinalizeSolutionStep()

    # processes to be executed before witting the output
    for process in list_of_processes:
        process.ExecuteBeforeOutputStep()
    
    # write output results GiD: (frequency writing is controlled internally)
    if(gid_output.IsOutputStep()):
        gid_output.PrintOutput()
                      
    # processes to be executed after witting the output
    for process in list_of_processes:
        process.ExecuteAfterOutputStep()            

# --------------------------------------------------------------------------
def FinalizeKSMProcess():
    for process in list_of_processes:
        process.ExecuteFinalize()
    gid_output.ExecuteFinalize()

# ======================================================================================================================================
# Optimization part
# ======================================================================================================================================

# --------------------------------------------------------------------------
def analyzeDesignAndRespondToOptimizer(X, controls, opt_itr, response):

    # Compute primals
    print("\n> Starting calculation of response values")
    start_time = timer.time()

    if(controls["strain_energy"]["calc_value"]):

        # Advance time iterator of main_model_part
        step = float(opt_itr)
        main_model_part.CloneTimeStep(step)

        # Udate mesh coordinates of main_model_part
        for node_id in X.keys():
            node = main_model_part.Nodes[node_id]
            node.X0 = node.X0 + X[node_id][0]
            node.Y0 = node.Y0 + X[node_id][1]
            node.Z0 = node.Z0 + X[node_id][2]

        interim_time = timer.time()
        print("> Time needed for updating the mesh = ",round(interim_time - start_time,2),"s")

        # Solve structural problem
        print("\n> Start SolidMechanicsApplication to solve structure")
        solve_structure(opt_itr)

        stop_time = timer.time()
        print("> Time needed for calculating structural response = ",round(stop_time - interim_time,2),"s")        

        # Calculate objective function value of current design and store in response container
        response_function_solver["strain_energy"].calculate_value()
        response["strain_energy"]["value"] = response_function_solver["strain_energy"].get_value()

    # Compute gradients
    print("\n> Start calculation of gradients")
    start_time = timer.time()

    if(controls["strain_energy"]["calc_gradient"]):        

        response_function_solver["strain_energy"].calculate_gradient()
        dFdX = response_function_solver["strain_energy"].get_gradient()

        # Get gradient information on design surface
        dFdXs = {}
        for node_id in X.keys():
            dFdXs[node_id] = dFdX[node_id]

        # Store gradient in response container
        response["strain_energy"]["gradient"] = dFdXs

    stop_time = timer.time()
    print("> Time needed for calculating gradients = ",round(stop_time - start_time,2),"s")

# --------------------------------------------------------------------------
optimizer.getAnalyzer().setAnalyzeFunction( analyzeDesignAndRespondToOptimizer )
optimizer.readInputModelPart()
optimizer.prepareOptimization()

print("\n> ==============================================================================================================")
print("> Starting optimization")
print("> ==============================================================================================================\n")

optimizer.optimize()

print("\n> ==============================================================================================================")
print("> Finished shape optimization!")
print("> ==============================================================================================================\n")

FinalizeKSMProcess()

# ======================================================================================================================================