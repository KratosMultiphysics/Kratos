from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#import kratos core and applications
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.ShapeOptimizationApplication import *
from KratosMultiphysics.ALEApplication import *

# For time measures
import time as timer

# ======================================================================================================================================
# Kratos CSM part
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
# Import mesh-motion solver and add solution variables
import mesh_solver_structural_similarity as mesh_solver_class
mesh_solver_class.AddVariables(main_model_part)

# Create solver for all response functions specified in the optimization settings 
# Note that internally variables related to the individual functions are added to the model part
import optimization_settings as opt_settings
import response_function_factory
response_function_solver = response_function_factory.CreateSolver( main_model_part, opt_settings )

# --------------------------------------------------------------------------

#read model_part (note: the buffer_size is set here) (restart can be read here)
CSM_solver.ImportModelPart()

#add dofs (always after importing the model part) (it must be integrated in the ImportModelPart)
# if we integrate it in the model part we cannot use combined solvers
CSM_solver.AddDofs()

# --------------------------------------------------------------------------
# Add Dofs for mesh solver
mesh_solver_class.AddDofs(main_model_part)

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
# Set and initialize mesh-solver
reform_dofs_at_each_step = False
compute_reactions = True
mesh_solver = mesh_solver_class.MeshSolverStructuralSimilarity(main_model_part,reform_dofs_at_each_step,compute_reactions)
mesh_solver.Initialize()

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
# Mesh motion part
# ======================================================================================================================================

def MoveMesh(X, opt_itr):
    
    # Extract surface nodes
    sub_model_part_name = "surface_nodes"     
    GeometryUtilities(main_model_part).extract_surface_nodes(sub_model_part_name)

    # Apply shape update as boundary condition for computation of mesh displacement 
    for node in main_model_part.GetSubModelPart(sub_model_part_name).Nodes:
        if(node.Id in X.keys()):
            node.Fix(MESH_DISPLACEMENT_X)
            node.Fix(MESH_DISPLACEMENT_Y)
            node.Fix(MESH_DISPLACEMENT_Z)              
            disp = Vector(3)
            disp[0] = X[node.Id][0]
            disp[1] = X[node.Id][1]
            disp[2] = X[node.Id][2]
            node.SetSolutionStepValue(MESH_DISPLACEMENT,0,disp)

    # Solve for mesh-update
    mesh_solver.Solve()

    # Update reference mesh (Since shape updates are imposed as incremental quantities)
    mesh_solver.UpdateReferenceMesh()

# --------------------------------------------------------------------------
def ComputeAndAddMeshDerivatives(dFdXs, dFdX):
   
    # Here we solve the pseudo-elastic mesh-motion system again using modified BCs
    # The contributions from the mesh derivatives appear as reaction forces

    for node in main_model_part.Nodes:

        # Apply dirichlet conditions
        if node.Id in dFdXs.keys():
            node.Fix(MESH_DISPLACEMENT_X)
            node.Fix(MESH_DISPLACEMENT_Y)
            node.Fix(MESH_DISPLACEMENT_Z)              
            xs = Vector(3)
            xs[0] = 0.0
            xs[1] = 0.0
            xs[2] = 0.0
            node.SetSolutionStepValue(MESH_DISPLACEMENT,0,xs)
        # Apply RHS conditions       
        else:
            rhs = Vector(3)
            rhs[0] = dFdX[node.Id][0]
            rhs[1] = dFdX[node.Id][1]
            rhs[2] = dFdX[node.Id][2]
            node.SetSolutionStepValue(MESH_RHS,0,rhs)

    # Solve mesh-motion problem with previously modified BCs
    mesh_solver.Solve()

    # Compute and add gradient contribution from mesh motion
    for node_id in dFdXs.keys():
        node = main_model_part.Nodes[node_id]
        sens_contribution = Vector(3)
        sens_contribution = node.GetSolutionStepValue(MESH_REACTION)
        dFdXs[node.Id] = dFdXs[node_id] + sens_contribution

# ======================================================================================================================================
# Optimization part
# ======================================================================================================================================

# Import kratos shape tools
import optimizer_factory

# --------------------------------------------------------------------------
def analyzer(X, controls, opt_itr, response):

    # Compute primals
    print("\n> Starting calculation of response values")
    start_time = timer.time()

    if(controls["strain_energy"]["calc_value"]):

        # Advance time iterator of main_model_part (for proper computation of mesh & state)
        step = float(opt_itr)
        main_model_part.CloneTimeStep(step)
    
        # Move mesh
        print("\n> Start ALEApplication to move mesh")
        MoveMesh(X, opt_itr)

        interim_time = timer.time()
        print("> Time needed for moving the mesh = ",round(interim_time - start_time,2),"s")

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

        # If contribution from mesh-motion to gradient shall be considered
        # ComputeAndAddMeshDerivatives(dFdXs, dFdX) 

        # Store gradient in response container
        response["strain_energy"]["gradient"] = dFdXs

    stop_time = timer.time()
    print("> Time needed for calculating gradients = ",round(stop_time - start_time,2),"s")

# --------------------------------------------------------------------------
# Initalize model_part here to have it available for further use in this main script
design_surface = ModelPart("design_surface")

# Create an optimizer object 
optimizer = optimizer_factory.CreateOptimizer(design_surface,opt_settings,analyzer)

print("\n> ==============================================================================================================")
print("> Starting optimization")
print("> ==============================================================================================================\n")

# Call the solve function in the optimizer
optimizer.optimize()

print("\n> ==============================================================================================================")
print("> Finished shape optimization!")
print("> ==============================================================================================================\n")

# Finish output of CSM
FinalizeKSMProcess()

# ======================================================================================================================================