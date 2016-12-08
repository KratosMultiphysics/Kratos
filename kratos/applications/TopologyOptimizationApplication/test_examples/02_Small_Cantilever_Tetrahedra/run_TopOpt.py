from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#### TIME MONITORING START ####

# time control starts
import time as timer
print(timer.ctime())

#### TIME MONITORING END ####

#import kratos core and applications
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.TopologyOptimizationApplication import *
from KratosMultiphysics.ExternalSolversApplication import *

######################################################################################
######################################################################################
######################################################################################

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
solver = solver_module.CreateSolver(main_model_part, ProjectParameters["solver_settings"])

#add variables (always before importing the model part) (it must be integrated in the ImportModelPart)
# if we integrate it in the model part we cannot use combined solvers
solver.AddVariables()

#read model_part (note: the buffer_size is set here) (restart can be read here)
solver.ImportModelPart()

#add dofs (always after importing the model part) (it must be integrated in the ImportModelPart)
# if we integrate it in the model part we cannot use combined solvers
solver.AddDofs()

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

#obtain the list of the processes to be applied

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

computing_model_part = solver.GetComputingModelPart()


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

# restart write included in gid IO ??

#### output settings end ####

## Sets strategies, builders, linear solvers, schemes and solving info, and fills the buffer
solver.Initialize()

# Start process
for process in list_of_processes:
    process.ExecuteBeforeSolutionLoop()
    
## Set results when are written in a single file
gid_output.ExecuteBeforeSolutionLoop()

def solve_structure(opt_itr):    

    # measure process time
    t0p = timer.clock()
    # measure wall time
    t0w = timer.time()

    # Iterate solver time
    time = opt_itr
    step = opt_itr
    main_model_part.CloneTimeStep(time)

     # processes to be executed at the begining of the solution step
    for process in list_of_processes:
        process.ExecuteInitializeSolutionStep()

    gid_output.ExecuteInitializeSolutionStep()
        
    # Actual solution
    solver.Solve()
       
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

    # measure process time
    tfp = timer.clock()
    # measure wall time
    tfw = timer.time()

    print("  KSM simulation finished                    [ Computing Time = (%.2f" % (tfp - t0p)," seconds process time) ( %.2f" % (tfw - t0w)," seconds wall time) ]")        

def FinalizeKSMProcess():

    for process in list_of_processes:
        process.ExecuteFinalize()
    
    # ending the problem (time integration finished)
    gid_output.ExecuteFinalize()

    # End simulation
    print("::[KSM]:: Analysis -END- ")
    print(" ")
    print(timer.ctime())

# ======================================================================================================================================
# Optimization part
# ======================================================================================================================================

def Analyzer(controls, response, opt_itr):
     
    # Create object to analyze structure response functions if required
    response_analyzer = StructureResponseFunctionUtilities(main_model_part)

    # Create object to analyze sensitivities if required (dummy linear solver for future adjoint sensitivity analysis)
    import linear_solver_factory
    linear_solver = linear_solver_factory.ConstructSolver(ProjectParameters["solver_settings"]["linear_solver_settings"])
    sensitivity_solver = StructureAdjointSensitivityStrategy(main_model_part, linear_solver,ProjectParameters["problem_data"]["domain_size"].GetInt())

    # Compute objective function value Call the Solid Mechanics Application to compute objective function value
    if(controls["strain_energy"]["calc_func"]):

        # Compute structure solution to get displacement field u
        print("::[START Kratos Solid Mechanics Application]::")
        solve_structure(opt_itr)

        # Calculate objective function value based on u and save into container
        print("\n::[Objective Function Calculation]::")
        response["strain_energy"]["func"] = response_analyzer.ComputeStrainEnergy()

    # Compute constraint function value
    if(controls["volume_fraction"]["calc_func"]):
        print("\n::[Constraint Function Calculation]::")
        target_volume_fraction = opt_parameters.initial_volume_fraction
        response["volume_fraction"]["func"] = response_analyzer.ComputeVolumeFraction() - target_volume_fraction

    # Compute sensitivities of objective function
    if(controls["strain_energy"]["calc_grad"]):        
        print("\n::[Objective sensitivity computation]::")
        sensitivity_solver.ComputeStrainEnergySensitivities()

    # Compute sensitivities of constraint function
    if(controls["volume_fraction"]["calc_grad"]):  
        print("\n::[Constraint sensitivity computation]::")
        sensitivity_solver.ComputeVolumeFractionSensitivities()

# Construct optimizer
import topology_optimizer_factory
import OptimizationParameters as opt_parameters
optimizer = topology_optimizer_factory.ConstructOptimizer(main_model_part, opt_parameters, Analyzer)

# Preprocess model (Set SOLID and VOID elements)
#for element_i in model_part.Elements:
#    center_coord_x  = element_i.GetValue(ELEM_CENTER_X)
#    center_coord_y  = element_i.GetValue(ELEM_CENTER_Y)
#    center_coord_z  = element_i.GetValue(ELEM_CENTER_Z)
#    if center_coord_y < 20 and center_coord_y > 19:
#      element_i.SetValue(SOLID_VOID, 1)

# Start optimization
optimizer.optimize()

# ======================================================================================================================================
# End of optimization part
# ======================================================================================================================================

FinalizeKSMProcess()
