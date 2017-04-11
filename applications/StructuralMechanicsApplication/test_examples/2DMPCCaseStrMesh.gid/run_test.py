from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#### TIME MONITORING START ####

# time control starts
import time as timer
print(timer.ctime())
# measure process time
t0p = timer.clock()
# measure wall time
t0w = timer.time()

#
def StartTimeMeasuring():
    # measure process time
    time_ip = timer.clock()
    return time_ip

def StopTimeMeasuring(time_ip, process, report):
    # measure process time
    time_fp = timer.clock()
    if( report ):
        used_time = time_fp - time_ip
        print("::[KSM Simulation]:: [ %.2f" % round(used_time,2),"s", process," ] ")

#### TIME MONITORING END ####

#import kratos core and applications
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *
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

#add variables
solver.AddVariables()

#read model_part (note: the buffer_size is set here) (restart can be read here)
solver.ImportModelPart()

#add dofs 
solver.AddDofs()

##get the list of the submodel part in the object Model
for i in range(ProjectParameters["solver_settings"]["processes_sub_model_part_list"].size()):
    part_name = ProjectParameters["solver_settings"]["processes_sub_model_part_list"][i].GetString()
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
list_of_processes = process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["constraints_process_list"] )

list_of_processes += process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["loads_process_list"] )
            
#print list of constructed processes
if(echo_level>1):
    for process in list_of_processes:
        print(process)

#TODO: decide which is the correct place to initialize the processes 
for process in list_of_processes:
    process.ExecuteInitialize()

#### processes settings end ####

#### START SOLUTION ####

#### output settings start ####

problem_path = os.getcwd()
problem_name = ProjectParameters["problem_data"]["problem_name"].GetString()

# initialize GiD  I/O (gid outputs, file_lists)
from gid_output_process import GiDOutputProcess
output_settings = ProjectParameters["output_configuration"]
gid_output = GiDOutputProcess(main_model_part,
                              problem_name,
                              output_settings)

gid_output.ExecuteInitialize()
#### output settings end ####

## Sets strategies, builders, linear solvers, schemes and solving info, and fills the buffer
solver.Initialize()

print(" ")
print("::[KSM Simulation]:: Analysis -START- ")

for process in list_of_processes:
    process.ExecuteBeforeSolutionLoop()
    
## Set results when are written in a single file
gid_output.ExecuteBeforeSolutionLoop()

## Stepping and time settings (get from process info or solving info)
#delta time
delta_time = ProjectParameters["problem_data"]["time_step"].GetDouble()
#start step
step       = 0
#start time
time       = ProjectParameters["problem_data"]["start_time"].GetDouble()
#end time
end_time   = ProjectParameters["problem_data"]["end_time"].GetDouble()

main_model_part.Nodes[10].SetSolutionStepValue(DISPLACEMENT_X, 0, 0.01)    
main_model_part.Nodes[10].SetSolutionStepValue(DISPLACEMENT_Y, 0, 0.0) 
main_model_part.Nodes[10].Fix(DISPLACEMENT_X) 
main_model_part.Nodes[10].Fix(DISPLACEMENT_Y)


main_model_part.Nodes[11].SetSolutionStepValue(DISPLACEMENT_X, 0, 0.01)    
main_model_part.Nodes[11].SetSolutionStepValue(DISPLACEMENT_Y, 0, 0.1)
main_model_part.Nodes[11].Fix(DISPLACEMENT_X) 
main_model_part.Nodes[11].Fix(DISPLACEMENT_Y)


main_model_part.Nodes[12].SetSolutionStepValue(DISPLACEMENT_X, 0, 0.01)    
main_model_part.Nodes[12].SetSolutionStepValue(DISPLACEMENT_Y, 0, 0.01)
main_model_part.Nodes[12].Fix(DISPLACEMENT_X) 
main_model_part.Nodes[12].Fix(DISPLACEMENT_Y)

constraintMaker = AddMultiPointConstraint() 

# solving the problem (time integration)
while(time <= end_time):

    time = time + delta_time
    step = step + 1
    main_model_part.CloneTimeStep(time)

   
    # print process info
    ## 
    print("##########################")
    print("Time :: ",time)
    print("Step :: ",step)
    
    for process in list_of_processes:
        process.ExecuteInitializeSolutionStep()

    gid_output.ExecuteInitializeSolutionStep()
    
    print("########################## Solving START")
    solver.Solve()
    print("########################## Solving END")
    
    if(step == 1):
        main_model_part.Nodes[6].SetValue(IS_SLAVE, True)
        main_model_part.Nodes[7].SetValue(IS_SLAVE, True)
        main_model_part.Nodes[9].SetValue(IS_SLAVE, True)
        print("Setting MPC Constraints .. !!") 
       
        constraintMaker.ApplyConstraint(main_model_part.Nodes[16],DISPLACEMENT_Y, main_model_part.Nodes[6], DISPLACEMENT_Y,  1.0 )        
        constraintMaker.ApplyConstraint(main_model_part.Nodes[16],DISPLACEMENT_X, main_model_part.Nodes[6], DISPLACEMENT_X,  1.0 )   
        constraintMaker.ApplyConstraint(main_model_part.Nodes[17],DISPLACEMENT_X, main_model_part.Nodes[7], DISPLACEMENT_X,  1.0 )
        constraintMaker.ApplyConstraint(main_model_part.Nodes[17],DISPLACEMENT_Y, main_model_part.Nodes[7], DISPLACEMENT_Y,  1.0 )                   
        constraintMaker.ApplyConstraint(main_model_part.Nodes[18],DISPLACEMENT_X, main_model_part.Nodes[9], DISPLACEMENT_X,  1.0 )   
        constraintMaker.ApplyConstraint(main_model_part.Nodes[18],DISPLACEMENT_Y, main_model_part.Nodes[9], DISPLACEMENT_Y,  1.0 )   

        print("########################## Solving START")
        solver.Solve()
        print("########################## Solving END")
       
    for process in list_of_processes:
        process.ExecuteFinalizeSolutionStep()
    
    gid_output.ExecuteFinalizeSolutionStep()

    for process in list_of_processes:
        process.ExecuteBeforeOutputStep()
    
    if gid_output.IsOutputStep():
        gid_output.PrintOutput()
                      
    for process in list_of_processes:
        process.ExecuteAfterOutputStep()


for process in list_of_processes:
    process.ExecuteFinalize()
    
# ending the problem (time integration finished)
gid_output.ExecuteFinalize()

print("::[KSM Simulation]:: Analysis -END- ")
print(" ")

#### END SOLUTION ####

# measure process time
tfp = timer.clock()
# measure wall time
tfw = timer.time()

print("::[KSM Simulation]:: [ Computing Time = (%.2f" % (tfp - t0p)," seconds process time) ( %.2f" % (tfw - t0w)," seconds wall time) ]")

print(timer.ctime())






