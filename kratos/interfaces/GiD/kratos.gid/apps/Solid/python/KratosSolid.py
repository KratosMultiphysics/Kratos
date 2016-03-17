from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#### TIME MONITORING START ####

# time control starts
from time import *
print(ctime())
# measure process time
t0p = clock()
# measure wall time
t0w = time()

#
def StartTimeMeasuring():
    # measure process time
    time_ip = clock()
    return time_ip

def StopTimeMeasuring(time_ip, process, report):
    # measure process time
    time_fp = clock()
    if( report ):
        used_time = time_fp - time_ip
        print("::[KSM Simulation]:: [ %.2f" % round(used_time,2),"s", process," ] ")

#### TIME MONITORING END ####

#import kratos core and applications
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *

parallel = OpenMPUtils()
print("::[KSM Simulation]:: OMP USING",parallel.GetNumThreads(),"THREADS ]")

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

#add variables (always before importing the model part)
solver.AddVariables()

#read model_part (note: the buffer_size is set here)
solver.ImportModelPart()

#add dofs (always after importing the model part)
solver.AddDofs()

#build sub_model_parts or submeshes (rearrange parts for the application of custom processes)
##TODO: replace MODEL for the Kratos one ASAP
##get the list of the submodel part in the object Model
for i in range(ProjectParameters["solver_settings"]["skin_parts"].size()):
    part_name = ProjectParameters["solver_settings"]["skin_parts"][i].GetString()
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
list_of_processes = process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["boundary_conditions_process_list"] )

#list_of_processes = []
#process_definition = ProjectParameters["boundary_conditions_process_list"]
#for i in range(process_definition.size()):
#    item = process_definition[i]
#    module = __import__(item["implemented_in_module"].GetString())
#    interface_file = __import__(item["implemented_in_file"].GetString())
#    p = interface_file.Factory(item, Model)
#    list_of_processes.append( p )
#    print("done ",i)
            
#print list of constructed processes
if(echo_level>1):
    for process in constructed_processes:
        print(process)

#TODO: decide which is the correct place to initialize the processes 
for process in list_of_processes:
    process.ExecuteInitialize()

#### processes settings end ####


#### output settings start ####

# initialize GiD  I/O (gid outputs, file_lists, output_frequency)
from gid_output import GiDOutput
output_settings = ProjectParameters["output_configuration"]
gid_io = GiDOutput(output_settings["output_filename"].GetString(),
                   output_settings["volume_output"].GetBool(),
                   output_settings["gid_post_mode"].GetString(),
                   output_settings["gid_multi_file_flag"].GetString(),
                   output_settings["gid_write_mesh_flag"].GetBool(),
                   output_settings["gid_write_conditions_flag"].GetBool())
output_time = output_settings["output_time"].GetDouble()

# restart write included in gid IO ¿?

#### output settings start ####


#### START SOLUTION ####

#TODO: think if there is a better way to do this
computing_model_part = solver.GetComputeModelPart()


## Sets strategies, builders, linear solvers, schemes and solving info, and fills the buffer
solver.Initialize()

## Set results when are written in a single file
gid_io.initialize_results(computing_model_part)

for process in list_of_processes:
    process.ExecuteBeforeSolutionLoop()

## Stepping and time settings (get from process info or solving info)
#delta time
delta_time = computing_model_part.ProcessInfo[DELTA_TIME]
#start step
step       = computing_model_part.ProcessInfo[TIME_STEPS] 
#start time
time       = computing_model_part.ProcessInfo[TIME]
#end time
end_time   = ProjectParameters["problem_data"]["end_time"].GetDouble()


# writing a initial state results file (if no restart)
gid_io.write_results(time, computing_model_part)

# solving the problem (time integration)
while(time <= end_time):

    #TODO: this must be done by a solving_info utility in the solver
    # store previous time step
    main_model_part.ProcessInfo[PREVIOUS_DELTA_TIME] = delta_time
    # set new time step ( it can change when solve is called )
    delta_time = model_part.ProcessInfo[DELTA_TIME]

    time = time + delta_time
    step = step + 1
    main_model_part.CloneTimeStep(time)

    # print process info
    ##
    
    for process in list_of_processes:
        process.ExecuteInitializeSolutionStep()
        
    solver.Solve()
        
    for process in list_of_processes:
        process.ExecuteFinalizeSolutionStep()

    #TODO: decide if it shall be done only when output is processed or not (boundary_conditions_processes ¿?)
    for process in list_of_processes:
        process.ExecuteBeforeOutputStep()
    
    # write results and restart files: (frequency writing is controlled internally by the gid_io)
    gid_io.write_results(time, computing_model_part)
              
    #TODO: decide if it shall be done only when output is processed or not
    for process in list_of_processes:
        process.ExecuteAfterOutputStep()


# ending the problem (time integration finished)
gid_io.finalize_results()

for process in list_of_processes:
    process.ExecuteFinalize()

print("::[KSM Simulation]:: Analysis -END- ")
print(" ")

# check solving information for any problem
solver.info_check()

#### END SOLUTION ####

# measure process time
tfp = clock()
# measure wall time
tfw = time()

print("::[KSM Simulation]:: [ Computing Time = (%.2f" % (tfp - t0p)," seconds process time) ( %.2f" % (tfw - t0w)," seconds wall time) ]")

print(ctime())

# to create a benchmark: add standard benchmark files and decomment next two lines 
# rename the file to: run_test.py
#from run_test_benchmark_results import *
#WriteBenchmarkResults(model_part)






