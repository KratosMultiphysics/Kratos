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
        print("::[IGA Simulation]:: [ %.2f" % round(used_time,2),"s", process," ] ")

#### TIME MONITORING END ####

# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *
from KratosMultiphysics.IGAStructuralMechanicsApplication import *
import KratosMultiphysics.ExternalSolversApplication 

import json
import iga_structural_mechanics_solver

#region Set Up Model Part / Parsing the Parameters
#### PARSING THE PARAMETERS ####

#import define_output
parameter_file = open("ProjectParameters.json",'r')
ProjectParameters = Parameters( parameter_file.read())

#set echo level
echo_level = ProjectParameters["solver_settings"]["echo_level"].GetInt()

model_part = ModelPart(ProjectParameters["problem_data"]["model_part_name"].GetString())
# Not yet implemented in the parser
model_part.ProcessInfo.SetValue(DOMAIN_SIZE, ProjectParameters["problem_data"]["domain_size"].GetInt())


###TODO replace this "model" for real one once available in kratos core
Model = {ProjectParameters["problem_data"]["model_part_name"].GetString() : model_part}


#construct the solver (main setting methods are located in the solver_module)
solver_module = __import__(ProjectParameters["solver_settings"]["solver_type"].GetString())
solver = solver_module.CreateSolver(model_part, ProjectParameters["solver_settings"])


import read_materials_process
read_materials_process.Factory(ProjectParameters,Model)


solver.AddVariables()
solver.ImportModelPart(ProjectParameters)
solver.AddDofs()

#build sub_model_parts or submeshes (rearrange parts for the application of custom processes)
##TODO: replace MODEL for the Kratos one ASAP
##get the list of the submodel part in the object Model
for i in range(ProjectParameters["solver_settings"]["processes_sub_model_part_list"].size()):
    part_name = ProjectParameters["solver_settings"]["processes_sub_model_part_list"][i].GetString()
    Model.update({part_name: model_part.GetSubModelPart(part_name)})



#print model_part and properties
if(echo_level>0):
    print("")
    print(model_part)
    for properties in model_part.Properties:
        print(properties)
#### model_part settings end ####
#endregion
#region Process Settings
#--- processes settings start ####

import process_factory
list_of_processes = process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["constraints_process_list"] )
list_of_processes += process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["loads_process_list"] )

#region hard coded boundary conditions
#for node in model_part.Nodes:
	##node.Fix(DISPLACEMENT_X)
	#node.SetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER_X, 0, 0)
	#node.SetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER_Y, 0, 0)
	#node.SetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER_Z, 0, 0)
	#node.SetSolutionStepValue(DISPLACEMENT_Z, 0, 0)
	#node.Fix(DISPLACEMENT_Z)
	#node.SetSolutionStepValue(DISPLACEMENT_Y, 0, 0)
	#node.SetSolutionStepValue(DISPLACEMENT_Z, 0, 0)
	#node.Fix(DISPLACEMENT_Y)
	#node.Fix(DISPLACEMENT_Z)
	#if (node.X < 0.01):
	#	node.SetSolutionStepValue(DISPLACEMENT_X, 0, 0)
	#	node.SetSolutionStepValue(DISPLACEMENT_Y, 0, 0)
	#	node.SetSolutionStepValue(DISPLACEMENT_Z, 0,0)
	#	node.Fix(DISPLACEMENT_X)
	#	node.Fix(DISPLACEMENT_Y)
	#	node.Fix(DISPLACEMENT_Z)
	#if (node.X > 9.9):
	#	node.SetSolutionStepValue(KratosMultiphysics.SolidMechanicsApplication.POINT_LOAD_X, 0, 3.3333)
#endregion

##print list of constructed processes
if(echo_level>1):
    for process in list_of_processes:
        print(process)

##TODO: decide which is the correct place to initialize the processes 
for process in list_of_processes:
    process.ExecuteInitialize()

#### processes settings end ####
#endregion
#region Start Solution / Initialize / GiD I/O
#### START SOLUTION ####

#TODO: think if there is a better way to do this
computing_model_part = solver.GetComputeModelPart()


#### output settings start ####

problem_path = os.getcwd()
problem_name = ProjectParameters["problem_data"]["problem_name"].GetString()

# initialize GiD  I/O (gid outputs, file_lists)
import gid_output_process
from gid_output_process import GiDOutputProcess
output_settings = ProjectParameters["output_configuration"]
gid_output = GiDOutputProcess(computing_model_part, ProjectParameters["problem_data"]["problem_name"].GetString(), ProjectParameters["output_configuration"])


gid_output.ExecuteInitialize()

# restart write included in gid IO ??

#### output settings end ####


solver.Initialize()

for process in list_of_processes:
	process.ExecuteBeforeSolutionLoop()

## Set results when are written in a single file
gid_output.ExecuteBeforeSolutionLoop()

    
# write results and restart files: (frequency writing is controlled internally by the gid_io)
if gid_output.IsOutputStep():
    gid_output.PrintOutput()

#endregion
#region Time Integration
## Stepping and time settings (get from process info or solving info)
#delta time
delta_time = ProjectParameters["problem_data"]["time_step"].GetDouble()
#start step
step = 0
#start time
time = ProjectParameters["problem_data"]["start_time"].GetDouble()
#end time
end_time = ProjectParameters["problem_data"]["end_time"].GetDouble()



# solving the problem (time integration)
while(time <= end_time):

    #TODO: this must be done by a solving_info utility in the solver
    # store previous time step
    #~ computing_model_part.ProcessInfo[PREVIOUS_DELTA_TIME] = delta_time
    # set new time step ( it can change when solve is called )
    #~ delta_time = computing_model_part.ProcessInfo[DELTA_TIME]

    time = time + delta_time
    step = step + 1
    model_part.ProcessInfo[TIME_STEPS] = step
    model_part.CloneTimeStep(time)

    # print process info
    ##
    
    for process in list_of_processes:
        process.ExecuteInitializeSolutionStep()

    gid_output.ExecuteInitializeSolutionStep()
        
    solver.Solve()
       
    for process in list_of_processes:
        process.ExecuteFinalizeSolutionStep()
    
    gid_output.ExecuteFinalizeSolutionStep()

    #TODO: decide if it shall be done only when output is processed or not (boundary_conditions_processes ??)
    for process in list_of_processes:
        process.ExecuteBeforeOutputStep()
    
    # write results and restart files: (frequency writing is controlled internally by the gid_io)
    if gid_output.IsOutputStep():
        gid_output.PrintOutput()
                      
    #TODO: decide if it shall be done only when output is processed or not
    for process in list_of_processes:
        process.ExecuteAfterOutputStep()


for process in list_of_processes:
	process.ExecuteFinalize()

# ending the problem (time integration finished)
gid_output.ExecuteFinalize()
print("::[IGA Simulation]:: Analysis -END- ")
print(" ")

# check solving information for any problem
#~ solver.InfoCheck() # InfoCheck not implemented yet.
#endregion


#### END SOLUTION ####

# measure process time
tfp = timer.clock()
# measure wall time
tfw = timer.time()

print("::[IGA Simulation]:: [ Computing Time = (%.2f" % (tfp - t0p)," seconds process time) ( %.2f" % (tfw - t0w)," seconds wall time) ]")

print(timer.ctime())

for node in model_part.Nodes:
	print(node.GetSolutionStepValue(DISPLACEMENT_X))
	print(node.GetSolutionStepValue(DISPLACEMENT_Y))
	print(node.GetSolutionStepValue(DISPLACEMENT_Z))
