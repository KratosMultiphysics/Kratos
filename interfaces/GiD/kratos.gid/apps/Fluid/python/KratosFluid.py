from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MeshingApplication import *

######################################################################################
######################################################################################
######################################################################################
##PARSING THE PARAMETERS
#import define_output

parameter_file = open("ProjectParameters.json",'r')
ProjectParameters = Parameters( parameter_file.read())


## defining a model part for the fluid
main_model_part = ModelPart(ProjectParameters["problem_data"]["model_part_name"].GetString())
main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, ProjectParameters["problem_data"]["domain_size"].GetInt())

###TODO replace this "model" for real one once available
Model = {ProjectParameters["problem_data"]["model_part_name"].GetString() : main_model_part}

#construct the solver
solver_module = __import__(ProjectParameters["solver_settings"]["solver_type"].GetString())
solver = solver_module.CreateSolver(main_model_part, ProjectParameters["solver_settings"])

solver.AddVariables()

#obtain the list of the processes to be applied
process_definition = Parameters(ProjectParameters["boundary_conditions_process_list"])

###read the model - note that SetBufferSize is done here
solver.ImportModelPart()

###add AddDofs
solver.AddDofs()

# initialize GiD  I/O
from gid_output_process import GiDOutputProcess
gid_output = GiDOutputProcess(solver.GetComputeModelPart(),
                              ProjectParameters["problem_data"]["problem_name"].GetString() ,
                              ProjectParameters["output_configuration"])

gid_output.ExecuteInitialize()

##here all of the allocation of the strategies etc is done
solver.Initialize()


##TODO: replace MODEL for the Kratos one ASAP
##get the list of the submodel part in the object Model
for i in range(ProjectParameters["solver_settings"]["skin_parts"].size()):
    part_name = ProjectParameters["solver_settings"]["skin_parts"][i].GetString()
    Model.update({part_name: main_model_part.GetSubModelPart(part_name)})

#TODO: decide which is the correct place to initialize the processes 
list_of_processes = []
process_definition = ProjectParameters["boundary_conditions_process_list"]
for i in range(process_definition.size()):
    item = process_definition[i]
    module = __import__(item["implemented_in_module"].GetString())
    interface_file = __import__(item["implemented_in_file"].GetString())
    p = interface_file.Factory(item, Model)
    list_of_processes.append( p )
    print("done ",i)
            

for process in list_of_processes:
    print("a")
    process.ExecuteInitialize()
    print("b")

#TODO: think if there is a better way to do this
fluid_model_part = solver.GetComputeModelPart()



gid_output.ExecuteBeforeSolutionLoop()
for process in list_of_processes:
    process.ExecuteBeforeSolutionLoop()

## Stepping and time settings
Dt = ProjectParameters["problem_data"]["time_step"].GetDouble()
end_time = ProjectParameters["problem_data"]["end_time"].GetDouble()


time = 0.0
step = 0
out = 0.0
while(time <= end_time):

    time = time + Dt
    step = step + 1
    main_model_part.CloneTimeStep(time)

    print("STEP = ", step)
    print("TIME = ", time)

    if(step >= 3):
        for process in list_of_processes:
            process.ExecuteInitializeSolutionStep()
        
        gid_output.ExecuteInitializeSolutionStep()
        
        solver.Solve()
        
        gid_output.ExecuteFinalizeSolutionStep()
        
        for process in list_of_processes:
            process.ExecuteFinalizeSolutionStep()

        #TODO: decide if it shall be done only when output is processed or not
        for process in list_of_processes:
            process.ExecuteBeforeOutputStep()
    
        if gid_output.IsOutputStep():
            gid_output.PrintOutput()
        
        for process in list_of_processes:
            process.ExecuteAfterOutputStep()

        out = out + Dt

gid_output.ExecuteFinalize()

for process in list_of_processes:
    process.ExecuteFinalize()













