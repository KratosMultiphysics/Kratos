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
    module = __import__(item["kratos_module"].GetString())
    interface_file = __import__(item["python_module"].GetString())
    p = interface_file.Factory(item, Model)
    list_of_processes.append( p )
    print("done ",i)
            

for process in list_of_processes:
    print("a")
    process.ExecuteInitialize()
    print("b")

#TODO: think if there is a better way to do this
fluid_model_part = solver.GetComputingModelPart()


# initialize GiD  I/O
from gid_output import GiDOutput
output_settings = ProjectParameters["output_configuration"]
gid_io = GiDOutput(output_settings["output_filename"].GetString(),
                   output_settings["volume_output"].GetBool(),
                   output_settings["gid_post_mode"].GetString(),
                   output_settings["gid_multi_file_flag"].GetString(),
                   output_settings["gid_write_mesh_flag"].GetBool(),
                   output_settings["gid_write_conditions_flag"].GetBool())
output_time = output_settings["output_time"].GetDouble()

gid_io.initialize_results(fluid_model_part)

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

    for process in list_of_processes:
        process.ExecuteInitializeSolutionStep()
        
    solver.Solve()
        
    for process in list_of_processes:
        process.ExecuteFinalizeSolutionStep()

    for process in list_of_processes:
        process.ExecuteBeforeOutputStep()
    
    if(output_time <= out):
        #TODO: following lines shall not be needed once the gid_io is adapted to using the parameters
        nodal_results = []
        for i in range(output_settings["nodal_results"].size()):
            nodal_results.append(output_settings["nodal_results"][i].GetString())
        gauss_points_results = []
        for i in range(output_settings["gauss_points_results"].size()):
            gauss_points_results.append(output_settings["gauss_points_results"][i].GetString())
            
        gid_io.write_results(
            time,
            fluid_model_part,
            nodal_results,
            gauss_points_results)
        out = 0
        
    for process in list_of_processes:
        process.ExecuteAfterOutputStep()

    out = out + Dt

gid_io.finalize_results()

for process in list_of_processes:
    process.ExecuteFinalize()