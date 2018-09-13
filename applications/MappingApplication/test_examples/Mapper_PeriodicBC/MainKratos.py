from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MeshingApplication import *
import KratosMultiphysics.MappingApplication as KratosMapping

# "Abusing"" the Mapper to enforce periodic boundary conditions
# This only works because he channel is straight and nodes an be projected
# between Inlet and Outlet
# The whole domain has to be initialized with an initial velocity
# Mapping Velocity from Outlet to Inlet

######################################################################################
######################################################################################
######################################################################################
##PARSING THE PARAMETERS
#import define_output

parameter_file = open("ProjectParameters.json",'r')
ProjectParameters = Parameters( parameter_file.read())


## Fluid model part definition 
main_model_part = ModelPart(ProjectParameters["problem_data"]["model_part_name"].GetString())
main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, ProjectParameters["problem_data"]["domain_size"].GetInt())

###TODO replace this "model" for real one once available
Model = {ProjectParameters["problem_data"]["model_part_name"].GetString() : main_model_part}

## Solver construction
solver_module = __import__(ProjectParameters["solver_settings"]["solver_type"].GetString())
solver = solver_module.CreateSolver(main_model_part, ProjectParameters["solver_settings"])

solver.AddVariables()

## Read the model - note that SetBufferSize is done here
solver.ImportModelPart()

## Add AddDofs
solver.AddDofs()

## Initialize GiD  I/O
from gid_output_process import GiDOutputProcess
gid_output = GiDOutputProcess(solver.GetComputingModelPart(),
                              ProjectParameters["problem_data"]["problem_name"].GetString() ,
                              ProjectParameters["output_configuration"])

gid_output.ExecuteInitialize()

##here all of the allocation of the strategies etc is done
solver.Initialize()


##TODO: replace MODEL for the Kratos one ASAP
## Get the list of the skin submodel parts in the object Model
for i in range(ProjectParameters["solver_settings"]["skin_parts"].size()):
    skin_part_name = ProjectParameters["solver_settings"]["skin_parts"][i].GetString()
    Model.update({skin_part_name: main_model_part.GetSubModelPart(skin_part_name)})

## Get the list of the initial conditions submodel parts in the object Model
for i in range(ProjectParameters["initial_conditions_process_list"].size()):
    initial_cond_part_name = ProjectParameters["initial_conditions_process_list"][i]["Parameters"]["model_part_name"].GetString()
    Model.update({initial_cond_part_name: main_model_part.GetSubModelPart(initial_cond_part_name)})
    
## Get the gravity submodel part in the object Model
for i in range(ProjectParameters["gravity"].size()):   
    gravity_part_name = ProjectParameters["gravity"][i]["Parameters"]["model_part_name"].GetString()
    Model.update({gravity_part_name: main_model_part.GetSubModelPart(gravity_part_name)})


## Processes construction    
import process_factory
# "list_of_processes" contains all the processes already constructed (boundary conditions, initial conditions and gravity) 
# Note that the conditions are firstly constructed. Otherwise, they may overwrite the BCs information.
list_of_processes = process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["initial_conditions_process_list"] )
list_of_processes += process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["boundary_conditions_process_list"] )
list_of_processes += process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["gravity"] )

## Processes initialization
for process in list_of_processes:
    process.ExecuteInitialize()

#TODO: think if there is a better way to do this
fluid_model_part = solver.GetComputingModelPart()

mapper_settings_file = open("MapperSettings.json",'r')
Projectparameters_Mapper = Parameters( mapper_settings_file.read())["mapper_settings"][0]

periodic_bc_mapper = KratosMapping.MapperFactory.CreateMapper(main_model_part,
                                                 main_model_part,
                                                 Projectparameters_Mapper)

## Stepping and time settings
Dt = ProjectParameters["problem_data"]["time_step"].GetDouble()
end_time = ProjectParameters["problem_data"]["end_time"].GetDouble()

time = 0.0
step = 0
out = 0.0

gid_output.ExecuteBeforeSolutionLoop()

for process in list_of_processes:
    process.ExecuteBeforeSolutionLoop()

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
        periodic_bc_mapper.Map(VELOCITY,VELOCITY) # Map Velocity from Outlet to Inlet
        solver.Solve()
        
        for process in list_of_processes:
            process.ExecuteFinalizeSolutionStep()
            
        gid_output.ExecuteFinalizeSolutionStep()

        #TODO: decide if it shall be done only when output is processed or not
        for process in list_of_processes:
            process.ExecuteBeforeOutputStep()
    
        if gid_output.IsOutputStep():
            gid_output.PrintOutput()
        
        for process in list_of_processes:
            process.ExecuteAfterOutputStep()

        out = out + Dt

for process in list_of_processes:
    process.ExecuteFinalize()
    
gid_output.ExecuteFinalize()













