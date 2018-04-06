from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MeshingApplication import *
import KratosMultiphysics.MappingApplication as KratosMapping

# In this example two domains are solved, a coarse background mesh and a fine mesh around
# an obstacle. The fine domain receives the values from the coarse domain as input on it's boundary

######################################################################################
######################################################################################
######################################################################################
##PARSING THE PARAMETERS
#import define_output

parameter_file_background = open("ProjectParameters_Background.json",'r')
Projectparameters_BG = Parameters( parameter_file_background.read())

parameter_file_bodyfitted = open("ProjectParameters_BodyFitted.json",'r')
Projectparameters_BF = Parameters( parameter_file_bodyfitted.read())

## Fluid model part definition 
main_model_part_bg = ModelPart(Projectparameters_BG["problem_data"]["model_part_name"].GetString())
main_model_part_bg.ProcessInfo.SetValue(DOMAIN_SIZE, Projectparameters_BG["problem_data"]["domain_size"].GetInt())

main_model_part_bf = ModelPart(Projectparameters_BF["problem_data"]["model_part_name"].GetString())
main_model_part_bf.ProcessInfo.SetValue(DOMAIN_SIZE, Projectparameters_BF["problem_data"]["domain_size"].GetInt())

###TODO replace this "model" for real one once available
Model_BG = {Projectparameters_BG["problem_data"]["model_part_name"].GetString() : main_model_part_bg}
Model_BF = {Projectparameters_BF["problem_data"]["model_part_name"].GetString() : main_model_part_bf}

## Solver construction
solver_module = __import__(Projectparameters_BG["solver_settings"]["solver_type"].GetString())
solver_bg = solver_module.CreateSolver(main_model_part_bg, Projectparameters_BG["solver_settings"])

solver_bg.AddVariables()

solver_module = __import__(Projectparameters_BF["solver_settings"]["solver_type"].GetString())
solver_bf = solver_module.CreateSolver(main_model_part_bf, Projectparameters_BF["solver_settings"])

solver_bf.AddVariables()

## Read the model - note that SetBufferSize is done here
solver_bg.ImportModelPart()
solver_bf.ImportModelPart()

## Add AddDofs
solver_bg.AddDofs()
solver_bf.AddDofs()

## Initialize GiD  I/O
from gid_output_process import GiDOutputProcess
gid_output_bg = GiDOutputProcess(solver_bg.GetComputingModelPart(),
                              Projectparameters_BG["problem_data"]["problem_name"].GetString() ,
                              Projectparameters_BG["output_configuration"])

gid_output_bg.ExecuteInitialize()

gid_output_bf = GiDOutputProcess(solver_bf.GetComputingModelPart(),
                              Projectparameters_BF["problem_data"]["problem_name"].GetString() ,
                              Projectparameters_BF["output_configuration"])

gid_output_bf.ExecuteInitialize()

##here all of the allocation of the strategies etc is done
solver_bg.Initialize()
solver_bf.Initialize()

##TODO: replace MODEL for the Kratos one ASAP
## Get the list of the skin submodel parts in the object Model
for i in range(Projectparameters_BG["solver_settings"]["skin_parts"].size()):
    skin_part_name = Projectparameters_BG["solver_settings"]["skin_parts"][i].GetString()
    Model_BG.update({skin_part_name: main_model_part_bg.GetSubModelPart(skin_part_name)})
    
for i in range(Projectparameters_BF["solver_settings"]["skin_parts"].size()):
    skin_part_name = Projectparameters_BF["solver_settings"]["skin_parts"][i].GetString()
    Model_BF.update({skin_part_name: main_model_part_bf.GetSubModelPart(skin_part_name)})

## Get the list of the initial conditions submodel parts in the object Model
for i in range(Projectparameters_BF["initial_conditions_process_list"].size()):
    initial_cond_part_name = Projectparameters_BF["initial_conditions_process_list"][i]["Parameters"]["model_part_name"].GetString()
    Model_BF.update({initial_cond_part_name: main_model_part_bf.GetSubModelPart(initial_cond_part_name)})


## Processes construction    
import process_factory
# "list_of_processes_bg" contains all the processes already constructed (boundary conditions, initial conditions and gravity) 
# Note that the conditions are firstly constructed. Otherwise, they may overwrite the BCs information.
list_of_processes_bg = process_factory.KratosProcessFactory(Model_BG).ConstructListOfProcesses( Projectparameters_BG["initial_conditions_process_list"] )
list_of_processes_bg += process_factory.KratosProcessFactory(Model_BG).ConstructListOfProcesses( Projectparameters_BG["boundary_conditions_process_list"] )

list_of_processes_bf = process_factory.KratosProcessFactory(Model_BF).ConstructListOfProcesses( Projectparameters_BF["initial_conditions_process_list"] )
list_of_processes_bf += process_factory.KratosProcessFactory(Model_BF).ConstructListOfProcesses( Projectparameters_BF["boundary_conditions_process_list"] )

## Processes initialization
for process in list_of_processes_bg:
    process.ExecuteInitialize()
    
for process in list_of_processes_bf:
    process.ExecuteInitialize()

# Mapper initialization
mapper_settings_file = open("MapperSettings.json",'r')
Projectparameters_Mapper = Parameters( mapper_settings_file.read())["mapper_settings"]

inlet_mapper = KratosMapping.MapperFactory.CreateMapper(main_model_part_bg,
                                            main_model_part_bf,
                                            Projectparameters_Mapper[0])

sides_mapper = KratosMapping.MapperFactory.CreateMapper(main_model_part_bg,
                                            main_model_part_bf,
                                            Projectparameters_Mapper[1])

outlet_mapper = KratosMapping.MapperFactory.CreateMapper(main_model_part_bg,
                                            main_model_part_bf,
                                            Projectparameters_Mapper[2])

## Stepping and time settings
Dt = Projectparameters_BG["problem_data"]["time_step"].GetDouble()
end_time = Projectparameters_BG["problem_data"]["end_time"].GetDouble()

time = 0.0
step = 0
out = 0.0

gid_output_bg.ExecuteBeforeSolutionLoop()
gid_output_bf.ExecuteBeforeSolutionLoop()

for process in list_of_processes_bg:
    process.ExecuteBeforeSolutionLoop()
    
for process in list_of_processes_bf:
    process.ExecuteBeforeSolutionLoop()

while(time <= end_time):

    time = time + Dt
    step = step + 1
    main_model_part_bg.CloneTimeStep(time)
    main_model_part_bf.CloneTimeStep(time)

    print("STEP = ", step)
    print("TIME = ", time)

    if(step >= 3):
        for process in list_of_processes_bg:
            process.ExecuteInitializeSolutionStep()
        for process in list_of_processes_bf:
            process.ExecuteInitializeSolutionStep()
        
        gid_output_bg.ExecuteInitializeSolutionStep()
        gid_output_bf.ExecuteInitializeSolutionStep()
        
        solver_bg.Solve()
        inlet_mapper.Map(VELOCITY, VELOCITY)
        sides_mapper.Map(VELOCITY, VELOCITY)
        outlet_mapper.Map(VELOCITY, VELOCITY)
        solver_bf.Solve()
        
        for process in list_of_processes_bg:
            process.ExecuteFinalizeSolutionStep()
        for process in list_of_processes_bf:
            process.ExecuteFinalizeSolutionStep()
        
        gid_output_bg.ExecuteFinalizeSolutionStep()
        gid_output_bf.ExecuteFinalizeSolutionStep()

        #TODO: decide if it shall be done only when output is processed or not
        for process in list_of_processes_bg:
            process.ExecuteBeforeOutputStep()
        for process in list_of_processes_bf:
            process.ExecuteBeforeOutputStep()
        
    
        if gid_output_bg.IsOutputStep():
            gid_output_bg.PrintOutput()
            gid_output_bf.PrintOutput()
        
        for process in list_of_processes_bg:
            process.ExecuteAfterOutputStep()
        for process in list_of_processes_bf:
            process.ExecuteAfterOutputStep()
        

        out = out + Dt

for process in list_of_processes_bg:
    process.ExecuteFinalize()
for process in list_of_processes_bf:
    process.ExecuteFinalize()
    
gid_output_bg.ExecuteFinalize()
gid_output_bf.ExecuteFinalize()













