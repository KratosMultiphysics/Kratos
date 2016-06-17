from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
#~ from KratosMultiphysics.IncompressibleFluidApplication import *
#~ from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.FSIApplication import *
#~ from KratosMultiphysics.ExternalSolversApplication import *
#~ from KratosMultiphysics.MeshingApplication import *

######################################################################################
######################################################################################
######################################################################################
##PARSING THE PARAMETERS
#import define_output

parameter_file = open("ProjectParameters.json",'r')
ProjectParameters = Parameters( parameter_file.read())


## Fluid-Structure model parts definition 
structure_main_model_part = ModelPart(ProjectParameters["structure_solver_settings"]["problem_data"]["model_part_name"].GetString())
structure_main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, ProjectParameters["structure_solver_settings"]["problem_data"]["domain_size"].GetInt())

fluid_main_model_part = ModelPart(ProjectParameters["fluid_solver_settings"]["problem_data"]["model_part_name"].GetString())
fluid_main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, ProjectParameters["fluid_solver_settings"]["problem_data"]["domain_size"].GetInt())

###TODO replace this "model" for real one once available
FluidModel = {ProjectParameters["structure_solver_settings"]["problem_data"]["model_part_name"].GetString() : structure_main_model_part}
SolidModel = {ProjectParameters["fluid_solver_settings"]["problem_data"]["model_part_name"].GetString() : fluid_main_model_part}

## Solver construction
solver_module = __import__("partitioned_fsi_solver") # Currently there is only one FSI solver up to date
solver = solver_module.CreateSolver(structure_main_model_part, fluid_main_model_part, ProjectParameters)

solver.AddVariables()

## Read the model - note that SetBufferSize is done here
solver.ImportModelPart()

## Add AddDofs
solver.AddDofs()

## Initialize GiD  I/O
from gid_output_process import GiDOutputProcess

gid_output_structure = GiDOutputProcess(solver.structure_solver.GetComputeModelPart(),
                                    ProjectParameters["structure_solver_settings"]["problem_data"]["problem_name"].GetString()+"_structure",
                                    ProjectParameters["structure_solver_settings"]["output_configuration"])

gid_output_fluid = GiDOutputProcess(solver.fluid_solver.GetComputeModelPart(),
                                    ProjectParameters["fluid_solver_settings"]["problem_data"]["problem_name"].GetString()+"_fluid",
                                    ProjectParameters["fluid_solver_settings"]["output_configuration"])

gid_output_structure.ExecuteInitialize()
gid_output_fluid.ExecuteInitialize()

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
fluid_model_part = solver.GetComputeModelPart()

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













