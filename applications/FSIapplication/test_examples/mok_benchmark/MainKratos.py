from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
from KratosMultiphysics.ALEApplication import *
from KratosMultiphysics.FSIApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.SolidMechanicsApplication import *

import numpy

######################################################################################
######################################################################################
######################################################################################

## Parsing ProjectParameters.json
parameter_file = open("ProjectParameters.json",'r')
ProjectParameters = Parameters( parameter_file.read())


## Fluid-Structure model parts definition 
structure_main_model_part = ModelPart(ProjectParameters["structure_solver_settings"]["problem_data"]["model_part_name"].GetString())
structure_main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, ProjectParameters["structure_solver_settings"]["problem_data"]["domain_size"].GetInt())

fluid_main_model_part = ModelPart(ProjectParameters["fluid_solver_settings"]["problem_data"]["model_part_name"].GetString())
fluid_main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, ProjectParameters["fluid_solver_settings"]["problem_data"]["domain_size"].GetInt())

###TODO replace this "model" for real one once available
FluidModel = {ProjectParameters["fluid_solver_settings"]["problem_data"]["model_part_name"].GetString() : fluid_main_model_part}
SolidModel = {ProjectParameters["structure_solver_settings"]["problem_data"]["model_part_name"].GetString() : structure_main_model_part}

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

gid_output_structure = GiDOutputProcess(solver.structure_solver.GetComputingModelPart(),
                                    ProjectParameters["structure_solver_settings"]["problem_data"]["problem_name"].GetString()+"_structure",
                                    ProjectParameters["structure_solver_settings"]["output_configuration"])

gid_output_fluid = GiDOutputProcess(solver.fluid_solver.GetComputingModelPart(),
                                    ProjectParameters["fluid_solver_settings"]["problem_data"]["problem_name"].GetString()+"_fluid",
                                    ProjectParameters["fluid_solver_settings"]["output_configuration"])

gid_output_structure.ExecuteInitialize()
gid_output_fluid.ExecuteInitialize()

#~ ##here all of the allocation of the strategies etc is done
#~ solver.Initialize()


##TODO: replace MODEL for the Kratos one ASAP
## Get the list of the skin submodel parts in the object Model (FLUID)
for i in range(ProjectParameters["fluid_solver_settings"]["solver_settings"]["skin_parts"].size()):
    skin_part_name = ProjectParameters["fluid_solver_settings"]["solver_settings"]["skin_parts"][i].GetString()
    FluidModel.update({skin_part_name: fluid_main_model_part.GetSubModelPart(skin_part_name)})
    
## Get the list of the no-skin submodel parts in the object Model (FLUID)
for i in range(ProjectParameters["fluid_solver_settings"]["solver_settings"]["no_skin_parts"].size()):
    no_skin_part_name = ProjectParameters["fluid_solver_settings"]["solver_settings"]["no_skin_parts"][i].GetString()
    FluidModel.update({no_skin_part_name: fluid_main_model_part.GetSubModelPart(no_skin_part_name)})

## Get the list of the initial conditions submodel parts in the object Model (FLUID)
for i in range(ProjectParameters["fluid_solver_settings"]["initial_conditions_process_list"].size()):
    initial_cond_part_name = ProjectParameters["fluid_solver_settings"]["initial_conditions_process_list"][i]["Parameters"]["model_part_name"].GetString()
    FluidModel.update({initial_cond_part_name: fluid_main_model_part.GetSubModelPart(initial_cond_part_name)})
    
## Get the gravity submodel part in the object Model (FLUID)
for i in range(ProjectParameters["fluid_solver_settings"]["gravity"].size()):   
    gravity_part_name = ProjectParameters["fluid_solver_settings"]["gravity"][i]["Parameters"]["model_part_name"].GetString()
    FluidModel.update({gravity_part_name: fluid_main_model_part.GetSubModelPart(gravity_part_name)})

## Get the list of the submodel part in the object Model (STRUCTURE)
for i in range(ProjectParameters["structure_solver_settings"]["solver_settings"]["processes_sub_model_part_list"].size()):
    part_name = ProjectParameters["structure_solver_settings"]["solver_settings"]["processes_sub_model_part_list"][i].GetString()
    SolidModel.update({part_name: structure_main_model_part.GetSubModelPart(part_name)})


## Processes construction    
import process_factory
# "list_of_processes" contains all the processes already constructed (boundary conditions, initial conditions and gravity) 
# Note that the conditions are firstly constructed. Otherwise, they may overwrite the BCs information.

# FLUID DOMAIN PROCESSES
list_of_processes = process_factory.KratosProcessFactory(FluidModel).ConstructListOfProcesses( ProjectParameters["fluid_solver_settings"]["initial_conditions_process_list"] )
list_of_processes += process_factory.KratosProcessFactory(FluidModel).ConstructListOfProcesses( ProjectParameters["fluid_solver_settings"]["boundary_conditions_process_list"] )
list_of_processes += process_factory.KratosProcessFactory(FluidModel).ConstructListOfProcesses( ProjectParameters["fluid_solver_settings"]["gravity"] )

# SOLID DOMAIN PROCESSES
list_of_processes += process_factory.KratosProcessFactory(SolidModel).ConstructListOfProcesses( ProjectParameters["structure_solver_settings"]["constraints_process_list"] )
list_of_processes += process_factory.KratosProcessFactory(SolidModel).ConstructListOfProcesses( ProjectParameters["structure_solver_settings"]["loads_process_list"] )


## Processes initialization
for process in list_of_processes:
    process.ExecuteInitialize()
    
    
# Solver initialization moved after the processes initialization, otherwise the flag INTERFACE is not set 
solver.Initialize()


## Stepping and time settings
Dt = ProjectParameters["fluid_solver_settings"]["problem_data"]["time_step"].GetDouble()
end_time = ProjectParameters["fluid_solver_settings"]["problem_data"]["end_time"].GetDouble()

time = 0.0
step = 0
out = 0.0

gid_output_structure.ExecuteBeforeSolutionLoop()
gid_output_fluid.ExecuteBeforeSolutionLoop()

for process in list_of_processes:
    process.ExecuteBeforeSolutionLoop()

while(time <= end_time):

    time = time + Dt
    step = step + 1
    
    # Custom velocity profile for Mok benchmark
    if time <=10:
        v_bar = 0.5*0.06067*(1.0-numpy.cos(0.1*numpy.pi*time))
    elif time>10:
        v_bar = 0.06067
    
    for node in solver.fluid_solver.main_model_part.GetSubModelPart("Inlet2D_Inlet").Nodes:
        
        vel = Vector(3)
        vel[0] = 4*v_bar*node.Y*(1-node.Y)
        vel[1] = 0.0
        vel[2] = 0.0
        
        node.SetSolutionStepValue(VELOCITY,0,vel)
        node.Fix(VELOCITY_X)
        node.Fix(VELOCITY_Y)
        node.Fix(VELOCITY_Z)

    solver.SetTimeStep(step)
    
    structure_main_model_part.CloneTimeStep(time)    
    fluid_main_model_part.CloneTimeStep(time)   

    print("STEP = ", step)
    print("TIME = ", time)

    #~ if(step >= 3):
    for process in list_of_processes:
        process.ExecuteInitializeSolutionStep()
    
    gid_output_structure.ExecuteInitializeSolutionStep()
    gid_output_fluid.ExecuteInitializeSolutionStep()
    
    print("Time step ",step," resolution starts...")
    solver.Solve()
    print("Time step ",step," solved.")
    
    for process in list_of_processes:
        process.ExecuteFinalizeSolutionStep()
        
    gid_output_structure.ExecuteFinalizeSolutionStep()
    gid_output_fluid.ExecuteFinalizeSolutionStep()

    #TODO: decide if it shall be done only when output is processed or not
    for process in list_of_processes:
        process.ExecuteBeforeOutputStep()

    if gid_output_structure.IsOutputStep():
        gid_output_structure.PrintOutput()
        
    if gid_output_fluid.IsOutputStep():
        gid_output_fluid.PrintOutput()
    
    for process in list_of_processes:
        process.ExecuteAfterOutputStep()

    out = out + Dt

for process in list_of_processes:
    process.ExecuteFinalize()
    
gid_output_structure.ExecuteFinalize()
gid_output_fluid.ExecuteFinalize()
