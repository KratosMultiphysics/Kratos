from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
from KratosMultiphysics.ALEApplication import *
from KratosMultiphysics.FSIApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *

######################################################################################
######################################################################################
######################################################################################

## Parsing ProjectParameters.json
parameter_file = open("ProjectParameters.json",'r')
ProjectParameters = Parameters( parameter_file.read())

## Get echo level and parallel type
verbosity_fluid = ProjectParameters["fluid_solver_settings"]["problem_data"]["echo_level"].GetInt()
verbosity_structure = ProjectParameters["structure_solver_settings"]["problem_data"]["echo_level"].GetInt()
verbosity = max(verbosity_fluid, verbosity_structure)
parallel_type = ProjectParameters["coupling_solver_settings"]["problem_data"]["parallel_type"].GetString()

## Import KratosMPI if needed
if (parallel_type == "MPI"):
    import KratosMultiphysics.mpi as KratosMPI

## Fluid-Structure model parts definition
structure_main_model_part = ModelPart(ProjectParameters["structure_solver_settings"]["problem_data"]["model_part_name"].GetString())
structure_main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, ProjectParameters["structure_solver_settings"]["problem_data"]["domain_size"].GetInt())

fluid_main_model_part = ModelPart(ProjectParameters["fluid_solver_settings"]["problem_data"]["model_part_name"].GetString())
fluid_main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, ProjectParameters["fluid_solver_settings"]["problem_data"]["domain_size"].GetInt())

## Solver construction
solver_module = __import__(ProjectParameters["coupling_solver_settings"]["solver_settings"]["solver_type"].GetString())
solver = solver_module.CreateSolver(structure_main_model_part, fluid_main_model_part, ProjectParameters)

solver.AddVariables()

## Read the model - note that SetBufferSize is done here
solver.ImportModelPart()

## Add AddDofs
solver.AddDofs()

## Initialize GiD  I/O
if (parallel_type == "OpenMP"):
    from gid_output_process import GiDOutputProcess
    gid_output_structure = GiDOutputProcess(solver.structure_solver.GetComputingModelPart(),
                                        ProjectParameters["structure_solver_settings"]["problem_data"]["problem_name"].GetString()+"_structure",
                                        ProjectParameters["structure_solver_settings"]["output_configuration"])
    gid_output_fluid = GiDOutputProcess(solver.fluid_solver.GetComputingModelPart(),
                                        ProjectParameters["fluid_solver_settings"]["problem_data"]["problem_name"].GetString()+"_fluid",
                                        ProjectParameters["fluid_solver_settings"]["output_configuration"])

elif (parallel_type == "MPI"):
    from gid_output_process_mpi import GiDOutputProcessMPI
    gid_output_structure = GiDOutputProcessMPI(solver.structure_solver.GetComputingModelPart(),
                                               ProjectParameters["structure_solver_settings"]["problem_data"]["problem_name"].GetString()+"_structure",
                                               ProjectParameters["structure_solver_settings"]["output_configuration"])

    gid_output_fluid = GiDOutputProcessMPI(solver.fluid_solver.GetComputingModelPart(),
                                           ProjectParameters["fluid_solver_settings"]["problem_data"]["problem_name"].GetString()+"_fluid",
                                           ProjectParameters["fluid_solver_settings"]["output_configuration"])

gid_output_structure.ExecuteInitialize()
gid_output_fluid.ExecuteInitialize()

## Creation of Kratos models
FluidModel = Model()
FluidModel.AddModelPart(fluid_main_model_part)
SolidModel = Model()
SolidModel.AddModelPart(structure_main_model_part)

## Get the list of the skin submodel parts in the object Model
for i in range(ProjectParameters["fluid_solver_settings"]["solver_settings"]["skin_parts"].size()):
    skin_part_name = ProjectParameters["fluid_solver_settings"]["solver_settings"]["skin_parts"][i].GetString()
    FluidModel.AddModelPart(fluid_main_model_part.GetSubModelPart(skin_part_name))

## Get the list of the no-skin submodel parts in the object Model (results processes and no-skin conditions)
for i in range(ProjectParameters["fluid_solver_settings"]["solver_settings"]["no_skin_parts"].size()):
    no_skin_part_name = ProjectParameters["fluid_solver_settings"]["solver_settings"]["no_skin_parts"][i].GetString()
    FluidModel.AddModelPart(fluid_main_model_part.GetSubModelPart(no_skin_part_name))

## Get the list of the initial conditions submodel parts in the object Model
for i in range(ProjectParameters["fluid_solver_settings"]["initial_conditions_process_list"].size()):
    initial_cond_part_name = ProjectParameters["fluid_solver_settings"]["initial_conditions_process_list"][i]["Parameters"]["model_part_name"].GetString()
    FluidModel.AddModelPart(fluid_main_model_part.GetSubModelPart(initial_cond_part_name))

## Get the gravity submodel part in the object Model
for i in range(ProjectParameters["fluid_solver_settings"]["gravity"].size()):
    gravity_part_name = ProjectParameters["fluid_solver_settings"]["gravity"][i]["Parameters"]["model_part_name"].GetString()
    FluidModel.AddModelPart(fluid_main_model_part.GetSubModelPart(gravity_part_name))

## Get the list of the submodel part in the object Model (STRUCTURE)
for i in range(ProjectParameters["structure_solver_settings"]["solver_settings"]["processes_sub_model_part_list"].size()):
    part_name = ProjectParameters["structure_solver_settings"]["solver_settings"]["processes_sub_model_part_list"][i].GetString()
    SolidModel.AddModelPart(structure_main_model_part.GetSubModelPart(part_name))

# Print model_parts and properties
if(verbosity_fluid > 1):
    print("")
    print(fluid_main_model_part)
    for properties in fluid_main_model_part.Properties:
        print(properties)

if(verbosity_structure > 1):
    print("")
    print(structure_main_model_part)
    for properties in structure_main_model_part.Properties:
        print(properties)

## Processes construction
import process_factory
# "list_of_processes" contains all the processes already constructed (boundary conditions, initial conditions and gravity)
# Note that the conditions are firstly constructed. Otherwise, they may overwrite the BCs information.

# FLUID DOMAIN PROCESSES
list_of_processes  = process_factory.KratosProcessFactory(FluidModel).ConstructListOfProcesses( ProjectParameters["fluid_solver_settings"]["gravity"] )
list_of_processes += process_factory.KratosProcessFactory(FluidModel).ConstructListOfProcesses( ProjectParameters["fluid_solver_settings"]["initial_conditions_process_list"] )
list_of_processes += process_factory.KratosProcessFactory(FluidModel).ConstructListOfProcesses( ProjectParameters["fluid_solver_settings"]["boundary_conditions_process_list"] )
list_of_processes += process_factory.KratosProcessFactory(FluidModel).ConstructListOfProcesses( ProjectParameters["fluid_solver_settings"]["auxiliar_process_list"] )

# SOLID DOMAIN PROCESSES
list_of_processes += process_factory.KratosProcessFactory(SolidModel).ConstructListOfProcesses( ProjectParameters["structure_solver_settings"]["constraints_process_list"] )
list_of_processes += process_factory.KratosProcessFactory(SolidModel).ConstructListOfProcesses( ProjectParameters["structure_solver_settings"]["loads_process_list"] )

if(verbosity > 1):
    for process in list_of_processes:
        print(process)

## Processes initialization
for process in list_of_processes:
    process.ExecuteInitialize()

# Solver initialization moved after the processes initialization, otherwise the flag INTERFACE is not set
solver.Initialize()

## Stepping and time settings
end_time = ProjectParameters["fluid_solver_settings"]["problem_data"]["end_time"].GetDouble()
time = 0.0
step = 0
out = 0.0

gid_output_structure.ExecuteBeforeSolutionLoop()
gid_output_fluid.ExecuteBeforeSolutionLoop()

for process in list_of_processes:
    process.ExecuteBeforeSolutionLoop()

## Writing the full ProjectParameters file before solving
if ((parallel_type == "OpenMP") or (KratosMPI.mpi.rank == 0)) and (verbosity > 0):
    f = open("ProjectParametersOutput.json", 'w')
    f.write(ProjectParameters.PrettyPrintJsonString())
    f.close()

while(time <= end_time):

    Dt = solver.ComputeDeltaTime()
    time = time + Dt
    step = step + 1

    solver.SetTimeStep(step)

    structure_main_model_part.CloneTimeStep(time)
    fluid_main_model_part.CloneTimeStep(time)

    if (parallel_type == "OpenMP") or (KratosMPI.mpi.rank == 0):
        print("STEP = ", step)
        print("TIME = ", time)

    for process in list_of_processes:
        process.ExecuteInitializeSolutionStep()

    gid_output_structure.ExecuteInitializeSolutionStep()
    gid_output_fluid.ExecuteInitializeSolutionStep()

    if (parallel_type == "OpenMP") or (KratosMPI.mpi.rank == 0):
        print("Time step ",step," resolution starts...")

    solver.Solve()

    if (parallel_type == "OpenMP") or (KratosMPI.mpi.rank == 0):
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
