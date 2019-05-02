from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MeshingApplication import *

# import own files
from python_solver.structure.structure_beam import *
from python_solver.mapper.mapping import *
from python_solver.convergence.Residual import *

# for FSI - mesh moving
from KratosMultiphysics.ALEApplication import *

######################################################################################
######################################################################################
######################################################################################
##PARSING THE PARAMETERS
#import define_output

parameter_file = open("ProjectParameters_Custom.json",'r')
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
gid_output = GiDOutputProcess(main_model_part,
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
list_of_processes += process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["compute_total_force"] )
list_of_processes += process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["compute_level_force"] )
list_of_processes += process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["aerodynamic_forces"] )
list_of_processes += process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["perimeter_results"] )
list_of_processes += process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["top_results"] )
( ProjectParameters["list_of_line_outputs"] )
list_of_processes += process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["list_of_point_outputs"] )
## Processes initialization
for process in list_of_processes:
    process.ExecuteInitialize()

# ALE boundary conditions
for i in range(ProjectParameters["solver_settings"]["skin_parts"].size()):
    for node in main_model_part.GetSubModelPart(ProjectParameters["solver_settings"]["skin_parts"][i].GetString()).Nodes:
        node.Fix(MESH_DISPLACEMENT_X)
        node.Fix(MESH_DISPLACEMENT_Y)
        node.Fix(MESH_DISPLACEMENT_Z)

#TODO: think if there is a better way to do this
#fluid_model_part = solver.GetComputeModelPart()

structure = Structure(ProjectParameters)
mapper = Mapper(main_model_part.GetSubModelPart("NoSlip3D_structure"), structure)
solution = Convergence(structure)

## Stepping and time settings
Dt = ProjectParameters["problem_data"]["time_step"].GetDouble()
end_time = ProjectParameters["problem_data"]["end_time"].GetDouble()

time = 0.0
step = 0
out = 0.0

gid_output.ExecuteBeforeSolutionLoop()

for process in list_of_processes:
    process.ExecuteBeforeSolutionLoop()

import sys
structure.predict_displacement()
while(time <= end_time):

    time = time + Dt
    step = step + 1
    main_model_part.CloneTimeStep(time)

    print("STEP = ", step)
    print("TIME = ", time)
    sys.stdout.flush()

    if(step >= 3):
        for process in list_of_processes:
            process.ExecuteInitializeSolutionStep()

        gid_output.ExecuteInitializeSolutionStep()
        initial_residual = []
        for k in range(0, structure.properties.fsi_max_iter):

            # Set Mesh displacement from Structure
            mapper.set_mesh_displacement()
            print("MESH MAPPING DONE!")

            # Solve mesh motion
            solver.SolveMeshMotion()
            print("MESH MOTION SOLVED!")

            # Apply Mesh Velocity to fluid solver
            mapper.set_mesh_velocity_to_fluid()
            print("MAPPING MESH TO FLUID DONE!")

             # Solve Fluid
            solver.SolveFluid()
            print("FLUID SOLVE DONE!")
            # Apply Neumann B.C.'s from fluid to structure
            mapper.extract_forces()
            print("MAPPER EXTRACT FORCES DOEN!")
            mapper.map_forces_to_structure()
            print("FORCES MAPPED TO STRUCTURE DONE!")
            # Solve structure
            structure.solve(mapper.mapped_forces)
            print("STRUCURAL SOLVER DONE")
            # Get solution from structure
            structure.get_displacement()

            # Calculate residual
            solution.cal_residual(structure)
            initial_residual.append(np.linalg.norm(solution.residual))
            # Check convergence critera
            if k==0 and (np.linalg.norm(solution.residual) <= solution.abs_residual):
                print("CONVERGENCE AT INTERFACE ACHIEVED")
                # Update structural results for convergence
                structure.update_result()
                break
            if k>0 and (np.linalg.norm(solution.residual) <= solution.abs_residual or np.linalg.norm(solution.residual) < initial_residual[0]*solution.rel_residual ):
                print("CONVERGENCE AT INTERFACE ACHIEVED")
                # Update structural results for convergence
                structure.update_result()
                break
            else:
                # solution.aitken_relaxation(k) # compute new relaxation coefficient
                print("RELAXATION COEFFICIENT: ", solution.relax_coef)
                solution.cal_relaxation(structure)
                print(
                    'ITERATION [', k, ']: RESIDUAL = ', np.linalg.norm(solution.residual))

                # Update structural results for convergence
                structure.update_relaxed_result(solution.relaxed_solution)

        # Print structural results
        structure.print_support_output(time)
        # Get back the reactions
        structure.get_forces_back(time)
        # Update structure time
        structure.update_structure_time()

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













