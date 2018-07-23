from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from math import sin,cos

#Definition of functions
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

movementList = ["TRANSLATION_X", "TRANSLATION_Y", "BENDING", "ROTATION"]
# =============================================================
movement = 2             # set 0-3
amplificationFactor = 2  # factor to amplify the movement
solveMesh = True         # set to "True" to solve the mesh
# =============================================================


#Apply displacement boundary conditions

def ApplyDisplacementConditions(main_model_part):
    for node in main_model_part.GetSubModelPart("Inlet2D_Inlet").Nodes:
        node.Fix(MESH_DISPLACEMENT_X)
        node.Fix(MESH_DISPLACEMENT_Y)
    for node in main_model_part.GetSubModelPart("Outlet2D_Outlet").Nodes:
        node.Fix(MESH_DISPLACEMENT_Y)
        node.Fix(MESH_DISPLACEMENT_X)
    for node in main_model_part.GetSubModelPart("Slip2D_Walls").Nodes:
        node.Fix(MESH_DISPLACEMENT_Y)
        node.Fix(MESH_DISPLACEMENT_X)
    for node in main_model_part.GetSubModelPart("NoSlip2D_Interface").Nodes:
        node.Fix(MESH_DISPLACEMENT_X)
        node.Fix(MESH_DISPLACEMENT_Y)

def DisplacementToMesh(fluid_interface, time, movement, ampFac):
    time *= ampFac
    for node in fluid_interface:
        if movement == "TRANSLATION_X":
            valueX = 0.5 * time
            valueY = 0
        elif movement == "TRANSLATION_Y":
            valueX = 0
            valueY = 0.2 * time
        elif movement == "BENDING":
            valueX = 0
            valueY = 3 * node.X0 * node.X0 * time
        elif movement == "ROTATION":
            xOld = node.X0
            yOld = node.Y0
            valueX = (cos(time*.6)*xOld + sin(time*.6)*yOld) - xOld
            valueY = (-sin(time*.6)*xOld + cos(time*.6)*yOld) - yOld

        else:
            wait = input("Wrong type of movement specified, please correct input")
        # set the prescribed values
        node.SetSolutionStepValue(MESH_DISPLACEMENT_X,0,valueX)
        node.SetSolutionStepValue(MESH_DISPLACEMENT_Y,0,valueY)


#Set number of threads

def SetParallelSize(num_threads):
    parallel = OpenMPUtils()
    print("Num Threads = ", num_threads)
    parallel.SetNumThreads(int(num_threads))

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
from KratosMultiphysics import *
from KratosMultiphysics.MeshMovingApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.MeshingApplication import *

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
#solver_module = __import__(ProjectParameters["solver_settings"]["solver_type"].GetString())
#solver = solver_module.CreateSolver(main_model_part, ProjectParameters["solver_settings"])

import python_solvers_wrapper_mesh_motion
solver = python_solvers_wrapper_mesh_motion.CreateSolver(main_model_part, ProjectParameters)

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


## Stepping and time settings
Dt = ProjectParameters["solver_settings"]["time_stepping"]["time_step"].GetDouble()
end_time = ProjectParameters["problem_data"]["end_time"].GetDouble()

time = 0.0
step = 0
out = 0.0

gid_output.ExecuteBeforeSolutionLoop()


ApplyDisplacementConditions(main_model_part)

while(time <= end_time):

    time = time + Dt
    step = step + 1
    main_model_part.CloneTimeStep(time)

    print("STEP = ", step)
    print("TIME = ", time)

    if(step >= 3):
        #for process in list_of_processes:
        #    process.ExecuteInitializeSolutionStep()

        #gid_output.ExecuteInitializeSolutionStep()
        DisplacementToMesh(main_model_part.GetSubModelPart("NoSlip2D_Interface").Nodes, time, movementList[movement], amplificationFactor)

        solver.Solve()

        #for process in list_of_processes:
        #    process.ExecuteFinalizeSolutionStep()

        gid_output.ExecuteFinalizeSolutionStep()

        #TODO: decide if it shall be done only when output is processed or not
        #for process in list_of_processes:
        #    process.ExecuteBeforeOutputStep()

        if gid_output.IsOutputStep():
            gid_output.PrintOutput()

        #for process in list_of_processes:
        #    process.ExecuteAfterOutputStep()

        out = out + Dt

#for process in list_of_processes:
#    process.ExecuteFinalize()

gid_output.ExecuteFinalize()













