from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import problem_settings

import time as timer
import sys

from KratosMultiphysics import *
from KratosMultiphysics.MeshlessApplication import *
from KratosMultiphysics.ExternalSolversApplication import *


start = timer.time()
print("Timer started")


# defining a model part
solid_model_part = ModelPart("SolidPart")
#
domain_size = problem_settings.domain_size
problem_name = problem_settings.problem_name


# importing the solver files and adding the variables
SolverType = problem_settings.SolverType

if(SolverType == "weakly_compressible_SPH"):
    import meshless_solverWCSPHnew
    meshless_solverWCSPHnew.AddVariables(solid_model_part)
elif(SolverType == "incompressible_SPH"):
    import meshless_solverISPHnew
    meshless_solverISPHnew.AddVariables(solid_model_part)
else:
    raise "solver type not supported: options are WCSPH - ISPH"


# reading the Model Part
gid_mode = GiDPostMode.GiD_PostBinary    # or GiDPostMode.GiD_PostAscii
use_multi_file = MultiFileFlag.MultipleFiles    # or MultiFileFlag.SingleFile
deformed_mesh_flag = WriteDeformedMeshFlag.WriteDeformed    # or WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly   # or WriteConditionsFlag.WriteConditions

gid_io = GidIO(problem_name, gid_mode, use_multi_file, deformed_mesh_flag, write_conditions)
model_part_io_solid = ModelPartIO(problem_name)
model_part_io_solid.ReadModelPart(solid_model_part)

# assign the property variables to the nodes
gravity = Vector(3)  # (0.0,-9.81,0.0)
gravity[0] = problem_settings.gravity_x
gravity[1] = problem_settings.gravity_y
gravity[2] = problem_settings.gravity_z

# The element
element_choice = problem_settings.element_choice

if element_choice == 1:
    real_search_radius = 1.001 * problem_settings.search_radius

if element_choice == 2:
    real_search_radius = 2.001 * problem_settings.search_radius

if element_choice == 3:
    real_search_radius = 3.001 * problem_settings.search_radius

if element_choice == 4:
    real_search_radius = 2.001 * problem_settings.search_radius

if element_choice == 5:
    real_search_radius = 3.001 * problem_settings.search_radius

if element_choice == 6:
    real_search_radius = 2.001 * problem_settings.search_radius

for node in solid_model_part.Nodes:
    node.SetSolutionStepValue(BODYFORCE_ACC, 0, gravity)
    node.SetSolutionStepValue(VISCOSITY, 0, problem_settings.Viscosity)
    node.SetSolutionStepValue(SEARCH_RADIUS, 0, real_search_radius)  # important!!!!
    node.SetSolutionStepValue(EFFECTIVE_RADIUS, 0, problem_settings.search_radius)  # important!!!!
    node.SetSolutionStepValue(DENSITY, 0, problem_settings.ReferenceDensity)
    node.SetSolutionStepValue(RADIUS, 0, 0)


solid_model_part.ProcessInfo.SetValue(SOUND_VELOCITY, problem_settings.soundspeed)


# Write the very first mesh
mesh_name = 0.0
gid_io.InitializeMesh(mesh_name)
gid_io.WriteMesh((solid_model_part).GetMesh())
gid_io.FinalizeMesh()


# setting the limits of the bounding box
box_corner1 = Vector(3)
box_corner1[0] = problem_settings.min_x
box_corner1[1] = problem_settings.min_y
box_corner1[2] = problem_settings.min_z
box_corner2 = Vector(3);
box_corner2[0] = problem_settings.max_x;
box_corner2[1] = problem_settings.max_y;
box_corner2[2] = problem_settings.max_z;


# setting up the buffer size: SHOULD BE DONE AFTER READING!!!
solid_model_part.SetBufferSize(3)

# Add Dofs

if(SolverType == "weakly_compressible_SPH"):
    meshless_solverWCSPHnew.AddDofs(solid_model_part)
elif(SolverType == "incompressible_SPH"):
    meshless_solverISPHnew.AddDofs(solid_model_part)
else:
    raise "solver type not supported: options are WCSPH - ISPH"


# The creation of the SPH particles
SetSPH = CreateSPHParticle(solid_model_part, domain_size, element_choice)
SetSPH.Execute()


# variables for the problem
output_step_count = problem_settings.output_step_count
time_step = problem_settings.time_step

solid_model_part.ProcessInfo.SetValue(DELTA_TIME, time_step)


final_time = problem_settings.max_time


soundspeed = problem_settings.soundspeed
ReferenceDensity = problem_settings.ReferenceDensity
ini_spacing = problem_settings.ini_spacing


# Load the solver python file
if(SolverType == "weakly_compressible_SPH"):
    solver = meshless_solverWCSPHnew.MeshlessSolverWCSPHnew(solid_model_part, domain_size, box_corner1, box_corner2)
elif(SolverType == "incompressible_SPH"):
    solver = meshless_solverISPHnew.MeshlessSolverISPHnew(solid_model_part, domain_size, box_corner1, box_corner2)
else:
    raise "solver type not supported: options are WCSPH - ISPH - PCISPH"


solver.Initialize()


print("11111111111111111111111")


# THE LOOP IN TIME##

time = 0.0
step = 1

if(SolverType == "weakly_compressible_SPH"):
    while(time < final_time):

        corrected_dt = solver.EstimateDeltaTime(time_step);
        print("Time step is:", corrected_dt)
        time = time + corrected_dt
        solid_model_part.CloneTimeStep(time)

        solver.Solve()

        if step % output_step_count == 1:
            solver.OutputToGID(gid_io, time)
            print("step is : ", step)

        print("End of Step number:", step)

        print("We are at time= ", time)
        step = step + 1
        if step == 500000000:
            break


elif(SolverType == "incompressible_SPH"):
    while(time < final_time):

        #~ solver.MarkOutsideNodes()
        #~ CheckBox=NodeAndElementEraseProcess(solid_model_part)
        #~ CheckBox.Execute()

        if step % output_step_count == 1:
            solver.OutputToGID(gid_io, time)
            print("step is : ", step)

        corrected_dt = solver.EstimateDeltaTime(time_step)
        print("Time step is:", corrected_dt)
        time = time + corrected_dt
        solid_model_part.CloneTimeStep(time)

        solver.Solve()

        print("End of Step number:", step)

        print("We are at time= ", time)
        step = step + 1
        if step == 100000:
            break


print("Solution Complete ! ")

print("Time Elapsed in seconds ")

print(str(timer.time() - start))
