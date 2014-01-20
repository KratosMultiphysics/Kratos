from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# -*- coding: utf-8 -*-
import re
import math
import sys
import ProjectParameters


#
#
# setting the domain size for the problem to be solved
domain_size = ProjectParameters.domain_size
# importing Kratos main library

from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.MeshingApplication import *

#
#

# defining a model part for the fluid and one for the structure
fluid_model_part = ModelPart("FluidPart")
#


# importing the solvers needed
import edgebased_levelset_substep_solver
edgebased_levelset_substep_solver.AddVariables(fluid_model_part)
fluid_model_part.AddNodalSolutionStepVariable(Y_WALL)
fluid_model_part.AddNodalSolutionStepVariable(DENSITY)

# introducing input file name
input_file_name = ProjectParameters.problem_name

# reading the fluid part
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.SingleFile
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly

gid_io = GidIO(
    input_file_name +
    "_F_k",
    gid_mode,
    multifile,
    deformed_mesh_flag,
    write_conditions)

model_part_io_fluid = ModelPartIO(input_file_name)
model_part_io_fluid.ReadModelPart(fluid_model_part)

# setting up the buffer size: SHOULD BE DONE AFTER READING!!!
fluid_model_part.SetBufferSize(2)

# adding dofs
edgebased_levelset_substep_solver.AddDofs(fluid_model_part)


# we assume here that all of the internal nodes are marked with a negative distance
# set the distance of all of the internal nodes to a small value
small_value = 0.0001
n_active = 0
for node in fluid_model_part.Nodes:
    dist = node.GetSolutionStepValue(DISTANCE)
    if(dist < 0.0):
        n_active = n_active + 1
        node.SetSolutionStepValue(DISTANCE, 0, -small_value)
    else:
        node.SetSolutionStepValue(DISTANCE, 0, small_value)

if(n_active == 0):
    raise "ERROR. At least one node has to be initialized with a distance lesser than 0"

# make sure that the porosity is not zero on any node (set by default to
# fluid only)
for node in fluid_model_part.Nodes:
    if(node.GetSolutionStepValue(POROSITY) == 0.0):
        node.SetSolutionStepValue(POROSITY, 0, 1.0)
    if(node.GetSolutionStepValue(DIAMETER) == 0.0):
        node.SetSolutionStepValue(DIAMETER, 0, 1.0)


# constructing the solver
body_force = Vector(3)
body_force[0] = ProjectParameters.body_force_x
body_force[1] = ProjectParameters.body_force_y
body_force[2] = ProjectParameters.body_force_z
viscosity = ProjectParameters.viscosity
density = ProjectParameters.density
fluid_solver = edgebased_levelset_substep_solver.EdgeBasedLevelSetSolver(
    fluid_model_part, domain_size, body_force, viscosity, density)
fluid_solver.redistance_frequency = ProjectParameters.redistance_frequency
fluid_solver.extrapolation_layers = ProjectParameters.extrapolation_layers
fluid_solver.stabdt_pressure_factor = ProjectParameters.stabdt_pressure_factor
fluid_solver.stabdt_convection_factor = ProjectParameters.stabdt_convection_factor
fluid_solver.compute_porous_resistance_law = 0
fluid_solver.pressure_linear_solver = BICGSTABSolver(1e-6, 5000)

fluid_solver.Initialize()

print("***********fluid solver created****************")

if(ProjectParameters.wall_law_y > 1e-10):
    fluid_solver.fluid_solver.ActivateWallResistance(
        ProjectParameters.wall_law_y)

#


print("fluid solver created")

# settings to be changed
max_Dt = ProjectParameters.Dt
initial_Dt = 0.001 * max_Dt
final_time = ProjectParameters.max_time
output_dt = ProjectParameters.output_time
safety_factor = ProjectParameters.safety_factor

number_of_inital_steps = ProjectParameters.number_of_initial_steps
initial_time_step = initial_Dt
out = 0

original_max_dt = max_Dt

# mesh to be printed
mesh_name = 0.0
gid_io.InitializeMesh(mesh_name)
gid_io.WriteMesh(fluid_model_part.GetMesh())
gid_io.FinalizeMesh()
gid_io.Flush()

gid_io.InitializeResults(mesh_name, (fluid_model_part).GetMesh())

normal_calculator = NormalCalculationUtils()
normal_calculator.CalculateOnSimplex(fluid_model_part.Conditions, 3)


for condition in fluid_model_part.Conditions:
    if(condition.GetValue(IS_INLET) > 0.00):
        normal = condition.GetValue(NORMAL)
        normal_size = math.sqrt(normal * normal)
        velocity = ProjectParameters.inlet_velocity * normal / normal_size
        for node in condition.GetNodes():
            node.SetSolutionStepValue(VELOCITY, velocity)
            node.SetSolutionStepValue(
                TEMPERATURE,
                ProjectParameters.FLUID_TEMPERATURE)
            node.Fix(VELOCITY_X)
            node.Fix(VELOCITY_Y)
            node.Fix(VELOCITY_Z)
            node.Fix(DISTANCE)
            node.Fix(TEMPERATURE)
            node.SetValue(IS_VISITED, 10.0)
            # print condition.Id, node.GetSolutionStepValue(VELOCITY)

inlet_nodes_list = []
for node in fluid_model_part.Nodes:
    if (node.IsFixed(VELOCITY_X)):
        inlet_nodes_list.append(node)
        node.Fix(DISTANCE)

fixed_dist_nodes = []
for node in fluid_model_part.Nodes:
    if(node.IsFixed(DISTANCE)):
        fixed_dist_nodes.append(node)


def WettenNodes(nodes):
    for node in nodes:
        if(node.GetSolutionStepValue(DISTANCE) > 0):
            node.SetSolutionStepValue(DISTANCE, 0, -0.001)


def FreeFixedInletValues(model_part):
    for node in model_part.Nodes:
        node.Free(VELOCITY_X)
        node.Free(VELOCITY_Y)
        node.Free(VELOCITY_Z)
        node.Free(TEMPERATURE)


max_safety_factor = safety_factor

time = 0.0
step = 0
next_output_time = output_dt
screen_output_dt = 0.1
next_screen_output = screen_output_dt
volume_correction_step = 1

print("Process Information")
print("---------------------------------------------------------------")
print("Max time          :", final_time)
print("Max delta time    :", max_Dt)
print("Output delta time :", output_dt)
print("Safety factor     :", safety_factor)


print("Filled %        current time     delta time      mass ratio")
print("---------------------------------------------------------------")
sys.stdout.flush()

measured_volume = fluid_solver.fluid_solver.ComputeWetVolume()
expected_volume = measured_volume

time1 = time
# AssignEnvironmentCondition.AssignCondition()
switch = 1.0
temp_time = 0.0
# mean_prerssure_file = open("mean_prerssure_file.out", "w")
while((time1 < final_time)):

    if(step < number_of_inital_steps):
        max_Dt = initial_time_step
    else:
        max_Dt = original_max_dt
        # progressively increment the safety factor
        # in the steps that follow a reduction of it
        safety_factor = safety_factor * 1.2
        if(safety_factor > max_safety_factor):
            safety_factor = max_safety_factor

    Dt = fluid_solver.EstimateTimeStep(safety_factor, max_Dt)
    time = time + Dt
    time1 = time1 + Dt
    percent_done = 100.00 * (time1 / final_time)

    fluid_model_part.CloneTimeStep(time)

    if(step >= 3):
        # print "time=",time," Dt = ",Dt
        WettenNodes(fixed_dist_nodes)
        fluid_solver.Solve()

    measured_volume = fluid_solver.measured_volume
    vol_variation = fluid_solver.fluid_solver.ComputeVolumeVariation()
    expected_volume = fluid_solver.expected_volume

    if(percent_done >= next_screen_output):
        print()
        print("Filled %.0f" % percent_done, "% \t", "%e" % time1, "\t", "%e" % Dt, "\t", measured_volume / expected_volume)
        sys.stdout.flush()
        next_screen_output += screen_output_dt

# I have to change this part to not duplicate the redistance. Pooyan.
    # if(volume_correction_step > ProjectParameters.redistance_frequency):
        # max_volume_error = 0.999
        # if(measured_volume / expected_volume < max_volume_error):
            # vol_variation =  fluid_solver.fluid_solver.ContinuousVolumeCorrection(expected_volume, measured_volume)
            # fluid_solver.Redistance();
            # volume_correction_step = 1
        # if(measured_volume / expected_volume > 1.00):
            # time1 = time * measured_volume / expected_volume # This is NOT
            # what I like to do! Pooyan.
    # volume_correction_step += 1
    if(time1 >= next_output_time):

        gid_io.WriteNodalResults(DISTANCE, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(VELOCITY, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(PRESSURE, fluid_model_part.Nodes, time, 0)

        gid_io.Flush()
        sys.stdout.flush()

        next_output_time += output_dt
        out = 0

    out = out + 1
    step = step + 1

gid_io.FinalizeResults()

print("Num Steps:", step)
print("END OF RUN EXECUTION")
