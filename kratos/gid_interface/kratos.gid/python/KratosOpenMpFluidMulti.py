from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# -*- coding: utf-8 -*-
import sys
import math
#import time as timer
import time as time_measure
import ProjectParameters

print(sys.argv)
if len(sys.argv) < 2:
    print("ERROR: The configuration file is missing")
    print("")
    print("Syntax:")
    print("    KratosFastFill.exe configfile.conf")
    print("")
    print("example:")
    print("    KratosFastFill.exe ProjectParameters.conf")
    sys.exit()
# including kratos path
exec(compile(open(sys.argv[1]).read(), sys.argv[1], 'exec'))
# setting the domain size for the problem to be solved
domain_size = ProjectParameters.domain_size

start_time = time_measure.time()

# importing Kratos main library
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ThermoMechanicalApplication import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
##################################################################
# defining a model part for the fluid and one for the structure
fluid_model_part = ModelPart("FluidPart")
#############################################
# importing the solvers needed
import Kratos_dpg_monolithic_solver_eulerian_test
Kratos_dpg_monolithic_solver_eulerian_test.AddVariables(fluid_model_part)

# For thermal analysis
my_settings = ConvectionDiffusionSettings()

my_settings.SetDensityVariable(DENSITY)
my_settings.SetDiffusionVariable(CONDUCTIVITY)
my_settings.SetMeshVelocityVariable(MESH_VELOCITY)
my_settings.SetConvectionVariable(VELOCITY)


fluid_model_part.AddNodalSolutionStepVariable(MESH_VELOCITY)
fluid_model_part.AddNodalSolutionStepVariable(NODAL_MASS)
fluid_model_part.AddNodalSolutionStepVariable(REACTION)


# introducing input file name
input_file_name = ProjectParameters.problem_name

# reading the fluid part
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.SingleFile
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
##gid_io = GidIO(input_file_name,gid_mode,multifile,deformed_mesh_flag, write_conditions)
gid_io = GidIO(input_file_name + "_F_k", gid_mode, multifile, deformed_mesh_flag, write_conditions)
model_part_io_fluid = ModelPartIO(input_file_name)
model_part_io_fluid.ReadModelPart(fluid_model_part)

# setting up the buffer size: SHOULD BE DONE AFTER READING!!!
fluid_model_part.SetBufferSize(3)

# adding dofs
Kratos_dpg_monolithic_solver_eulerian_test.AddDofs(fluid_model_part)


# creating the solvers
# fluid solver wall_law_y
CFL = 10.0
volume_correction_switch = True
y_wall_val = 1.0 * ProjectParameters.wall_law_y
y_wall_fac = 1.0
inlet_velocity = 10.0


##air_temp = 0.5*(ProjectParameters.FLUID_TEMPERATURE + ProjectParameters.SOLID_TEMPERATURE)
for node in fluid_model_part.Nodes:
    node.SetValue(Y_WALL, y_wall_val)
    node.SetSolutionStepValue(BODY_FORCE_X, 0, ProjectParameters.body_force_x)
    node.SetSolutionStepValue(BODY_FORCE_Y, 0, ProjectParameters.body_force_y)
    node.SetSolutionStepValue(BODY_FORCE_Z, 0, ProjectParameters.body_force_z)
    node.SetSolutionStepValue(DISTANCE, 0, 1000.0)
    node.SetValue(IS_STRUCTURE, 0)
    node.SetSolutionStepValue(IS_SLIP, 0, 0.0)
    # node.SetSolutionStepValue(VELOCITY_X,0,0.0)
    # node.SetSolutionStepValue(VELOCITY_Y,0,0.0)
    # node.SetSolutionStepValue(VELOCITY_Z,0,0.0)
    # if(node.Has(VELOCITY_X)):
      # print node.GetSolutionStepValue(VELOCITY_X)
# for condition in fluid_model_part.Conditions:
   # To compute Normals and slip
   # for node in condition.GetNodes():
         # node.SetSolutionStepValue(IS_SLIP,0,1.0)
         # node.SetSolutionStepValue(IS_STRUCTURE,0,1.0)
         # node.SetValue(IS_STRUCTURE,1.0)
         # if(condition.GetValue(IS_INLET)> 0.0):
           # node.SetSolutionStepValue(DISTANCE,0,-100.0)
           # node.Fix(DISTANCE)


# importing the solvers needed
SolverSettings = ProjectParameters.FluidSolverConfiguration
#solver_module = import_solver(SolverSettings)

# Creating the fluid solver
fluid_solver = Kratos_dpg_monolithic_solver_eulerian_test.CreateSolver(fluid_model_part, SolverSettings)


def OpenAirExitInFarDryZone(model_part, bx, by, bz):
    body_force = Vector(3)
    body_force[0] = bx
    body_force[1] = by
    body_force[2] = bz
    zero = Vector(3)
    zero[0] = 0.0
    zero[1] = 0.0
    zero[2] = 0.0

    number_of_exits = 0
    max_d = 0.0
    for node in model_part.Nodes:
        d = node.GetSolutionStepValue(DISTANCE)
        if(d > max_d):
            max_d = d

    out_skin_list = []

    for node in model_part.Nodes:
        slip_flag = node.GetSolutionStepValue(IS_SLIP)
        dist = node.GetSolutionStepValue(DISTANCE)
        node.SetSolutionStepValue(IS_FREE_SURFACE, 0, 0.0)
        if(dist >= 0.9 * max_d):
            node.SetSolutionStepValue(BODY_FORCE, 0, zero)
            node.SetSolutionStepValue(IS_FREE_SURFACE, 0, 10.0)
            # node.SetSolutionStepValue(VELOCITY,0,zero)
            # node.SetSolutionStepValue(VELOCITY,1,zero)
            # node.SetSolutionStepValue(VELOCITY,2,zero)
            number_of_exits += 1

            node.Fix(PRESSURE)
            node.SetSolutionStepValue(PRESSURE, 0, 0.0)
            if (slip_flag != 0):

                node.SetValue(IS_STRUCTURE, 0.0)
        else:
            node.SetSolutionStepValue(BODY_FORCE, 0, body_force)
            node.Free(PRESSURE)
            if (slip_flag == 10.0 or slip_flag == 20.0 or slip_flag == 30.0 or slip_flag == 40.0):
                if(dist > 0.0):
                    if(slip_flag != 10):
                        node.Fix(PRESSURE)
                        node.SetSolutionStepValue(PRESSURE, 0, 0.0)
                        node.SetValue(IS_STRUCTURE, 0.0)
                    else:
                        out_skin_list.append(node)
                        node.SetValue(IS_STRUCTURE, 1.0)
                if(dist < 0.0):  # and slip_flag == 30.0
                    if(slip_flag != 40):
                        node.SetValue(IS_STRUCTURE, 1.0)
                    else:
                        node.SetValue(IS_STRUCTURE, 0.0)
                    if(slip_flag != 10 and slip_flag != 40):
                        node.SetValue(Y_WALL, 100.0 * y_wall_val)

    if(number_of_exits == 0 and len(out_skin_list) != 0):
        max_d = 0.0
        for node in out_skin_list:
            d = node.GetSolutionStepValue(DISTANCE)
            if(d > max_d):
                max_d = d

        dlimit = 0.5 * max_d

        for node in out_skin_list:
            d = node.GetSolutionStepValue(DISTANCE)
            if(d > dlimit):
                node.Fix(PRESSURE)
                node.SetSolutionStepValue(PRESSURE, 0, 0.0)
                node.SetValue(IS_STRUCTURE, 0.0)
                node.SetSolutionStepValue(IS_FREE_SURFACE, 0, 20.0)

        #farthest_node = 0
        #max_d = 0.0
        # for node in out_skin_list:
            #d = node.GetSolutionStepValue(DISTANCE)
            # if(d > max_d):
                #max_d = d
                #farthest_node = node

        # farthest_node.Fix(PRESSURE)
        #farthest_node.SetSolutionStepValue(PRESSURE, 0, 0.0)
        #farthest_node.SetValue(IS_STRUCTURE, 0.0)


# elif(SolverType == "monolithic_solver_eulerian"):
#fluid_solver = Kratos_dpg_monolithic_solver_eulerian.MonolithicSolver(fluid_model_part,domain_size)
#fluid_solver.oss_switch =ProjectParameters.use_orthogonal_subscales
#fluid_solver.dynamic_tau_fluid =ProjectParameters.use_dt_in_stabilization
fluid_solver.dynamic_tau_levelset = 0.001
#fluid_solver.max_iter =  ProjectParameters.max_iterations
fluid_solver.echo_level = 0
#fluid_solver.CalculateReactionFlag = ProjectParameters.Calculate_reactions
fluid_solver.ReformDofSetAtEachStep = False
#fluid_solver.rel_vel_tol = ProjectParameters.velocity_relative_tolerance
#fluid_solver.abs_vel_tol = ProjectParameters.velocity_absolute_tolerance
#fluid_solver.rel_pres_tol = ProjectParameters.pressure_relative_tolerance
#fluid_solver.abs_pres_tol = ProjectParameters.pressure_absolute_tolerance
#fluid_solver.linear_solver =  monolithic_linear_solver
fluid_solver.use_slip_conditions = True
fluid_solver.volume_correction_switch = volume_correction_switch
fluid_solver.redistance_frequency = 1
fluid_solver.CFL = CFL
fluid_solver.rho1 = 1000
fluid_solver.rho2 = 1
fluid_solver.mu = 3.0e-3
fluid_solver.Initialize()


if(fluid_model_part.NumberOfMeshes() > 1):
    for mesh_number in range(2, fluid_model_part.NumberOfMeshes()):
        #mesh_nodes = model_part.GetMesh(mesh_number).Nodes
        if(fluid_model_part.Properties[mesh_number][IS_SLIP] == 0.0):
            mesh_nodes = fluid_model_part.GetMesh(mesh_number).Nodes
            is_calc = 0
            for node in mesh_nodes:
                node.Fix(VELOCITY_X)
                node.Fix(VELOCITY_Y)
                node.Fix(VELOCITY_Z)
                node.SetValue(IS_STRUCTURE, 0)
                node.SetSolutionStepValue(IS_SLIP, 0, 0.0)
                node.SetSolutionStepValue(DISTANCE, 0, -100.0)
                node.Fix(DISTANCE)
                node.SetSolutionStepValue(FLAG_VARIABLE, 0, 11.0)
                if(is_calc == 0):
                    vx = node.GetSolutionStepValue(VELOCITY_X)
                    vy = node.GetSolutionStepValue(VELOCITY_Y)
                    vz = node.GetSolutionStepValue(VELOCITY_Z)
                    inlet_velocity = vx * vx + vy * vy + vz * vz
                    inlet_velocity = math.sqrt(inlet_velocity)
                    is_calc = 1

        if(fluid_model_part.Properties[mesh_number][IS_SLIP] == 1.0):
            mesh_nodes = fluid_model_part.GetMesh(mesh_number).Nodes
            for node in mesh_nodes:
                node.SetValue(IS_STRUCTURE, 1)
                node.SetSolutionStepValue(FLAG_VARIABLE, 0, 11.0)
    # assign free IS_SLIP different than 10 and 20 and 30 to recognize free BC
    for condition in fluid_model_part.Conditions:
        for node in condition.GetNodes():
            flg = node.GetSolutionStepValue(FLAG_VARIABLE)
            if(flg != 11.0):
                node.SetSolutionStepValue(IS_SLIP, 0, 40.0)


# Stepping and time settings
Dt = ProjectParameters.Dt
max_Dt = 2.0 * Dt
Nsteps = ProjectParameters.nsteps
final_time = ProjectParameters.max_time
output_time = ProjectParameters.output_time
time = 0.0  # ProjectParameters.Start_time
out = 0
step = 0

# volume correction
net_volume = 0.0
fluid_model_part.ProcessInfo.SetValue(NET_INPUT_MATERIAL, 0.0)


# mesh to be printed
mesh_name = 0.0
gid_io.InitializeMesh(mesh_name)
gid_io.WriteMesh(fluid_model_part.GetMesh())
gid_io.FinalizeMesh()
gid_io.Flush()
gid_io.InitializeResults(mesh_name, (fluid_model_part).GetMesh())


time = time + max_Dt * 0.0001
fluid_model_part.CloneTimeStep(time)


# print "fluid solver created"
while(time <= final_time):
    Dt = EstimateTimeStep3D().ComputeDt(fluid_model_part, 2.0 * fluid_solver.max_edge_size, CFL, 0.03 * max_Dt, max_Dt)

    time = time + Dt

    fluid_model_part.CloneTimeStep(time)

    print("STEP = ", step)
    print("TIME = ", time)
    #elapsed_time = timer.clock() - start_time
    elapsed_time = time_measure.time() - start_time
    print("DT= ", Dt, "max_dt= ", max_Dt, "Total time= ", elapsed_time)

    BiphasicFillingUtilities().ComputeNetInletVolume(fluid_model_part)
    OpenAirExitInFarDryZone(fluid_model_part, ProjectParameters.body_force_x, ProjectParameters.body_force_y, ProjectParameters.body_force_z)
    fluid_solver.Solve(step)

    max_acceptable_acc_norm = 6.0 * inlet_velocity / Dt
    BiphasicFillingUtilities().ApplyVelocityLimitation(fluid_model_part, max_acceptable_acc_norm)

    # if(out == output_step):
    if(output_time <= out or 1 == 1):
        gid_io.WriteNodalResults(PRESSURE, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(VELOCITY, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(DISTANCE, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(BODY_FORCE, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(IS_SLIP, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(FLAG_VARIABLE, fluid_model_part.Nodes, time, 0)

        last_print = time
        gid_io.Flush()
        sys.stdout.flush()
        out = 0

    out = out + Dt
    step = step + 1
# write the last step results
gid_io.FinalizeResults()
