domain_size = 3

from KratosMultiphysics import *
from KratosMultiphysics.BloodFlowApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.MeshingApplication import *


def GetNodeAfter(table, prop):
    return table[prop - 1][4]


import math
import time
import full_nodes_table

# defining a model part for the fluid and one for the structure
model_part1D = ModelPart("FluidPart")
model_part3D = ModelPart("FluidPart3D")

#proc_info = model_part1D.ProcessInfo

# defining a model part for the fluid and one for the structure
model_part1D.AddNodalSolutionStepVariable(RADIUS)
model_part1D.AddNodalSolutionStepVariable(BETA)
model_part1D.AddNodalSolutionStepVariable(C0)
model_part1D.AddNodalSolutionStepVariable(NODAL_AREA)
model_part1D.AddNodalSolutionStepVariable(NODAL_MASS)
model_part1D.AddNodalSolutionStepVariable(RHS)
model_part1D.AddNodalSolutionStepVariable(WORK)
model_part1D.AddNodalSolutionStepVariable(FLOW)
model_part1D.AddNodalSolutionStepVariable(VELOCITY)
model_part1D.AddNodalSolutionStepVariable(THICKNESS)
model_part1D.AddNodalSolutionStepVariable(YOUNG_MODULUS)
model_part1D.AddNodalSolutionStepVariable(POISSON_RATIO)
model_part1D.AddNodalSolutionStepVariable(DENSITY)
model_part1D.AddNodalSolutionStepVariable(TERMINAL_RESISTANCE)
model_part1D.AddNodalSolutionStepVariable(FLAG_VARIABLE)
model_part1D.AddNodalSolutionStepVariable(PRESSURE)

import vms_fractional_step_solver as solver

solver.AddVariables(model_part3D)
model_part3D.AddNodalSolutionStepVariable(FLAG_VARIABLE)

# ARCHIVE TO SET ::::::::::::::::::::::::::: >>>>> VARIABLES
import config
import removal_tool

var_aux = True
ascii = config.ascii_results
relative_path_3D = config.relative_path_3D
relative_path_1D = config.relative_path_1D
name_model_3D_2 = config.name
cardiac_cycle = config.nro_cardiac_cycles
time_cardiac_cycle = 1.63 # total_time (last value of the Cardiac_cycle
initial_pressure = 0  # Pa

Coupled_Simulation = True  # True-->Coupled_3d_1d False-->Only 1D
Sub_steping = False  # True-->Sub_step_control (only for the Coupled_3d_1d)
sub_step = 200 #config.sub_step
step_size_control = False  # True-->Activate	False-->fix (step_size)
step_size = 0.00001 #config.step_size

CardiacCycleConvergence = False

artery_type = config.artery_type[0]
cardiac_cycle_to_3D = cardiac_cycle  # 3D cardiac_Cycle will be running
final_time = cardiac_cycle * time_cardiac_cycle

name_model_3D = relative_path_3D + name_model_3D_2
print name_model_3D

if (cardiac_cycle <= 0):
    print "Number of Cardiac Cycles must be > 0"
    print "Please check config file"
    raw_input()
    error

if (cardiac_cycle_to_3D > cardiac_cycle):
    print "Please revise:: cardiac_cycle_to_3D"
    print "Please check config file"
    raw_input()
    error

if (artery_type == 1):
    input_file_name = "Left_Balanced_Dominant"
    input_file_name = relative_path_1D + "Left_Balanced_Dominant"
    print "1D Arterial model:::>>> ", input_file_name

if (artery_type == 2):
    input_file_name = "Left_LCA_Dominant"
    input_file_name = relative_path_1D + "Left_LCA_Dominant"
    print "1D Arterial model:::>>> ", input_file_name

if (artery_type == 3):
    input_file_name = "Left_RCA_Dominant"
    input_file_name = relative_path_1D + "Left_RCA_Dominant"
    print "1D Arterial model:::>>> ", input_file_name

if (artery_type == 4):
    input_file_name = "Left_Small_RCA_Dominant"
    input_file_name = relative_path_1D + "Left_Small_RCA_Dominant"
    print "1D Arterial model:::>>> ", input_file_name

if (artery_type == 5):
    input_file_name = "Right_Balanced_Dominant"
    input_file_name = relative_path_1D + "Right_Balanced_Dominant"
    #input_file_name = relative_path_1D + "alastruey_new_coronaries"
    print "1D Arterial model:::>>> ", input_file_name

if (artery_type == 6):
    input_file_name = "Right_LCA_Dominant"
    input_file_name = relative_path_1D + "Right_LCA_Dominant"
    print "1D Arterial model:::>>> ", input_file_name

if (artery_type == 7):
    input_file_name = "Right_RCA_Dominant"
    input_file_name = relative_path_1D + "Right_RCA_Dominant"
    print "1D Arterial model:::>>> ", input_file_name

if (artery_type == 8):
    input_file_name = "Right_Small_RCA_Dominant"
    input_file_name = relative_path_1D + "Right_Small_RCA_Dominant"
    print "1D Arterial model:::>>> ", input_file_name

# reading the fluid part
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteConditions
# WriteElementsOnly
# WriteConditions
model_part_io_fluid = ModelPartIO(input_file_name)
model_part_io_fluid.ReadModelPart(model_part1D)

# reading the fluid part
#model_part_io_3D = ModelPartIO("volume_mesh")
# Only when 1D problem is running, please Verify that you are saveing the
# results-->Variable StepCoupled
if (Coupled_Simulation == False):
    if (ascii == False):
        gid_io = GidIO(input_file_name, gid_mode, multifile, deformed_mesh_flag, write_conditions)
    print "Only 1D Problem is running"
    plot3d = False  # Plot 1D variables
    Sub_steping = False
    output_step = 1  # for each time step, ie: 1 save 1 Dt, 100 save 100 Dt
else:  # 3D-1D coupled problem is running
    # Plot3d=true---> Plot 3D variables.
    # Plot3d=False---> Plot 1D variables.
    plot3d = True
    if (plot3d == True):
        if (ascii == False):
            gid_io = GidIO(name_model_3D, gid_mode, multifile, deformed_mesh_flag, write_conditions)
    else:
        if (ascii == False):
            gid_io = GidIO(input_file_name, gid_mode, multifile, deformed_mesh_flag, write_conditions)
    print "Coupled 3D-1D Problem is running"
    removal_tool.DoRemoval(model_part1D)
    model_part_io_3D = ModelPartIO(name_model_3D)
    model_part_io_3D.ReadModelPart(model_part3D)
    model_part3D.SetBufferSize(3)
    solver.AddDofs(model_part3D)

integrator = ArteryTimeIntegrator()
integrator.Initialize(model_part1D)
inletconditiontable = model_part1D.GetTable(1)

minlength = 1e+12
minlength = integrator.Element_minLength(model_part1D)

# setting up the buffer size: SHOULD BE DONE AFTER READING!!!
model_part1D.SetBufferSize(2)
out = 1
# mesh to be printed

if (ascii == False):
    if (plot3d == False):
        mesh_name = 0.0
        gid_io.InitializeMesh(mesh_name)
        gid_io.WriteMesh(model_part1D.GetMesh())
        gid_io.WriteNodeMesh(model_part1D.GetMesh())
        gid_io.FinalizeMesh()
        gid_io.InitializeResults(mesh_name, (model_part1D).GetMesh())
        print "Writing 1D Mesh------------------------"
    else:
        mesh_name = 0.0
        gid_io.InitializeMesh(mesh_name)
        gid_io.WriteMesh(model_part3D.GetMesh())
        gid_io.FinalizeMesh()
        gid_io.InitializeResults(mesh_name, (model_part3D).GetMesh())
        print "Writing 3D Mesh------------------------"

# Initial Flow
q_inicial = inletconditiontable.GetValue(0)
# time_cardiac_cycle=inletconditiontable.GetValue()
print "Total Time:   ", time_cardiac_cycle

# Initial Values
for node in model_part1D.Nodes:
    if (node.IsFixed(FLOW) == False):
        node.SetSolutionStepValue(FLOW, 0, q_inicial)
        node.SetSolutionStepValue(PRESSURE, 0, initial_pressure)
    else:
        node.SetSolutionStepValue(FLOW, 0, q_inicial)
        node.SetSolutionStepValue(PRESSURE, 0, initial_pressure)

if (Coupled_Simulation == True):
    solver3D = solver.IncompressibleFluidSolver(model_part3D, 3)
    solver3D.max_val_its = 2
    solver3D.max_press_its = 2
    solver3D.predictor_corrector = False
    solver3D.vel_toll = 1e-3
    solver3D.press_toll = 1e-3
    solver3D.dynamic_tau = 1.0
    solver3D.compute_reactions = False
    solver3D.Initialize()
    import CouplingTools1D_3Dv3
    transfer_obj = CouplingTools1D_3Dv3.TransferTools(model_part1D, model_part3D)
    transfer_obj.Setting3d()
    transfer_obj.Initialize()
    #transfer_obj.Initial_Contitions()
    print "READY TO SIMULATE................................."
else:
    print "ONLY 1D SOLVER...................................."

time = 0.0
total_time = 0.0
step = 0.0
minlength = 1e+12
minlength = integrator.Element_minLength(model_part1D)
nro_cardiac_cycle = 1.0
inicial = 1
results = str(name_model_3D + ".cvpr")
f = open(results, 'w')

fixed_flow_nodes = []
for node in model_part1D.Nodes:
    if (node.IsFixed(FLOW) == True and node.GetSolutionStepValue(FLAG_VARIABLE) == 0):
        fixed_flow_nodes.append(node)
#import time as timer
#timer.sleep(1)
#import sys
#sys.stdout.flush()

total_time = total_time + step_size
##Flag_Converge=False

#Initial Delta Step
Dt = integrator.EstimateDeltaTime(model_part1D, 0.8, minlength)

# Aprox to save results
steps_eval = final_time / Dt
save_results = config.save_results
save_Steps = final_time / save_results
if (save_results < steps_eval):
    output_step = math.floor(steps_eval / save_Steps)
else:
    out_res = 1
    output_step = math.floor(steps_eval / out_res)

if (Sub_steping == True):
    control_sub_step = 1
    #sub_step = 200
    if (output_step < sub_step):
        output_step = sub_step
        #for each time step, ie: 1 save 1 Dt, 100 save 100 Dt
else:
    Sub_steping == True
    control_sub_step = 1
    #sub_step = 200
    if ((output_step < sub_step) and (Coupled_Simulation == True)):
        output_step = sub_step

while (total_time < final_time):

    model_part1D.CloneTimeStep(total_time)

    if (step < 3):
        model_part3D.CloneTimeStep(total_time)

    #for node in model_part1D.Nodes:
    for node in fixed_flow_nodes:
        if (Coupled_Simulation == False):
            if (node.IsFixed(FLOW) == True and node.GetSolutionStepValue(FLAG_VARIABLE) == 0):
                if (step > 1):
                    Q = inletconditiontable.GetValue(time)
                    Q = 2 * Q - node.GetSolutionStepValue(FLOW, 1)
                    node.SetSolutionStepValue(FLOW, 0, Q)
                else:
                    Q = inletconditiontable.GetValue(time)
                    node.SetSolutionStepValue(FLOW, 0, Q)
        else:
            if ((node.IsFixed(FLOW) == True) and node.GetSolutionStepValue(FLAG_VARIABLE) == 0):
                if (step > 1):
                    Q = inletconditiontable.GetValue(time)
                    Q = 2 * Q - node.GetSolutionStepValue(FLOW, 1)
                    node.SetSolutionStepValue(FLOW, 0, Q)
                else:
                    Q = inletconditiontable.GetValue(time)
                    node.SetSolutionStepValue(FLOW, 0, Q)

    #Q = inletconditiontable.GetValue(time)
    #Q = 2 * Q - node.GetSolutionStepValue(FLOW, 1)
    #for node in fixed_flow_nodes:
    #    if (step > 1):
    #        Q = inletconditiontable.GetValue(time)
    #        Q = 2 * Q - node.GetSolutionStepValue(FLOW, 1)
    #        node.SetSolutionStepValue(FLOW, 0, Q)
    #    else:
    #        Q = inletconditiontable.GetValue(time)
    #        node.SetSolutionStepValue(FLOW, 0, Q)

    if ((Coupled_Simulation == True) and (control_sub_step == sub_step) and (nro_cardiac_cycle == cardiac_cycle_to_3D)):
        if (step >= 3):
            print "...............................Solve 3D.................................."
            # integrator.SolveStep(model_part1D)
            # integrator.ComputePressure(model_part1D)
            model_part3D.CloneTimeStep(total_time)
            print "Total_time 3d:", total_time
            print "----------------------------Transfer 1d(tn) to 3d(tn+1)------------------"
            transfer_obj.Transfer1D_to_3D()
            print "Solve 3D para ", total_time
            solver3D.Solve()
            print "----------------------------Transfer 3d to 1d----------------------------"
            transfer_obj.Transfer3D_to_1D()
            integrator.SolveStep(model_part1D)
            integrator.ComputePressure(model_part1D)
            # removal_tool.ComputePressure(model_part1D)
            control_sub_step = 0
            var_aux = True
            out = out + 1
    else:
        if ((var_aux == True) and (Sub_steping == True)):
            print "Solve 1D ------------------------------> ", total_time
            print "Sub_step-->", sub_step
            var_aux = False
        else:
            print "Solve 1D ------------------------------> ", total_time

        integrator.SolveStep(model_part1D)
        # removal_tool.ComputePressure(model_part1D) --> This is only needed when the results are saveed.
        control_sub_step = control_sub_step + 1
        out = out + 1
        ##print "Step--->End"

    print math.fmod(step, output_step)
    #if((control_sub_step == output_step)):  # or (Coupled_Simulation == False)):
    if (math.fmod(step, output_step) == 0):
        if (plot3d == False):
        # removal_tool.ComputePressure(model_part1D)
            integrator.ComputePressure(model_part1D) # Only when I saved the results.
        if (CardiacCycleConvergence == True):
            print "Check Convergence"
            # Function to check the cardiac_Cycle convergence (move outside)
        time_cardiac_cycle = integrator.CheckCardiacCovergence(model_part1D, time_cardiac_cycle)

        if (ascii == True):
            f.write("Time_step \n")
            f.write(str(total_time))
            f.write("\n")
            # f.write("Number_of_nodes\n")
            f.write(
                "Nodal_data_table (Index-NodeId-PRESSURE-VELOCITY_X-RADIUS) \n")
            indextowrite = 1
            for node in model_part1D.Nodes:
                nodewrite = node.Id
                ToWrite = str(indextowrite) + " "
                ToWrite += str(nodewrite) + " " + str(
                    node.GetSolutionStepValue(PRESSURE)) + " "
                ToWrite += str(node.GetSolutionStepValue(VELOCITY_X)) + " " + str(
                    node.GetSolutionStepValue(RADIUS)) + "\n"
                indextowrite = indextowrite + 1
                f.write(ToWrite)
            print "Writing CVPR results 1D Results(time)------------------------>", total_time
        else:
            gid_io.WriteNodalResults(NODAL_AREA, model_part1D.Nodes, total_time, 0)
            gid_io.WriteNodalResults(NODAL_MASS, model_part1D.Nodes, total_time, 0)
            gid_io.WriteNodalResults(FLOW, model_part1D.Nodes, total_time, 0)
            gid_io.WriteNodalResults(RADIUS, model_part1D.Nodes, total_time, 0)
            gid_io.WriteNodalResults(RHS, model_part1D.Nodes, total_time, 0)
            gid_io.WriteNodalResults(VELOCITY, model_part1D.Nodes, total_time, 0)
            gid_io.WriteNodalResults(PRESSURE, model_part1D.Nodes, total_time, 0)
            gid_io.WriteNodalResults(WORK, model_part1D.Nodes, total_time, 0)
            print "Writing GiD results 1D Results(time)------------------------>", total_time
        if (Coupled_Simulation == False):
            control_sub_step = 0
    else:
        ##node_out=GetNodeAfter(full_nodes_table.table, config.inlets_1d[0][1])
        # print node_out
        # press_1d_out=model_part1D.Nodes[node_out].GetSolutionStepValue(PRESSURE)
        if (ascii == True):
            f.write("Time_step \n")
            f.write(str(total_time))
            f.write("\n")
            # f.write("Number_of_nodes\n")
            # f.write(xxxx)
            f.write(
                "Nodal_data_table (Index-NodeId-Pressure-Velocity_X-Velocity_Y-Velocity_Z) \n")
            indextowrite = 1
            for node in model_part3D.Nodes:
                nodewrite = node.Id
                pressure_rect = node.GetSolutionStepValue(PRESSURE)
                # +press_1d_out
                ToWrite = str(indextowrite) + " " + str(
                    nodewrite) + " " + str(pressure_rect) + " "
                ToWrite += str(node.GetSolutionStepValue(VELOCITY_X)) + " "
                ToWrite += str(node.GetSolutionStepValue(VELOCITY_Y)) + " "
                ToWrite += str(
                    node.GetSolutionStepValue(VELOCITY_Z)) + "\n"
                indextowrite = indextowrite + 1
                f.write(ToWrite)
            print "Writing CVPR results 1D Results(time)------------------------>", total_time
        else:
            gid_io.WriteNodalResults(VELOCITY, model_part3D.Nodes, total_time, 0)
            gid_io.WriteNodalResults(PRESSURE, model_part3D.Nodes, total_time, 0)
            gid_io.WriteNodalResults(VISCOSITY, model_part3D.Nodes, total_time, 0)
            gid_io.WriteNodalResults(DENSITY, model_part3D.Nodes, total_time, 0)
            gid_io.WriteNodalResults(FLAG_VARIABLE, model_part3D.Nodes, total_time, 0)
            out = 0
            print "Writing GiD results 3D results(step)------------------------>", total_time

    if (time > time_cardiac_cycle):
        time = time - time_cardiac_cycle
        nro_cardiac_cycle = nro_cardiac_cycle + 1
        control_sub_step = 0
        out = 0

    if (step_size_control == True):
        Dt = integrator.EstimateDeltaTime(model_part1D, 0.8, minlength)
        #print "Using Variable_Delta_Step:", Dt
    else:
        Dt = step_size

    time = time + Dt
    total_time = total_time + Dt
    out = out + 1
    step = step + 1
    # print "out_step", output_step
    # print "StepCoupled", StepCoupled
    # print "out", out

print "Time=", total_time
print "---------------------------------------------END-------------------------------------------"
gid_io.FinalizeResults()
