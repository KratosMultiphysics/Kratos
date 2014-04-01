import sys
output_file=open("python_output","w")
sys.stdout =output_file

domain_size = 3

from KratosMultiphysics import *
from KratosMultiphysics.BloodFlowApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.ExternalSolversApplication import *


def GetNodeAfter(table, prop):
    return table[prop - 1][4]

import math
import time
import Only1D


# defining a model part for the fluid and one for the structure
model_part1D = ModelPart("FluidPart")
model_part3D = ModelPart("FluidPart3D")

# proc_info = model_part1D.ProcessInfo

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
model_part1D.AddNodalSolutionStepVariable(PRESSURE_VENOUS)
model_part1D.AddNodalSolutionStepVariable(SYSTOLIC_PRESSURE)
model_part1D.AddNodalSolutionStepVariable(DYASTOLIC_PRESSURE)
model_part1D.AddNodalSolutionStepVariable(AVERAGE_PRESSURE)

import vms_fractional_step_solver as solver
solver.AddVariables(model_part3D)
model_part3D.AddNodalSolutionStepVariable(FLAG_VARIABLE)

# ARCHIVE TO SET ::::::::::::::::::::::::::: >>>>> VARIABLES
import config
import removal_tool

print("CONFIG TO VERSION 20_Noviembre_2013")


FitValues = True  # To get A-B of the 3D model
var_aux = True
only1Dtest = config.only1Dtest
Fit_control = False
FitRadius = config.FitRadius
ascii = config.ascii_results
relative_path_3D = config.relative_path_3D
relative_path_1D = config.relative_path_1D
name_model_3D_2 = config.name
cardiac_cycle = config.nro_cardiac_cycles
# total_time (last value of the Cardiac_cycle
time_cardiac_cycle = config.time_period

diastolic_pressure = config.diastolic_pressure  # Pa
systolic_pressure = config.systolic_pressure  # Pa
time_period = config.time_period

Coupled_Simulation = True  # True-->Coupled_3d_1d False-->Only 1D
Sub_steping = True  # True-->Sub_step_control (only for the Coupled_3d_1d)
sub_step = 50  # config.sub_step
step_size_control = True  # True-->Activate	False-->fix (step_size)
step_size = 0.00001  # config.step_size
CardiacCycleConvergence = False
nro_cardiac_cycle = 1.0
pressure_factor = config.pressure_factor

artery_type = config.artery_type[0]
cardiac_cycle_to_3D = cardiac_cycle  # 3D cardiac_Cycle will be running
final_time = cardiac_cycle * time_cardiac_cycle


diastolic_pressure = pressure_factor * diastolic_pressure
systolic_pressure = pressure_factor * systolic_pressure

InletProfileType = config.inlet_pressure_type
if ((InletProfileType == "coseno") or (InletProfileType == "parabolic") or (InletProfileType == "table")):
    InletPressure = True
else:
    InletPressure = False
# print "Coupled_Simulation", Coupled_Simulation
# print "FitValues", FitValues
# print "Fit_control", Fit_control

name_model_3D = relative_path_3D + name_model_3D_2
print("----------------------------------------------------------------------------------------------------------")
print("Coupled_Simulation ", Coupled_Simulation)
print("name_model_3D ", name_model_3D)
print("Total Time:   ", time_cardiac_cycle)

if (cardiac_cycle <= 0):
    print("Number of Cardiac Cycles must be > 0")
    print("Please check config file")
    input()
    error

if (cardiac_cycle_to_3D > cardiac_cycle):
    print("Please revise:: cardiac_cycle_to_3D")
    print("Please check config file")
    input()
    error

if (systolic_pressure < diastolic_pressure):
    print("Please revise:: systolic and diastolic pressure:")
    print("systolic pressure -> ", systolic_pressure, "diastolic_pressure-> ", diastolic_pressure)
    print("Please check config file")
    input()
    error

if (artery_type == 1):
    input_file_name1 = "Left_Balanced_Dominant"
    input_file_name = relative_path_1D + "Left_Balanced_Dominant"
    print("1D Arterial model:::>>> ", input_file_name)

if (artery_type == 2):
    input_file_name1 = "Left_LCA_Dominant"
    input_file_name = relative_path_1D + "Left_LCA_Dominant"  # _pressure"
    print("1D Arterial model:::>>> ", input_file_name)

if (artery_type == 3):
    input_file_name1 = "Left_RCA_Dominant"
    input_file_name = relative_path_1D + "Left_RCA_Dominant"
    print("1D Arterial model:::>>> ", input_file_name)

if (artery_type == 4):
    input_file_name1 = "Left_Small_RCA_Dominant"
    input_file_name = relative_path_1D + "Left_Small_RCA_Dominant"
    print("1D Arterial model:::>>> ", input_file_name)

if (artery_type == 5):
    input_file_name1 = "Right_Balanced_Dominant"
    input_file_name = relative_path_1D + "Right_Balanced_Dominant"
    print("1D Arterial model:::>>> ", input_file_name)

if (artery_type == 6):
    input_file_name1 = "Right_LCA_Dominant"
    input_file_name = relative_path_1D + "Right_LCA_Dominant"
    print("1D Arterial model:::>>> ", input_file_name)

if (artery_type == 7):
    input_file_name1 = "Right_RCA_Dominant"
    input_file_name = relative_path_1D + "Right_RCA_Dominant"
    print("FitValues1D Arterial model:::>>> ", input_file_name)

if (artery_type == 8):
    input_file_name1 = "Right_Small_RCA_Dominant"
    input_file_name = relative_path_1D + "Right_Small_RCA_Dominant"
    print("1D Arterial model:::>>> ", input_file_name)

if (artery_type == 9):
    input_file_name1 = "Right_test"
    input_file_name = relative_path_1D + "Right_test"
    print("1D Arterial model:::>>> ", input_file_name)

if ((artery_type < 5)):
    for i in range(0, len(config.deactivate_list)):
        artery_number_id = config.deactivate_list[i]
        print(artery_number_id)
        if (artery_number_id < 15):
            print("Please check your config file. Artery number = ", artery_number_id, ", dont belong to the ", input_file_name1)
            print("Please check your model")
            input()
            error
else:
    for i in range(0, len(config.deactivate_list)):
        artery_number_id = config.deactivate_list[i]
        if (artery_number_id > 16):
            print("Please check your config file. Artery number = ", artery_number_id, "dont belong to the ", input_file_name1)
            print("Please check your model")
            input()
            error

# AB_Matriz=[]
# row_AB=0
# for i in range(3): #nodes_number
    # AB_Matriz.append([0]*4)   #Variables to save

# ffit_A=[]
# ffit_A.append(0)
#results = str("FitValues_AB.txt")
#ffit_A[0] = open(results, 'w')

#ffit = []
# ffit.append(0)
# ffit.append(1)
# ffit.append(2)

# results = str("wqeqweqwewq.txt")3DQCA059_11
#ffit[0] = open(results, 'w')
#ffit[0].write("NODE-PRESSURE_1D-FLOW1D-NODAL_AREA_1D) \n")
#results = str("FitValues_2.txt")
#ffit[1] = open(results, 'w')
#ffit[1].write("NODE-PRESSURE_1D-FLOW1D-NODAL_AREA_1D) \n")
#results = str("FitValues_3.txt")
#ffit[2] = open(results, 'w')
#ffit[2].write("NODE-PRESSURE_1D-FLOW1D-NODAL_AREA_1D) \n")

#ffit3d =[]
# ffit3d.append(0)
#results = str("FitValues_3d_1.txt")
#ffit3d[0] = open(results, 'w')

# raw_input()
# reading the fluid part

gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
# WriteElementsOnly
# WriteConditions
model_part_io_fluid = ModelPartIO(input_file_name)
model_part_io_fluid.ReadModelPart(model_part1D)

# for prop in model_part1D.Properties:
        # R=prop.GetValue(RADIUS)
        # print R

# raw_input()
# reading the fluid part
# model_part_io_3D = ModelPartIO("volume_mesh")
# Only when 1D problem is running, please Verify that you are saveing the
# results-->Variable StepCouplednt

# Fit_control = This variable is only for doing test with the 1D model.
# CouplingTools1D_3Dv5 is needed to inizializate the radius
if((FitValues == True) and (Coupled_Simulation == False)):
    Coupled_Simulation = True
    Fit_control = True
    print("Fit_control", Fit_control)

plot1d = True
plot3d = True
if (Coupled_Simulation == False):
    if (ascii == False):
        gid_io = GidIO(input_file_name, gid_mode,
                       multifile, deformed_mesh_flag, write_conditions)
    else:
        results = str(input_file_name + ".cvpr")
        f1d = open(results, 'w')
    print("Only 1D Problem is running")
    plot3d = False  # Plot 1D variables
    Sub_steping = False
    output_step = 1  # for eachffit.append = 1
else:  # 3D-1D coupled problem is running
    # Plot3d=true---> Plot 3D variables. Plot3d=False---> Plot 1D variables.
    if (plot3d == True):
        if(ascii == False):
            gid_io_3d = GidIO(
                name_model_3D, gid_mode, multifile, deformed_mesh_flag, write_conditions)
        else:
            results = str(name_model_3D + ".cvpr")
            f3d = open(resuls, 'w')
    if (plot1d == True):
        if (ascii == False):
            gid_io = GidIO(
                input_file_name, gid_mode, multifile, deformed_mesh_flag, write_conditions)
        else:
            results = str(input_file_name + ".cvpr")
            f1d = open(results, 'w')
    print("Coupled 3D-1D Problem is running")
    #inlets_outles_nodes = []
    FFR_Inlet_NODES_Values = []
    FFR_Outlet_NODES_Values = []
    [FFR_Inlet_NODES_Values,
        FFR_Outlet_NODES_Values] = removal_tool.DoRemoval(model_part1D)
    # for i in range(0, len(inlets_outles_nodes[0])):
        # node_inlet=inlets_outles_nodes[i]
        # FFR_Inlet_NODES_Values.append(node_inlet)
    # for i in range(0, len(inlets_outles_nodes[1])):
        #node_outlet = inlets_outles_nodes[i]
        # FFR_Outlet_NODES_Values.append(node_outlet)

    # for j in range(0, len(FFR_Inlet_NODES_Values)):
        # print "FFR_Inlet_NODES_Values", FFR_Inlet_NODES_Values[j].Id
    # for j in range(0, len(FFR_Outlet_NODES_Values)):
        # print "FFR_Outlet_NODES_Values", FFR_Outlet_NODES_Values[j].Id
    # raw_input()
    model_part_io_3D = ModelPartIO(name_model_3D)
    model_part_io_3D.ReadModelPart(model_part3D)
    model_part3D.SetBufferSize(3)
    solver.AddDofs(model_part3D)

# setting up the buffer size: SHOULD BE DONE AFTER READING!!!
model_part1D.SetBufferSize(2)
out = 1
# mesh to be printed
if(ascii == False):
    if(plot1d == True):
        mesh_name = 0.0
        gid_io.InitializeMesh(mesh_name)
        gid_io.WriteMesh(model_part1D.GetMesh())
        gid_io.WriteNodeMesh(model_part1D.GetMesh())
        gid_io.FinalizeMesh()
        gid_io.InitializeResults(mesh_name, (model_part1D).GetMesh())
        print("Writing 1D Mesh------------------------")
    if(plot3d == True and only1Dtest == False):
        mesh_name = 0.0
        gid_io_3d.InitializeMesh(mesh_name)
        gid_io_3d.WriteMesh(model_part3D.GetMesh())
        gid_io_3d.FinalizeMesh()
        gid_io_3d.InitializeResults(mesh_name, (model_part3D).GetMesh())
        print("Writing 3D Mesh------------------------")

# Initial conditions for the 1D model.
integrator = ArteryTimeIntegrator()
integrator.Initialize(model_part1D)
inletconditiontable = model_part1D.GetTable(1)
minlength = 1e+12
minlength = integrator.Element_minLength(model_part1D)
# Initial Flow
#q_inicial = inletconditiontable.GetValue(0)
# print "Q=====>", q_inicial
# time_cardiac_cycle=inletconditiontable.GetValue()
for cond in model_part1D.Conditions:
    for node in cond.GetNodes():
    # cond.SetValue(PRESSURE_VENOUS,dyastolic_hypermia_pressure)
    #node.SetSolutionStepValue(PRESSURE_VENOUS, 0, dyastolic_hypermia_pressure)
        Resistence_factor = 1
        cond.SetValue(PRESSURE_DT, Resistence_factor)
        # print diastolic_pressure
        # print PRESSURE_VENOUS
        # print node.Id
        # print cond.Id
# Initial Values. Set the intial pressure as reference for the 1D model
# (all nodes take the systolic pressure and the initial flow)
if (InletPressure == True):
    if (InletProfileType == "parabolic"):
        print("Using Parabolic Pressure as input:")
        q_inicial = config.Q_initial
        pressure_parameter_1 = (
            (systolic_pressure - diastolic_pressure) * 4) / time_period
        pressure_parameter_2 = (-pressure_parameter_1) / time_period
        pressure_parameter_3 = diastolic_pressure
    elif (InletProfileType == "coseno"):
        q_inicial = config.Q_initial
        pressure_parameter_1 = (systolic_pressure + diastolic_pressure) / 2
        pressure_parameter_2 = (diastolic_pressure - systolic_pressure) / 2
    print("Pressure_Factor: ", pressure_factor)
    print("Diastolic pressure: ", diastolic_pressure)
    print("Systolic Pressure: ", systolic_pressure)
    # diastolic_pressure=inletconditiontable.GetValue(0)
    for node in model_part1D.Nodes:
        if(node.IsFixed(PRESSURE) == False):
            node.SetSolutionStepValue(FLOW, 0, q_inicial)
            node.SetSolutionStepValue(PRESSURE, 0, diastolic_pressure)
            node.SetSolutionStepValue(
                DYASTOLIC_PRESSURE, 0, diastolic_pressure)
            node.SetSolutionStepValue(PRESSURE_VENOUS, 0, diastolic_pressure)
        else:
            node.SetSolutionStepValue(FLOW, 0, q_inicial)
            node.SetSolutionStepValue(PRESSURE, 0, diastolic_pressure)
            node.SetSolutionStepValue(
                DYASTOLIC_PRESSURE, 0, diastolic_pressure)
            node.SetSolutionStepValue(PRESSURE_VENOUS, 0, diastolic_pressure)
else:
    q_inicial = inletconditiontable.GetValue(0)
    for node in model_part1D.Nodes:
        if(node.IsFixed(FLOW) == False):
            node.SetSolutionStepValue(FLOW, 0, q_inicial)
            node.SetSolutionStepValue(PRESSURE, 0, diastolic_pressure)
            node.SetSolutionStepValue(
                DYASTOLIC_PRESSURE, 0, diastolic_pressure)
            node.SetSolutionStepValue(PRESSURE_VENOUS, 0, diastolic_pressure)
        else:
            node.SetSolutionStepValue(FLOW, 0, q_inicial)
            node.SetSolutionStepValue(PRESSURE, 0, diastolic_pressure)
            node.SetSolutionStepValue(
                DYASTOLIC_PRESSURE, 0, diastolic_pressure)
            node.SetSolutionStepValue(PRESSURE_VENOUS, 0, diastolic_pressure)

if (Coupled_Simulation == True):
    class pressure_config:
        solver_type = "AMGCL"
        scaling = False
        preconditioner_type = None
        max_iteration = 500
        tolerance = 1e-6
        smoother_type = "ILU0"
        krylov_type = "CG"
    class velocity_config:
        solver_type = "AMGCL"
        scaling = False
        preconditioner_type = None
        max_iteration = 500
        tolerance = 1e-6
        smoother_type = "ILU0"
        krylov_type = "GMRES"    
    import linear_solver_factory
    
    solver3D = solver.IncompressibleFluidSolver(model_part3D, 3)
    solver3D.max_vel_its = 5
    solver3D.max_press_its = 20
    solver3D.predictor_corrector = True
    solver3D.vel_toll = 1e-3
    solver3D.press_toll = 1e-3
    solver3D.dynamic_tau = 1.0 #0.01
    solver3D.compute_reactions = False
    solver3D.activate_smagorinsky(0.2)
    solver3D.velocity_linear_solver = linear_solver_factory.ConstructSolver(velocity_config)
    solver3D.pressure_linear_solver = linear_solver_factory.ConstructSolver(pressure_config)
    solver3D.Initialize()
    import CouplingTools1D_3Dv5
    transfer_obj = CouplingTools1D_3Dv5.TransferTools(
        model_part1D, model_part3D)
    #transfer_obj.Setting3d(diastolic_pressure)
    #transfer_obj.Initialize()
    blood_viscosity = config.blood_viscosity
    blood_density = config.blood_density
    P_initial = config.diastolic_pressure
    transfer_obj.Setting3d(P_initial, blood_density, blood_viscosity)
    #transfer_obj.Setting3d(diastolic_pressure)
    transfer_obj.Initialize()
    #transfer_obj.Velocity_Initial_Contitions()  # Setting the velocity
    #transfer_obj.Initial_Contitions(diastolic_pressure)
    #transfer_obj.Initialize()
    transfer_obj.Velocity_Initial_Contitions()  # Setting the velocity
    #transfer_obj.Initial_Contitions(diastolic_pressure)
    print("READY TO SIMULATE.................................")
else:
    print("ONLY 1D SOLVER....................................")

print("Coupled_Simulation", Coupled_Simulation)
print(FitValues)
print(Fit_control)
if((FitValues == True) and (Fit_control == True)):
    Coupled_Simulation = False
    print("Coupled_Simulation", Coupled_Simulation)

print("----------------------------------------------------------------------------------------------------------")
time = 0.0
total_time = 0.0
step = 0.0
minlength = 1e+12
minlength = integrator.Element_minLength(model_part1D)
inicial = 1

fixed_flow_nodes = []
if (InletPressure == False):
    for node in model_part1D.Nodes:
        if(node.IsFixed(FLOW) == True and node.GetSolutionStepValue(FLAG_VARIABLE) == 0):
            fixed_flow_nodes.append(node)
            print(" NODE_fixed_FLOW (INLET) ", node.Id)
else:
    for node in model_part1D.Nodes:
        if(node.IsFixed(PRESSURE) == True and node.GetSolutionStepValue(FLAG_VARIABLE) == 0):
            fixed_flow_nodes.append(node)
            print(" NODE_fixed_PRESSURE (INLET)", node.Id)

total_time = total_time + step_size
# Flag_Converge=False

# Initial Delta Step
Dt = integrator.EstimateDeltaTime(model_part1D, 0.8, minlength)

# Aprox to save results
steps_eval = final_time / Dt
save_results = config.save_results
save_Steps = final_time / save_results
if(save_results < steps_eval):
    output_step = math.floor(steps_eval / save_Steps)
else:
    out_res = 1
    output_step = math.floor(steps_eval / out_res)

if(Sub_steping == True):
    control_sub_step = 1
    # sub_step = 200
    if (output_step < sub_step):
        output_step = sub_step
    # for each time step, ie: 1 save 1 Dt, 100 save 100 Dt
else:
    Sub_steping == True
    control_sub_step = 1
    # sub_step = 200
    if ((output_step < sub_step) and (Coupled_Simulation == True)):
        output_step = sub_step

for node in model_part3D.Nodes:
    if(node.IsFixed(PRESSURE)):
        node.SetValue(IS_STRUCTURE, False)

for cond in model_part3D.Conditions:
    if(cond.Properties.Id >= 1000):  # outlet
        cond.SetValue(IS_STRUCTURE, False)

kkkkkk = 0
ffit_test = []
ffit_test.append(kkkkkk)
file_test = str("file_test.txt")
ffit_test[kkkkkk] = open(file_test, 'w')

myTimer = Timer()
if(only1Dtest == False):
    while(total_time < final_time):

        model_part1D.CloneTimeStep(total_time)
        if(step < 3):
            model_part3D.CloneTimeStep(total_time)

        if (InletPressure == False):
            for node in fixed_flow_nodes:
            # if (Coupled_Simulation == False):
                if(node.IsFixed(FLOW) == True and node.GetSolutionStepValue(FLAG_VARIABLE) == 0):
                    if (step > 1):
                        Q = inletconditiontable.GetValue(time)
                        Q = 2 * Q - node.GetSolutionStepValue(FLOW, 1)
                        node.SetSolutionStepValue(FLOW, 0, Q)
                    else:
                        Q = inletconditiontable.GetValue(time)
                        node.SetSolutionStepValue(FLOW, 0, Q)
        else:
            for node in fixed_flow_nodes:
                # if (Coupled_Simulation == False):
                if(node.IsFixed(PRESSURE) == True and node.GetSolutionStepValue(FLAG_VARIABLE) == 0):
                    if (InletProfileType == "parabolic"):
                        print("Using Parabolic Pressure as input:")
                        instant_pressure = (pressure_parameter_2 * time * time) + (
                            pressure_parameter_1 * time) + pressure_parameter_3
                        print("Using Parabolic Pressure Profile as Inlet: ", instant_pressure)
                    elif (InletProfileType == "coseno"):
                        instant_pressure = (pressure_parameter_1) + (
                            (pressure_parameter_2) * (math.cos(2 * math.pi * total_time / time_period)))
                        print("Using Coseno Pressure Profile as Inlet: ", instant_pressure)
                    elif (InletProfileType == "table"):
                        instant_pressure = inletconditiontable.GetValue(time)
                        print("Using Table Pressure Profile(mdpa) as Inlet: ", instant_pressure)
                    # if (total_time <= 0.1):
                        # instant_pressure=diastolic_pressure
                    # elif (0.1 < total_time <0.3):
                        #instant_pressure = diastolic_pressure + total_time*systolic_pressure
                    # else:
                        # instant_pressure=diastolic_pressure
                    A0 = node.GetValue(NODAL_AREA)
                    A_aux = node.GetSolutionStepValue(NODAL_AREA, 1)
                    # Q=node.GetSolutionStepValue(FLOW,1)
                    A = (
                        math.sqrt(A0) + (((instant_pressure - diastolic_pressure) * A0) / node.GetSolutionStepValue(BETA, 1))) ** 2
                    # print A
                    A2 = pow(2 * pow(A, 0.25) - pow(A_aux, 0.25), 4)
                    # print A2
                    node.SetSolutionStepValue(NODAL_AREA, 0, A)
                    # node.SetSolutionStepValue(FLOW,0,Q)

        # and (nro_cardiac_cycle == cardiac_cycle_to_3D)):
        if ((Coupled_Simulation == True) and (control_sub_step == sub_step)):
            if(step >= 3):
                print("...............................Solve 3D..................................")
                integrator.SolveStep(model_part1D)
                # integrator.ComputePressure(model_part1D)
                integrator.ComputePressure(model_part1D, diastolic_pressure)
                model_part3D.CloneTimeStep(total_time)
                print("Total_time 3d:", total_time)
                print("----------------------------Transfer 1d(tn) to 3d(tn+1)------------------")
                transfer_obj.Transfer1D_to_3D(diastolic_pressure)
                print("Solve 3D para ", total_time)
                myTimer.Start("solver3D.Solve()")
                solver3D.Solve()
                myTimer.Stop("solver3D.Solve()")
                print(myTimer)
                print("----------------------------Transfer 3d to 1d----------------------------")
                # transfer_obj.Transfer3D_to_1D()
                # integrator.SolveStep(model_part1D)
                #integrator.ComputePressure(model_part1D, diastolic_pressure)
                # removal_tool.ComputePressure(model_part1D)
                control_sub_step = 0
                var_aux = True
                out = out + 1
        else:
            if ((var_aux == True) and (Sub_steping == True)):
                print("Solve 1D ------------------------------> ", total_time)
                print("Sub_step-->", sub_step)
                var_aux = False
            else:
                print("Solve 1D ------------------------------> ", total_time)

            integrator.SolveStep(model_part1D)
            # Only when I saved the results.
            integrator.ComputePressure(model_part1D, diastolic_pressure)
            # removal_tool.ComputePressure(model_part1D) --> This is only needed
            # when the results are saveed.
            control_sub_step = control_sub_step + 1
            out = out + 1
            # print "Step--->End"
            print(step)
        # and (nro_cardiac_cycle == cardiac_cycle_to_3D)):
        if((math.fmod(step, output_step) == 0) and (step >= 3)):
            # removal_tool.ComputePressure(model_part1D)

            if ((FitValues == True)):  # and (step>=sub_step)):
                transfer_obj.FitValues_1d(total_time)
                if ((Coupled_Simulation == True) and (step >= 3)):
                    transfer_obj.FitValues_3d(total_time, ffit_test, kkkkkk)
                    # raw_input()
            # transfer_obj.FitValues_Inlet()

            if (CardiacCycleConvergence == True):
                print("Check Convergence")
                # Function to check the cardiac_Cycle convergence (move
                # outside)
                time_cardiac_cycle = integrator.CheckCardiacCovergence(
                    model_part1D, time_cardiac_cycle)

            if (plot1d == True):
                if(ascii == True):
                    f1d.write("Time_step \n")
                    f1d.write(str(total_time))
                    f1d.write("\n")
                    # f.write("Number_of_nodes\n")
                    f1d.write(
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
                        f1d.write(ToWrite)
                    print("Writing CVPR results 1D Results(time)------------------------>", total_time)
                else:
                    #gid_io.WriteNodalResults(NODAL_AREA, model_part1D.Nodes, total_time, 0)
                    #gid_io.WriteNodalResults(NODAL_MASS, model_part1D.Nodes, total_time, 0)
                    gid_io.WriteNodalResults(
                        FLOW, model_part1D.Nodes, total_time, 0)
                    #gid_io.WriteNodalResults(RADIUS, model_part1D.Nodes, total_time, 0)
                    #gid_io.WriteNodalResults(RHS, model_part1D.Nodes, total_time, 0)
                    #gid_io.WriteNodalResults(FLAG_VARIABLE, model_part1D.Nodes, total_time, 0)
                    #gid_io.WriteNodalResults(VELOCITY, model_part1D.Nodes, total_time, 0)
                    gid_io.WriteNodalResults(
                        PRESSURE, model_part1D.Nodes, total_time, 0)
                    #gid_io.WriteNodalResults(WORK, model_part1D.Nodes, total_time, 0)
                    print("Writing GiD results 1D Results(time)------------------------>", total_time)

            if (Coupled_Simulation == False):
                control_sub_step = 0

            if(plot3d == True):
                # node_out=GetNodeAfter(full_nodes_table.table, config.inlets_1d[0][1])
                # print node_out
                # press_1d_out=model_part1D.Nodes[node_out].GetSolutionStepValue(PRESSURE)
                if(ascii == True):
                    f3d.write("Time_step \n")
                    f3d.write(str(total_time))
                    f3d.write("\n")
                    # f.write("Number_of_nodes\n")
                    # f.write(xxxx)
                    f3d.write(
                        "Nodal_data_table (Index-NodeId-Pressure-Velocity_X-Velocity_Y-Velocity_Z) \n")
                    indextowrite = 1
                    for node in model_part3D.Nodes:
                        nodewrite = node.Id
                        pressure_rect = node.GetSolutionStepValue(PRESSURE)
                        # +press_1d_out
                        ToWrite = str(indextowrite) + " " + str(
                            nodewrite) + " " + str(pressure_rect) + " "
                        ToWrite += str(
                            node.GetSolutionStepValue(VELOCITY_X)) + " "
                        ToWrite += str(
                            node.GetSolutionStepValue(VELOCITY_Y)) + " "
                        ToWrite += str(
                            node.GetSolutionStepValue(VELOCITY_Z)) + "\n"
                        indextowrite = indextowrite + 1
                        f3d.write(ToWrite)
                    print("Writing CVPR results 3D Results(time)------------------------>", total_time)
                else:
                    gid_io_3d.WriteNodalResults(
                        VELOCITY, model_part3D.Nodes, total_time, 0)
                    gid_io_3d.WriteNodalResults(
                        PRESSURE, model_part3D.Nodes, total_time, 0)
                    gid_io_3d.WriteNodalResults(
                        VISCOSITY, model_part3D.Nodes, total_time, 0)
                    gid_io_3d.WriteNodalResults(
                        DENSITY, model_part3D.Nodes, total_time, 0)
                    gid_io_3d.WriteNodalResults(
                        FLAG_VARIABLE, model_part3D.Nodes, total_time, 0)
                    out = 0
                    print("Writing GiD results 3D results(step)------------------------>", total_time)

        if(time > time_cardiac_cycle):
            time = time - time_cardiac_cycle
            nro_cardiac_cycle = nro_cardiac_cycle + 1
            control_sub_step = 0
            out = 0

        if (step_size_control == True):
            Dt = integrator.EstimateDeltaTime(model_part1D, 0.8, minlength)
            # print "Using Variable_Delta_Step:", Dt
        else:
            Dt = step_size

        time = time + Dt
        total_time = total_time + Dt
        out = out + 1
        step = step + 1
        # print "out_step", output_step
        # print "StepCoupled", StepCoupled
        # print "out", out

print("Finish 1D-3D")
print("Computing 3D Energy Losses")

if(only1Dtest == False):
    ab_3d_list = []
    if (FitValues == True):
        print("Fitting A-B values")
        ab_1d_list = transfer_obj.Fit_ABValues_1D()
        if (Coupled_Simulation == True):
            ab_3d_list = transfer_obj.Fit_ABValues_3D()
else:
    A = config.A
    B = config.B
    for j in range(0, len(FFR_Outlet_NODES_Values)):
        print("FFR_Outlet_NODES_Values", FFR_Outlet_NODES_Values[j].Id)
        ab_3d_list = []
        if (A == 0):
            # A = 257000000 # 54
            # A = 630000000 # 59
            # A = 184000000 # 34
            # A = 313000000 # 53
            # A = 198000000 # 32
            A = 0
        if (B == 0):
            # B= 4470000000000000 # 54
            # B= 561000000000000 # 59
            # B= 1940000000000000 #34
            # B= 2340000000000000 #53
            # B= 8510000000000000 # 32
            B = 0
        nodeAB = config.nodeAB
        if (nodeAB == 0):
            nodeAB = FFR_Outlet_NODES_Values[j].Id
        print(nodeAB)
        ab_3d_list.append([nodeAB, A, B])
        print(ab_3d_list)

if(Coupled_Simulation == False):
    for item in ab_1d_list:
        [node_id, A, B] = item
        for cond in model_part1D.Conditions:
            cond_of_interest = False
            for node in cond.GetNodes():
                if(node.Id == node_id):
                    cond_of_interest = True
                    print("Node and Conditions for the 1D FFR model")
                    print("Node", node_id)
                    print("Parameter A=", A)
                    print("Parameter B=", B)
            if(cond_of_interest == True):
                cond.SetValue(a, A)
                cond.SetValue(b, B)
                print("A and B setted_1")

for item in ab_3d_list:
    [node_id, A, B] = item
    for cond in model_part1D.Conditions:
        cond_of_interest = False
        for node in cond.GetNodes():
            if(node.Id == node_id):
                cond_of_interest = True
                print("Node and Conditions for the 1D FFR model")
                print("Node", node_id)
                print("Parameter A=", A)
                print("Parameter B=", B)
        if(cond_of_interest == True):
            cond.SetValue(a, A)
            cond.SetValue(b, B)
            print("A and B setted_2")

# for cond in model_part1D.Conditions:
    # for node in cond.GetNodes():
        # print node.GetValue(TERMINAL_RESISTANCE)
# print "TERMINAL_RESISTANCE"
# raw_input()
print("Setting Hypermenia Conditions")
dyastolic_hypermia_pressure = config.dyastolic_hypermia_pressure
print(dyastolic_hypermia_pressure)
if(only1Dtest == True):
    for cond in model_part1D.Conditions:
        for node in cond.GetNodes():
                # cond.SetValue(PRESSURE_VENOUS,dyastolic_hypermia_pressure)
            node.SetSolutionStepValue(
                PRESSURE_VENOUS, 0, dyastolic_hypermia_pressure)
            Resistence_factor = 1
            cond.SetValue(PRESSURE_DT, Resistence_factor)
            # print diastolic_pressure
            # print PRESSURE_VENOUS
            # print node.Id
            # print cond.Id


#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
# 3D reduced model is setted in transfer_obj.Fit_ABValues_3D
print("Finish Compute 3D Energy Losses")
print("Computing FFR Analysis")
Only1D.Only1D(model_part1D, config, input_file_name,
              FFR_Inlet_NODES_Values, FFR_Outlet_NODES_Values)
print("Finish FFR Analysis")
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------

# C
# m=node.GetSolutionStepValue(NODAL_MASS)
# node.SetSolutionStepValue(NODAL_MASS,0,m+C)

if (plot1d == True):
    if (ascii == False):
        gid_io.FinalizeResults()
    else:
        f1d.close()

if (plot3d == True):
    if (ascii == False):
        gid_io_3d.FinalizeResults()
    else:
        f3d.close()


print("Time=", total_time)
print("---------------------------------------------END-------------------------------------------")
