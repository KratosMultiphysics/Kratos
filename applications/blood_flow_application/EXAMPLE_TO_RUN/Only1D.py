#makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division
from KratosMultiphysics import *
from KratosMultiphysics.BloodFlowApplication import *
CheckForPreviousImport()
import math
import time


def Only1D(model_part1D, config, simulation_config, input_file_name, FFR_Inlet_NODES_Values, FFR_Outlet_NODES_Values):

    # model_part1D: 1D.MDPA
    # config: initial conditions
    # simulation_config: Simulation conditions
    # input_file_name: Name
    # FFR_Inlet_NODES_Values
    # FFR_Outlet_NODES_Values
    print("VERSION 11.12.2013")
    print("---------------------------------- FFR COMPUTATION ------------------------------------------")
    # Set initial_Radius
    # FitRadius=simulation_config.FitRadius
    blood_viscosity = simulation_config.blood_viscosity
    blood_density = simulation_config.blood_density

    # Set computational_time
    cardiac_cycle = config.nro_cardiac_cycles
    # total_time (last value of the Cardiac_cycle
    time_cardiac_cycle = config.time_period
    # True-->Sub_step_control (only for the Coupled_3d_1d)
    Sub_steping = simulation_config.Sub_steping
    sub_step = simulation_config.sub_step  #
    # True-->Activate	False-->fix (step_size)
    step_size_control = simulation_config.step_size_control
    step_size = simulation_config.step_size  # config.step_size
    CardiacCycleConvergence = simulation_config.CardiacCycleConvergence

    # Pressure conditions
    diastolic_pressure = config.diastolic_pressure  # Pa
    systolic_pressure = config.systolic_pressure  # Pa
    diastolic_hypermia_pressure = config.diastolic_hypermia_pressure  # Pa
    time_period = config.time_period
    InletProfileType = simulation_config.inlet_pressure_type

    # Set_Results (# Aprox to save results)
    cardiac_cycle_to_3D = cardiac_cycle  # 3D cardiac_Cycle will be running
    final_time = cardiac_cycle * time_cardiac_cycle
    ascii = simulation_config.ascii_results
    save_results = simulation_config.save_results

    # Set Aux_variables
    # This variable is only for doing test A-B with the 1D model.
    Fit_control = False
    var_aux = True
    pressure_factor = simulation_config.pressure_factor
    Resistence_factor = simulation_config.Resistence_factor
    Condition_Variable = simulation_config.Condition_Variable
    nro_cardiac_cycle = 1.0
    diastolic_pressure = pressure_factor * diastolic_pressure
    systolic_pressure = pressure_factor * systolic_pressure

    # Set Initial_conditions
    Q_initial = simulation_config.Q_initial
    P_initial = diastolic_pressure

    if ((InletProfileType == "coseno") or (InletProfileType == "parabolic") or (InletProfileType == "table")):
        InletPressure = True
    else:
        InletPressure = False

    # Select FFR Nodes
    node_IN = FFR_Inlet_NODES_Values[0].Id
    counter_FFR = 1
    FFR_Nodes = []
    FFR_VALUES = []
    Q_VALUES = []
    FFR_Nodes.append(node_IN)
    FFR_VALUES.append(0)
    Q_VALUES.append(0)
    for j in range(0, len(FFR_Outlet_NODES_Values)):
        FFR_Nodes.append(FFR_Outlet_NODES_Values[j].Id)
        FFR_VALUES.append(0)
        Q_VALUES.append(0)

    # Preparing GID Files
    gid_mode = GiDPostMode.GiD_PostBinary
    multifile = MultiFileFlag.MultipleFiles
    deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
    # WriteElementsOnly  # WriteConditions
    write_conditions = WriteConditionsFlag.WriteElementsOnly

    input_file_name = input_file_name + "_FFR"
    plot1d = True
    if (ascii == False):
        gid_io = GidIO(input_file_name, gid_mode,
                       multifile, deformed_mesh_flag, write_conditions)
    else:
        results = str(input_file_name + ".cvpr")
        f1d = open(results, 'w')
    print("1D_FFR Problem is running")
    output_step = 1  # for eachffit.append = 1

    model_part1D.SetBufferSize(2)
    out = 1
    # mesh to be printed
    if(ascii == False):
        if(plot1d):
            mesh_name = 1.0
            gid_io.InitializeMesh(mesh_name)
            gid_io.WriteMesh(model_part1D.GetMesh())
            gid_io.WriteNodeMesh(model_part1D.GetMesh())
            gid_io.FinalizeMesh()
            gid_io.InitializeResults(mesh_name, (model_part1D).GetMesh())
            print("Writing 1D Mesh------------------------")

    # Initial conditions for the 1D model HYPERMIA SITUATIONS
    integrator = ArteryTimeIntegrator()
    integrator.Initialize(model_part1D)
    inletconditiontable = model_part1D.GetTable(1)
    minlength = 1e+12
    minlength = integrator.Element_minLength(model_part1D)

    if(Condition_Variable):
        for cond in model_part1D.Conditions:
            for node in cond.GetNodes():
                # cond.SetValue(PRESSURE_VENOUS,diastolic_hypermia_pressure)
                node.SetSolutionStepValue(
                    PRESSURE_VENOUS, 0, diastolic_hypermia_pressure)
                cond.SetValue(PRESSURE_DT, Resistence_factor)

    # Initial Values. Set the intial pressure as reference for the 1D model
    # (all nodes take the systolic pressure and the initial flow)
    if (InletPressure):
        if (InletProfileType == "parabolic"):
            print("Using Parabolic Pressure as input:")
            pressure_parameter_1 = (
                (systolic_pressure - diastolic_pressure) * 4) / time_period
            pressure_parameter_2 = (-pressure_parameter_1) / time_period
            pressure_parameter_3 = diastolic_pressure
        elif (InletProfileType == "coseno"):
            pressure_parameter_1 = (systolic_pressure + diastolic_pressure) / 2
            pressure_parameter_2 = (diastolic_pressure - systolic_pressure) / 2
        print("Pressure_Factor: ", pressure_factor)
        print("Diastolic pressure: ", diastolic_pressure)
        print("Systolic Pressure: ", systolic_pressure)
        # diastolic_pressure=inletconditiontable.GetValue(0)
        for node in model_part1D.Nodes:
            if(node.IsFixed(PRESSURE) == False):
                node.SetSolutionStepValue(FLOW, 0, Q_initial)
                node.SetSolutionStepValue(PRESSURE, 0, P_initial)
                node.SetSolutionStepValue(
                    DYASTOLIC_PRESSURE, 0, diastolic_pressure)
                node.SetSolutionStepValue(
                    PRESSURE_VENOUS, 0, diastolic_hypermia_pressure)
            else:
                node.SetSolutionStepValue(FLOW, 0, Q_initial)
                node.SetSolutionStepValue(PRESSURE, 0, P_initial)
                node.SetSolutionStepValue(
                    DYASTOLIC_PRESSURE, 0, diastolic_pressure)
                node.SetSolutionStepValue(
                    PRESSURE_VENOUS, 0, diastolic_hypermia_pressure)
    else:
        Q_initial = inletconditiontable.GetValue(0)
        for node in model_part1D.Nodes:
            if(node.IsFixed(FLOW) == False):
                node.SetSolutionStepValue(FLOW, 0, Q_initial)
                node.SetSolutionStepValue(PRESSURE, 0, P_initial)
                node.SetSolutionStepValue(
                    DYASTOLIC_PRESSURE, 0, diastolic_pressure)
                node.SetSolutionStepValue(
                    PRESSURE_VENOUS, 0, diastolic_hypermia_pressure)
            else:
                node.SetSolutionStepValue(FLOW, 0, Q_initial)
                node.SetSolutionStepValue(PRESSURE, 0, P_initial)
                node.SetSolutionStepValue(
                    DYASTOLIC_PRESSURE, 0, diastolic_pressure)
                node.SetSolutionStepValue(
                    PRESSURE_VENOUS, 0, diastolic_hypermia_pressure)

    print("Hyperemia conditions set")
    print("----------------------------------------------------------------------------------------------------------")
    time = 0.0
    total_time = 0.0
    step = 0.0
    inicial = 1

    fixed_flow_nodes = []
    if (InletPressure == False):
        for node in model_part1D.Nodes:
            if(node.IsFixed(FLOW) and node.GetSolutionStepValue(FLAG_VARIABLE) == 0):
                fixed_flow_nodes.append(node)
                print(" NODE_fixed_FLOW (INLET) ", node)
    else:
        for node in model_part1D.Nodes:
            if(node.IsFixed(PRESSURE) and node.GetSolutionStepValue(FLAG_VARIABLE) == 0):
                fixed_flow_nodes.append(node)
                print(" NODE_fixed_PRESSURE (INLET)", node)

    total_time = total_time + step_size
    # Initial Delta Step
    print("Using Adaptative Time Step")
    Dt = integrator.EstimateDeltaTime(model_part1D, 0.8, minlength)
    steps_eval = final_time / Dt
    save_Steps = final_time / save_results
    print("----------------------------------------------------------------------------------------------------------")
    # Aprox to save results

    if(save_results < steps_eval):
        output_step = math.floor(steps_eval / save_Steps)
    else:
        out_res = 1
        output_step = math.floor(steps_eval / out_res)

    # Time controler
    myTimer = Timer()
    # START FFR COMPUTATION
    while(total_time < final_time):
        model_part1D.CloneTimeStep(total_time)
        if(step < 3):
            model_part3D.CloneTimeStep(total_time)

        if (InletPressure == False):
            for node in fixed_flow_nodes:
            # if (Coupled_Simulation == False):
                if(node.IsFixed(FLOW) and node.GetSolutionStepValue(FLAG_VARIABLE) == 0):
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
                if(node.IsFixed(PRESSURE) and node.GetSolutionStepValue(FLAG_VARIABLE) == 0):
                    if (InletProfileType == "parabolic"):
                        print("Using Parabolic Pressure as input:")
                        instant_pressure = (pressure_parameter_2 * time * time) + \
                            (pressure_parameter_1 * time) + \
                            pressure_parameter_3
                        print(
                            "Using Parabolic Pressure Profile as Inlet: ", instant_pressure)
                    elif (InletProfileType == "coseno"):
                        instant_pressure = (pressure_parameter_1) + ((
                            pressure_parameter_2) * (math.cos(2 * math.pi * total_time / time_period)))
                        print(
                            "Using Coseno Pressure Profile as Inlet: ", instant_pressure)
                    elif (InletProfileType == "table"):
                        instant_pressure = inletconditiontable.GetValue(time)
                        print(
                            "Using Table Pressure Profile(mdpa) as Inlet: ", instant_pressure)
                    # Transfer pressure in Area & Flow
                    A0 = node.GetValue(NODAL_AREA)
                    A_aux = node.GetSolutionStepValue(NODAL_AREA, 1)
                    A = (
                        math.sqrt(A0) + (((instant_pressure - diastolic_pressure) * A0) / node.GetSolutionStepValue(BETA, 1))) ** 2
                    A2 = pow(2 * pow(A, 0.25) - pow(A_aux, 0.25), 4)
                    node.SetSolutionStepValue(NODAL_AREA, 0, A)
                    # node.SetSolutionStepValue(FLOW,0,Q)

        myTimer.Start("solver_1D.Solve()")
        integrator.SolveStep(model_part1D)
        myTimer.Stop("solver_1D.Solve()")
        out = out + 1

        if((math.fmod(step, output_step) == 0)):
                # removal_tool.ComputePressure(model_part1D)
            # Only when I saved the results.
            integrator.ComputePressure(model_part1D, dyastolic_pressure)
            if (CardiacCycleConvergence):
                print("Check Convergence")
                # Function to check the cardiac_Cycle convergence (move
                # outside)
                time_cardiac_cycle = integrator.CheckCardiacCovergence(
                    model_part1D, time_cardiac_cycle)

            if (plot1d):
                if(ascii):
                    f1d.write("Time_step \n")
                    f1d.write(str(total_time))
                    f1d.write("\n")
                    # f.write("Number_of_nodes\n")
                    f1d.write(
                        "Nodal_data_table (Index-NodeId-PRESSURE-FLOW-RADIUS) \n")
                    indextowrite = 1
                    for node in model_part1D.Nodes:
                        nodewrite = node.Id
                        ToWrite = str(indextowrite) + " "
                        ToWrite += str(nodewrite) + " " + str(
                            node.GetSolutionStepValue(PRESSURE)) + " "
                        ToWrite += str(node.GetSolutionStepValue(FLOW)) + " " + str(
                            node.GetSolutionStepValue(RADIUS)) + "\n"
                        indextowrite = indextowrite + 1
                        f1d.write(ToWrite)
                    print(
                        "Writing CVPR results 1D_FFR Results(time)------------------------>", total_time)
                else:
                    gid_io.WriteNodalResults(
                        NODAL_AREA, model_part1D.Nodes, total_time, 0)
                    #gid_io.WriteNodalResults(NODAL_MASS, model_part1D.Nodes, total_time, 0)
                    gid_io.WriteNodalResults(
                        FLOW, model_part1D.Nodes, total_time, 0)
                    gid_io.WriteNodalResults(
                        RADIUS, model_part1D.Nodes, total_time, 0)
                    #gid_io.WriteNodalResults(RHS, model_part1D.Nodes, total_time, 0)
                    #gid_io.WriteNodalResults(FLAG_VARIABLE, model_part1D.Nodes, total_time, 0)
                    #gid_io.WriteNodalResults(VELOCITY, model_part1D.Nodes, total_time, 0)
                    gid_io.WriteNodalResults(
                        PRESSURE, model_part1D.Nodes, total_time, 0)
                    #gid_io.WriteNodalResults(WORK, model_part1D.Nodes, total_time, 0)
                    print(
                        "Writing GiD results 1D_FFR Results(time)------------------------>", total_time)

            counter_FFR = counter_FFR + 1
            # if (total_time > (final_time)):
            i_FFR = 0
            for i_FFR in range(0, len(FFR_Nodes)):
                # node=FFR_Nodes[i_FFR]
                node = model_part1D.Nodes[FFR_Nodes[i_FFR]]
                Pressure = node.GetSolutionStepValue(PRESSURE)
                # print Pressure
                FFR_VALUES[i_FFR] = FFR_VALUES[i_FFR] + Pressure
                Flow = node.GetSolutionStepValue(FLOW)
                Q_VALUES[i_FFR] = Q_VALUES[i_FFR] + Flow
                # print node
                # print "FFR_VALUES_TOTAL", FFR_VALUES[i_FFR]

            control_sub_step = 0

        if(time > time_cardiac_cycle):
            time = time - time_cardiac_cycle
            nro_cardiac_cycle = nro_cardiac_cycle + 1
            control_sub_step = 0
            out = 0

        if (step_size_control):
            Dt = integrator.EstimateDeltaTime(model_part1D, 0.8, minlength)
        else:
            Dt = step_size

        time = time + Dt
        total_time = total_time + Dt
        out = out + 1
        step = step + 1

    print("FINISH FFR COMPUTATION")
    print("----------------------------------------------------------------------------------------------------------")
    print(myTimer)

    if (plot1d):
        if (ascii == False):
            gid_io.FinalizeResults()
        else:
            f1d.close()

    k = 0
    FFR_RESULTS_FILE = []
    FFR_RESULTS_FILE.append(k)
    FFR_FILE = str("FFR_RESULTS_FILE.txt")
    FFR_RESULTS_FILE[k] = open(FFR_FILE, 'w')

    print("Time=", total_time)
    print("---------------------------------------------FFR VALUES-------------------------------------------")
    for i_FFR in range(0, len(FFR_Nodes)):
            # print "Mean Pressure", FFR_VALUES[i_FFR]
            # print counter_FFR
        print("Node::> ", FFR_Nodes[i_FFR])
        Q_VALUES[i_FFR] = Q_VALUES[i_FFR] / counter_FFR
        FFR_VALUES[i_FFR] = FFR_VALUES[i_FFR] / counter_FFR
        print("Mean Pressure", FFR_VALUES[i_FFR])
        # print "Mean Flow", Q_VALUES[i_FFR]

    mean_inlet_pressure = FFR_VALUES[0]
    mean_inlet_flow = Q_VALUES[0]
    result_Flow = 0
    result_FFR = 0
    for i_FFR in range(1, len(FFR_Nodes)):
        result_FFR = FFR_VALUES[i_FFR] / mean_inlet_pressure
        node = FFR_Nodes[i_FFR]
        print("FFR :::> ", node, "--", result_FFR)
        ToWrite = str(node) + " "
        ToWrite += str(result_FFR) + "\n"
        FFR_RESULTS_FILE[k].write(ToWrite)

        #result_Flow = Q_VALUES[i_FFR]/mean_inlet_flow
        # print "FLOW_Index :::> ", FFR_Nodes[i_FFR], "--",  result_Flow

    print("---------------------------------------------END-------------------------------------------")
