#makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division
from KratosMultiphysics import *
from KratosMultiphysics.BloodFlowApplication import *
CheckForPreviousImport()
import math


def Only1D(model_part1D, model_part3D, total_time, config, simulation_config, input_file_name, FFR_Inlet_NODES_Values, FFR_Outlet_NODES_Values,summary_file,mean_flow_1d_value,Aortic_inlet):

    import time
    import CouplingTools1D_3Dv5
    # model_part1D: 1D.MDPA
    # config: initial conditions
    # simulation_config: Simulation conditions
    # input_file_name: Name
    # FFR_Inlet_NODES_Values
    # FFR_Outlet_NODES_Values
    Config_version="ONLY1D.PY:VERSION 04_June_2014_local"
    transfer_obj = CouplingTools1D_3Dv5.TransferTools(model_part1D, model_part3D)
    transfer_obj.Initialize()
    print(Config_version)    
    print("---------------------------------- FFR COMPUTATION ------------------------------------------")
    # Set initial_Radius
    # Write Conditions Used
    #file_test = str("Summary_Case(Conditions Used for FFR COMPUTATION).txt")
    #summary_file = open(file_test, 'w')
    
    #ToWriteIn_Summary = "This case was running, day: " + "\n"
    #ToWriteIn_Summary += str(time.strftime("%H:%M:%S")) + "\n"
    #ToWriteIn_Summary += str(time.strftime("%d/%m/%y")) + "\n"
    ToWriteIn_Summary = str(Config_version) + "\n"
    summary_file.write(ToWriteIn_Summary)           
    
    FitRadius=simulation_config.FitRadius
    blood_viscosity = simulation_config.blood_viscosity
    blood_density = simulation_config.blood_density

    # Set computational_time
    cardiac_cycle = config.nro_cardiac_cycles
    # total_time (last value of the Cardiac_cycle
    time_cardiac_cycle = config.time_period
    # True-->Sub_step_control (only for the Coupled_3d_1d)
    #Sub_steping = simulation_config.Sub_steping
    sub_step = simulation_config.sub_step  #
    # True-->Activate	False-->fix (step_size)
    step_size_control = simulation_config.step_size_control
    step_size = simulation_config.step_size  # config.step_size
    CardiacCycleConvergence = simulation_config.CardiacCycleConvergence

    # Pressure conditions
    #diastolic_pressure = config.diastolic_pressure  # Pa
    systolic_hypermia_pressure = config.systolic_hypermia_pressure  # Pa
    diastolic_hypermia_pressure = config.diastolic_hypermia_pressure  # Pa
    time_period = config.time_period
    InletProfileType = simulation_config.inlet_pressure_type

    # Set_Results (# Aprox to save results)
    cardiac_cycle_to_3D = cardiac_cycle  # 3D cardiac_Cycle will be running
    nro_FFR_cardiac_cycles = config.nro_FFR_cardiac_cycles
    final_time = total_time + (nro_FFR_cardiac_cycles * time_cardiac_cycle)
    ascii = simulation_config.ascii_results
    save_results = simulation_config.save_results
    nro_FFR_cardiac_cycles_counter=nro_FFR_cardiac_cycles
    nro_FFR_cardiac_cycles_counter=1
    # Set Aux_variables
    # This variable is only for doing test A-B with the 1D model.
    Fit_control = False
    var_aux = True
    pressure_factor = simulation_config.pressure_factor
    Resistence_factor = simulation_config.Resistence_factor
    Condition_Variable = simulation_config.Condition_Variable
    nro_cardiac_cycle = 1.0
    diastolic_hypermia_pressure = pressure_factor * diastolic_hypermia_pressure
    systolic_hypermia_pressure = pressure_factor * systolic_hypermia_pressure

    # Set Initial_conditions
    Q_initial = simulation_config.Q_initial
    P_initial = diastolic_hypermia_pressure

    transfer_obj.Setting3d(P_initial, blood_density, blood_viscosity)
    if ((InletProfileType == "coseno") or (InletProfileType == "parabolic") or (InletProfileType == "table")):
        InletPressure = True
    else:
        InletPressure = False
    
    #ToWriteIn_Summary += "Diastolic_pressure  : " + str(diastolic_pressure) + "\n"
    #ToWriteIn_Summary += "Systolic_pressure  : " + str(systolic_pressure) + "\n"
    #ToWriteIn_Summary += "Initial Pressure  : " + str(P_initial) + "\n"
    #ToWriteIn_Summary += "Nro_cardiac_cycle  : " + str(nro_cardiac_cycle) + "\n"
    #ToWriteIn_Summary += "Total_time: ------------------------> " + str(total_time) + "\n"    
    #ToWriteIn_Summary += "Inlet Pressure type  : " + str(InletProfileType) + "\n"
    #ToWriteIn_Summary += "Blood_viscosity  : " + str(blood_viscosity) + "\n"
    #ToWriteIn_Summary += "Blood_density  : " + str(blood_density) + "\n"
    #ToWriteIn_Summary += "Name_model_1D: "
    #ToWriteIn_Summary += str(input_file_name)
    #ToWriteIn_Summary += "\n"
    #ToWriteIn_Summary += "Fit Radius: " + str(FitRadius) + "\n"
    #ToWriteIn_Summary += "Initial Flow: " + str(Q_initial) + "\n"
    #ToWriteIn_Summary += "Initial Pressure: " + str(P_initial) + "\n"

        
    # Select FFR Nodes
    node_IN = FFR_Inlet_NODES_Values[0].Id
    print (node_IN)
    counter_FFR = 1
    FFR_Nodes = []
    IFFR_Nodes = []
    FFR_VALUES = []
    Q_VALUES = []
    FFR_Nodes.append(node_IN)
    IFFR_Nodes.append(node_IN)
    FFR_VALUES.append(0)
    Q_VALUES.append(0)
    for j in range(0, len(FFR_Outlet_NODES_Values)):
        FFR_Nodes.append(FFR_Outlet_NODES_Values[j].Id)
        IFFR_Nodes.append(FFR_Outlet_NODES_Values[j].Id)
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
        gid_io = GidIO(input_file_name, gid_mode,multifile, deformed_mesh_flag, write_conditions)
    else:
        results = str(input_file_name + ".cvpr")
        f1d = open(results, 'w')

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
    
    ToWriteIn_Summary = "-----------------------------------------------" + "\n"
    ToWriteIn_Summary += "Setting Hypermenia Conditions"+ "\n"
    ToWriteIn_Summary +="Condition Variable to modify the TERMINAL_RESISTANCE: " + str(Condition_Variable) + "\n"
    if(Condition_Variable):	  
	  #ToWriteIn_Summary+="Pressure imposed: " + str(PRESSURE_DT) + "\n"
	  ToWriteIn_Summary+="Modify terminal resistence using factor of: " + str(Resistence_factor) + "\n"
	  for cond in model_part1D.Conditions:
		for node in cond.GetNodes():
			# cond.SetValue(PRESSURE_VENOUS,diastolic_hypermia_pressure)
			#node.SetSolutionStepValue(PRESSURE_VENOUS, 0, diastolic_hypermia_pressure)
			cond.SetValue(PRESSURE_DT, Resistence_factor)
			#print (str(cond.GetValue(PRESSURE_DT)))
    # Initial Values. Set the intial pressure as reference for the 1D model
    # (all nodes take the systolic pressure and the initial flow)
    #raw_input()
    ToWriteIn_Summary += "Sistolic_hypermia_pressure: "+ str(systolic_hypermia_pressure) + "\n"
    ToWriteIn_Summary += "Diastolic_hypermia_pressure: "+ str(diastolic_hypermia_pressure) + "\n"
    #summary_file.write(ToWriteIn_Summary)
	
    # Initial Values. Set the intial pressure as reference for the 1D model
    # (all nodes take the systolic pressure and the initial flow)
    venous_pressure=diastolic_hypermia_pressure*0.9
    if (InletPressure):
	    if (InletProfileType == "parabolic"):
		    print("Using Parabolic Pressure as input")
		    ToWriteIn_Summary += "Using Parabolic Pressure as input" + "\n"
		    pressure_parameter_1 = ((systolic_hypermia_pressure - diastolic_hypermia_pressure) * 4) / time_period
		    pressure_parameter_2 = (-pressure_parameter_1) / time_period
		    pressure_parameter_3 = diastolic_hypermia_pressure
	    elif (InletProfileType == "coseno"):
		    print("Using Coseno Pressure profile as input:")
		    ToWriteIn_Summary += "Using Coseno Pressure profile as input" + "\n"
		    pressure_parameter_1 = (systolic_hypermia_pressure + diastolic_hypermia_pressure) / 2
		    pressure_parameter_2 = (diastolic_hypermia_pressure - systolic_hypermia_pressure) / 2	    
	    # diastolic_pressure=inletconditiontable.GetValue(0)
	    for node in model_part1D.Nodes:
		    if(node.IsFixed(PRESSURE) == False):
			    node.SetSolutionStepValue(FLOW, 0, Q_initial)
			    node.SetSolutionStepValue(PRESSURE, 0, P_initial)
			    node.SetSolutionStepValue(DYASTOLIC_PRESSURE, 0, diastolic_hypermia_pressure)
			    node.SetSolutionStepValue(PRESSURE_VENOUS, 0, venous_pressure)
		    else:
			    node.SetSolutionStepValue(FLOW, 0, Q_initial)
			    node.SetSolutionStepValue(PRESSURE, 0, P_initial)
			    node.SetSolutionStepValue(DYASTOLIC_PRESSURE, 0, diastolic_hypermia_pressure)
			    node.SetSolutionStepValue(PRESSURE_VENOUS, 0, venous_pressure)
    else:
	    Q_initial = inletconditiontable.GetValue(0)
	    for node in model_part1D.Nodes:
		    if(node.IsFixed(FLOW) == False):
			    node.SetSolutionStepValue(FLOW, 0, Q_initial)
			    node.SetSolutionStepValue(PRESSURE, 0, P_initial)
			    node.SetSolutionStepValue(DYASTOLIC_PRESSURE, 0, diastolic_hypermia_pressure)
			    node.SetSolutionStepValue(PRESSURE_VENOUS, 0, venous_pressure)
		    else:
			    node.SetSolutionStepValue(FLOW, 0, Q_initial)
			    node.SetSolutionStepValue(PRESSURE, 0, P_initial)
			    node.SetSolutionStepValue(DYASTOLIC_PRESSURE, 0,diastolic_hypermia_pressure)
			    node.SetSolutionStepValue(PRESSURE_VENOUS, 0, venous_pressure)

    print("Hyperemia conditions set")
    print("----------------------------------------------------------------------------------------------------------")
    ToWriteIn_Summary += "Hyperemia conditions set"+ "\n"
    summary_file.write(ToWriteIn_Summary)
    time = 0.0
    #total_time = 0.0
    step = 0.0
    inicial = 1
    counter = 1.0

    fixed_flow_nodes = []
    fixed_pressure_nodes = []
    if (InletPressure == False):
	    for node in model_part1D.Nodes:
		    if(node.IsFixed(FLOW) and node.GetSolutionStepValue(FLAG_VARIABLE) == 0):
			    fixed_flow_nodes.append(node)
			    print(" NODE_fixed_FLOW (INLET) ", node.Id)
    else:
	    for node in model_part1D.Nodes:
		    if(node.IsFixed(PRESSURE) and node.GetSolutionStepValue(FLAG_VARIABLE) == 0):
			    fixed_pressure_nodes.append(node)
			    print(" NODE_fixed_PRESSURE (INLET)", node.Id)

    total_time = total_time + step_size
    # Initial Delta Step
    # Initial Delta Step
    if (step_size_control == True):
	    print("Using Adaptative Time Step")
	    Dt = integrator.EstimateDeltaTime(model_part1D, 0.6, minlength)
	    #ToWriteIn_Summary += "Using Adaptative Time Step: " + str(Dt) + " (aprox) Seconds" + "\n"
    else:
	    Dt = step_size
	    print("Using Fix Time Step")
	    print(str(step_size))
	    #ToWriteIn_Summary += "Using Fix Time Step: " + str(Dt) + " seconds" + "\n"

    #summary_file.write(ToWriteIn_Summary)

    steps_eval = final_time / Dt
    print("----------------------------------------------------------------------------------------------------------")
    if(save_results != 0.0):
	    save_Steps = final_time / save_results
	    if(save_results < steps_eval):
		    output_step = math.floor(steps_eval / save_Steps)
	    else:
		    out_res = 1
		    output_step = math.floor(steps_eval / out_res)
    else:
	    output_step = 1
    
    print("1D_FFR Problem is running")
    ToWriteIn_Summary = "1D_FFR Problem is running"+ "\n"
    ToWriteIn_Summary += "... ... ... ... ... ... ... ... ... ... " + "\n"
    #ToWriteIn_Summary += "Type of the profile input condition used: " + str(InletProfileType) + "\n"
    #ToWriteIn_Summary +="Type of the profile input condition used: "  + str(InletProfileType) + "  Units: Pa " + "\n"    	    
    #ToWriteIn_Summary += " InletPressure:	"  + str(InletPressure) + "\n"    
    summary_file.write(ToWriteIn_Summary)
    # Time controler
    myTimer = Timer()
    # START FFR COMPUTATION
    ToWriteIn_Summary = "FFR_Cardiac_cycle : ------------------------> " + str(nro_FFR_cardiac_cycles_counter) + "\n"
    summary_file.write(ToWriteIn_Summary)
    
    instant_pressure=0.0
    total_aortic_pressure=0.0
    Pressure_IN_TOTAL=0.0
    FLOW_IN_FFR_TOTAL=0.0
    Pressure_OUT_TOTAL=0.0
    
    mean_flow_1d_value=0.0
    mean_pressure_1d_value=0.0
    out_pressure_1d_value=0.0
    out_pressure_1d_value_FFR=0.0
    mean_flow_1d_inlet=0.0
    mean_flow_1d_FFR= 0.0
    mean_pressure_1d_FFR=0.0
    counter_FFR=0.0
    
    
    final_time=final_time+0.000001    
    while(total_time <= final_time):
        nro_FFR_cardiac_cycles_counter
        model_part1D.CloneTimeStep(total_time)  
        counter_FFR=counter_FFR+1
	#print(str(step))		
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
		for node in fixed_pressure_nodes:
			# if (Coupled_Simulation == False):
			if(node.IsFixed(PRESSURE) and node.GetSolutionStepValue(FLAG_VARIABLE) == 0):
				if (InletProfileType == "parabolic"):
					#print("Using Parabolic Pressure as input:")
					instant_pressure = (pressure_parameter_2 * time * time) + \
						(pressure_parameter_1 * time) + \
						pressure_parameter_3
					total_aortic_pressure=total_aortic_pressure+instant_pressure
					#print("Using Parabolic Pressure Profile as Inlet: ", instant_pressure)
				elif (InletProfileType == "coseno"):
					instant_pressure = (pressure_parameter_1) + ((pressure_parameter_2) * (math.cos(2 * math.pi * total_time / time_period)))
					total_aortic_pressure=total_aortic_pressure+instant_pressure
					#print("Using Coseno Pressure Profile as Inlet: ", instant_pressure)
				elif (InletProfileType == "table"):
					instant_pressure = inletconditiontable.GetValue(time)
					total_aortic_pressure=total_aortic_pressure+instant_pressure
					#print("Using Table Pressure Profile(mdpa) as Inlet: ", instant_pressure)
				# Transfer pressure in Area & Flow
				A0 = node.GetValue(NODAL_AREA)
				A_aux = node.GetSolutionStepValue(NODAL_AREA, 1)
				A = (math.sqrt(A0) + (((instant_pressure - diastolic_hypermia_pressure) * A0) / node.GetSolutionStepValue(BETA, 1))) ** 2
				A2 = pow(2 * pow(A, 0.25) - pow(A_aux, 0.25), 4)
				node.SetSolutionStepValue(NODAL_AREA, 0, A)
				# node.SetSolutionStepValue(FLOW,0,Q)

	#ToWriteIn_Summary = " Instant pressure imposed in the inlet node:	"  + str(instant_pressure) + " Pa " + "\n"
	#ToWriteIn_Summary += " Instant area imposed in the inlet node:	"  + str(A) + " mm " + "\n"	
	#summary_file.write(ToWriteIn_Summary)
	#ToWriteIn_summary_file_pressure_Flow_Aortic=str(Q)
	#ToWriteIn_summary_file_pressure_inlet1D=str(instant_pressure) + "\n"
	#summary_file_pressure_inlet1D.write(ToWriteIn_summary_file_pressure_inlet1D)
	#summary_file_pressure_Flow_Aortic(ToWriteIn_summary_file_pressure_Flow_Aortic)
		    
	myTimer.Start("solver_1D.Solve()")
	integrator.SolveStep(model_part1D)
	myTimer.Stop("solver_1D.Solve()")

	# if (total_time > (final_time)):
	# NO VALIDO PARA BIFURCACIONES
	###TESTTTTTSTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
	[mean_flow_1d_value,mean_pressure_1d_value,out_pressure_1d_value]=transfer_obj.Save_Values_For_different_situations(mean_flow_1d_value,mean_pressure_1d_value,out_pressure_1d_value)				
	
	i_FFR = 0
	Pressure_IN=0.0
	FLOW_IN_FFR=0.0
	Pressure_OUT=0.0
	for i_FFR in range(0, len(FFR_Inlet_NODES_Values)):
	    # node=FFR_Nodes[i_FFR]
	    node = FFR_Inlet_NODES_Values[i_FFR]
	    Pressure_IN = node.GetSolutionStepValue(PRESSURE)
	    Pressure_IN_TOTAL=Pressure_IN_TOTAL+Pressure_IN
	    FLOW_IN_FFR=node.GetSolutionStepValue(FLOW)
	    FLOW_IN_FFR_TOTAL=FLOW_IN_FFR_TOTAL+FLOW_IN_FFR
	    #print(node)
	    #print(Pressure_IN_TOTAL)
	    #print(FLOW_IN_FFR_TOTAL)
	    	    
	for i_FFR in range(0, len(FFR_Outlet_NODES_Values)):
	    node =FFR_Outlet_NODES_Values[i_FFR]
	    Pressure_OUT = node.GetSolutionStepValue(PRESSURE)
	    Pressure_OUT_TOTAL=Pressure_OUT_TOTAL+Pressure_OUT
	    #print(node)
	    #print(Pressure_OUT_TOTAL)

	
	if((math.fmod(step, output_step) == 0)):
		# removal_tool.ComputePressure(model_part1D)
	    # Only when I saved the results.
	    integrator.ComputePressure(model_part1D, diastolic_hypermia_pressure)
	    if (CardiacCycleConvergence):
		print("Check Convergence")
		# Function to check the cardiac_Cycle convergence (move
		# outside)
		time_cardiac_cycle = integrator.CheckCardiacCovergence(model_part1D, time_cardiac_cycle)

	    if (plot1d):
		if(ascii):
		    f1d.write("Time_step \n")
		    f1d.write(str(total_time))
		    f1d.write("\n")
		    # f.write("Number_of_nodes\n")
		    f1d.write("Nodal_data_table (Index-NodeId-PRESSURE-FLOW-RADIUS) \n")
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
		    print("Writing CVPR results 1D_FFR Results(time)------------------------>", total_time)
		else:
		    gid_io.WriteNodalResults(NODAL_AREA, model_part1D.Nodes, total_time, 0)
		    #gid_io.WriteNodalResults(NODAL_MASS, model_part1D.Nodes, total_time, 0)
		    gid_io.WriteNodalResults(FLOW, model_part1D.Nodes, total_time, 0)
		    gid_io.WriteNodalResults(RADIUS, model_part1D.Nodes, total_time, 0)
		    #gid_io.WriteNodalResults(RHS, model_part1D.Nodes, total_time, 0)
		    #gid_io.WriteNodalResults(FLAG_VARIABLE, model_part1D.Nodes, total_time, 0)
		    #gid_io.WriteNodalResults(VELOCITY, model_part1D.Nodes, total_time, 0)
		    gid_io.WriteNodalResults(PRESSURE, model_part1D.Nodes, total_time, 0)
		    #gid_io.WriteNodalResults(WORK, model_part1D.Nodes, total_time, 0)
		    print("Writing GiD results 1D_FFR Results(time)------------------------>", total_time)

	if(time >= time_cardiac_cycle):
	    time = time - time_cardiac_cycle
	    nro_FFR_cardiac_cycles_counter = nro_FFR_cardiac_cycles_counter + 1
	    print("FFR_Cardiac_cycle: ------------------------>", nro_FFR_cardiac_cycles_counter)
	    ToWriteIn_Summary = "FFR_Cardiac_cycle: ------------------------> " + str(nro_FFR_cardiac_cycles_counter) + "\n"
	    summary_file.write(ToWriteIn_Summary)
	    total_aortic_pressure=0.0
	    mean_pressure_1d_value=0.0
	    out_pressure_1d_value=0.0
	    mean_flow_1d_value=0.0
	    counter_FFR=0	   

	if (step_size_control):
	    Dt = integrator.EstimateDeltaTime(model_part1D, 0.6, minlength)
	else:
	    Dt = step_size

	time = time + Dt
	total_time = total_time + Dt
	step = step + 1

    print("FINISH FFR COMPUTATION")
    print("----------------------------------------------------------------------------------------------------------")
    ToWriteIn_Summary = "----------------------------------------------------------------------"+ "\n"
    ToWriteIn_Summary += "FINISH FFR COMPUTATION"+ "\n"
    summary_file.write(ToWriteIn_Summary)
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
    ToWriteIn_Summary = "---------------------------------------------FFR VALUES-------------------------------------------"+ "\n"
    summary_file.write(ToWriteIn_Summary)
    
    #counter_FFR=counter_FFR/nro_FFR_cardiac_cycles
    print(total_aortic_pressure)
    print(mean_flow_1d_value)
    print(mean_pressure_1d_value)
    print(out_pressure_1d_value)
		
	
    raw_input()
    #transfer_obj.ComputeFFR_Hypermia_Values(total_aortic_pressure,mean_pressure_1d_value,out_pressure_1d_value,mean_flow_1d_value,counter_FFR,summary_file)
    transfer_obj.ComputeFFR_Health_Values(total_aortic_pressure,mean_pressure_1d_value,out_pressure_1d_value,mean_flow_1d_value,counter_FFR,summary_file)
    raw_input()				
    #NO PARA BIFURCACION
    Mean_Pressure_OUT_TOTAL=0.0
    Mean_FLOW_IN_FFR_TOTAL=0.0 
    Mean_Pressure_IN_TOTAL=0.0
    Mean_total_aortic_pressure=0.0
    FFR_Value = 0.0
    #Mean_Pressure_IN_TOTAL = Pressure_IN_TOTAL / counter_FFR
    Mean_Pressure_IN_TOTAL = mean_pressure_1d_value / counter_FFR
    Mean_Pressure_OUT_TOTAL = out_pressure_1d_value / counter_FFR
    Mean_FLOW_IN_FFR_TOTAL=mean_flow_1d_value/counter_FFR
    Mean_total_aortic_pressure=total_aortic_pressure/counter_FFR

    FFR_Value_Aortic=Mean_Pressure_OUT_TOTAL/Mean_total_aortic_pressure
    result_FFR_3D = Mean_Pressure_OUT_TOTAL/Mean_Pressure_IN_TOTAL
    print("FFR_Value_Aortic: ", FFR_Value_Aortic)
    print("result_FFR_3D: ", result_FFR_3D)
    
    ToWrite = str(Aortic_inlet) + " "
    ToWrite += str(FFR_Value_Aortic) + "\n"
    FFR_RESULTS_FILE[k].write(ToWrite)
    ToWriteIn_Summary = "FFR_AORTIC_VALUE:  "+ str(FFR_Value_Aortic)+ "\n"
    summary_file.write(ToWriteIn_Summary)
    
    ToWrite = " 3D "
    ToWrite += str(result_FFR_3D) + "\n"
    FFR_RESULTS_FILE[k].write(ToWrite)
    ToWriteIn_Summary = "FFR_3D:  "+ str(result_FFR_3D)+ "\n"
    summary_file.write(ToWriteIn_Summary)       
    
    #print("---------------------------------------------FFR VALUES_2-------------------------------------------")
    #ToWriteIn_Summary = "---------------------------------------------FFR VALUES_2-------------------------------------------"+ "\n"
    raw_input()
    for i_FFR in range(1, len(FFR_Inlet_NODES_Values)):
	node=FFR_Outlet_NODES_Values[i_FFR].Id
	result_FFR = Mean_Pressure_OUT_TOTAL/Mean_Pressure_IN_TOTAL
	result_Flow = Mean_FLOW_IN_FFR_TOTAL/mean_flow_1d_value
	ToWrite = str(node) + " "
	ToWrite += str(result_FFR) + "\n"
        FFR_RESULTS_FILE[k].write(ToWrite)
	ToWriteIn_Summary = "FFR_MEAN_VALUE:  "+ str(result_FFR)+ "\n"
	summary_file.write(ToWriteIn_Summary)
	print(result_FFR)
	print(result_Flow)
	print(node)

    
    
    #print(counter_FFR)
    #for i_FFR in range(0, len(FFR_Nodes)):
            ## print "Mean Pressure", FFR_VALUES[i_FFR]
            ## print counter_FFR
        #print("Node::> ", FFR_Nodes[i_FFR])
        #Q_VALUES[i_FFR] = Q_VALUES[i_FFR] / counter_FFR
        #FFR_VALUES[i_FFR] = FFR_VALUES[i_FFR] / counter_FFR
        ##print("Mean Pressure", FFR_VALUES[i_FFR])
        ## print "Mean Flow", Q_VALUES[i_FFR]	
	##ToWriteIn_Summary = "Mean Pressure(" + str(FFR_Nodes[i_FFR]) + "):  " + str(FFR_VALUES[i_FFR])+ " Pa" + "\n"
	##summary_file.write(ToWriteIn_Summary)
    
    #mean_inlet_pressure = FFR_VALUES[0]
    ##mean_inlet_flow = Q_VALUES[0]
    #result_Flow = 0
    #result_FFR = 0
    #for i_FFR in range(1, len(FFR_Nodes)):
        #result_FFR = FFR_VALUES[i_FFR] / mean_inlet_pressure
        #node = FFR_Nodes[i_FFR]
        #print("FFR_MEAN_VALUE :::> ", node, "--", result_FFR)
        #ToWrite = str(node) + " "
        #ToWrite += str(result_FFR) + "\n"
        #FFR_RESULTS_FILE[k].write(ToWrite)
	#ToWriteIn_Summary = "FFR_MEAN_VALUE:  "+ str(result_FFR)+ "\n"
	#summary_file.write(ToWriteIn_Summary)
        #result_Flow = Q_VALUES[i_FFR]/mean_flow_1d_value
        ## print "FLOW_Index :::> ", FFR_Nodes[i_FFR], "--",  result_Flow

    #####result_FFR_max = 0   
    #####time = time_cardiac_cycle/2
    #####time = (nro_FFR_cardiac_cycles * time_cardiac_cycle) - time
    #####print (time)
    #####if (InletProfileType == "parabolic"):
	#####instant_pressure = (pressure_parameter_2 * time * time) + (pressure_parameter_1 * time) + pressure_parameter_3	    
    #####elif (InletProfileType == "coseno"):
	#####instant_pressure = (pressure_parameter_1) + ((pressure_parameter_2) * (math.cos(2 * math.pi * time / time_period)))
    #####for i_FFR in range(1, len(FFR_Nodes)):	
	#####node = FFR_Nodes[i_FFR]
	#####result_FFR =  FFR_VALUES[i_FFR] / instant_pressure 
	#####print("FFR_MAX_VALUE :::> ", node, "--", result_FFR)
	#####print("FFR_VALUES[i_FFR] :::> ", FFR_VALUES[i_FFR])
	#####print("instant_pressure :::> ", instant_pressure)
	#####ToWrite = str(node) + " "
	#####ToWrite += str(result_FFR) + "\n"
	#####FFR_RESULTS_FILE[k].write(ToWrite)
	#####ToWriteIn_Summary = "FFR_MAX_VALUE:  "+ str(result_FFR)+ "\n"
	#####summary_file.write(ToWriteIn_Summary)
	#result_Flow = Q_VALUES[i_FFR]/mean_inlet_flow
# print "FLOW_Index :::> ", FFR_Nodes[i_FFR], "--",  result_Flow

        
    return [Mean_FLOW_IN_FFR_TOTAL]
    print("---------------------------------------------END-------------------------------------------")
