# -*- coding: utf-8 -*-
import shutil, sys, os, time, json
sys.path.append("/work/piquee/Softwares/Kratos/applications/CoSimulationApplication/co_sim_io/python")
import CoSimIO
import PyPara, PyPrep, PySolv, PyDeform

with open('tau_settings.json') as json_file:
    tau_settings = json.load(json_file)

tau_path = tau_settings["tau_path"]
sys.path.append(tau_settings["kratos_path"])
sys.path.append(tau_path + "py_turb1eq/")
working_path = os.getcwd() + '/'

# tau_functions can only be imported after appending kratos' path
import tau_functions as TauFunctions

# Definition of the parameter file
para_path='airfoil_Structured.cntl'
para_path_mod = para_path + ".mod"
shutil.copy(para_path, para_path_mod)

# Initialize Tau python classes and auxiliary variable step
Para = PyPara.Parafile(para_path_mod)
Prep = PyPrep.Preprocessing(para_path_mod)
Solver = PySolv.Solver(para_path_mod)
Deform = PyDeform.Deformation(para_path_mod)
step = 0

def AdvanceInTime(current_time):
    # Preprocessing needs to be done before getting the time and time step
    TauFunctions.PrintBlockHeader("Start Preprocessing at time %s" %(str(time)))
    Prep.run(write_dualgrid=False,free_primgrid=False)
    TauFunctions.PrintBlockHeader("Stop Preprocessing at time %s" %(str(time)))
    TauFunctions.PrintBlockHeader("Initialize Solver at time %s" %(str(time)))
    Solver.init(verbose = 1, reset_steps = True, n_time_steps = 1)

    # Get current time and time step from tau
    tau_current_time = float(tau_solver_unsteady_get_physical_time())
    tau_time_step = float(Para.get_para_value('Unsteady physical time step size'))

    if tau_settings["echo_level"] > 0:
        print "TAU SOLVER AdvanceInTime"
        print 'tau_current_time = ', tau_current_time
        print 'tau_time_step = ', tau_time_step
    return tau_current_time + tau_time_step

def InitializeSolutionStep():
    if tau_settings["echo_level"] > 0:
        print("TAU SOLVER InitializeSolutionStep")

def SolveSolutionStep():
    if tau_settings["echo_level"] > 0:
        print("TAU SOLVER SolveSolutionStep")
    Solver.outer_loop()
    Solver.output()
    tau_plt_init_tecplot_params(para_path_mod)
    tau_solver_write_output_conditional()
    TauFunctions.ConvertOutputToDat(working_path, tau_path, step, para_path_mod)

def FinalizeSolutionStep():
    if tau_settings["echo_level"] > 0:
        print("TAU SOLVER FinalizeSolutionStep")

    tau_parallel_sync()
    Solver.finalize()
    tau_free_dualgrid()
    tau_free_prims()
    Para.free_parameters()
    global step
    step += 1

def ImportData(conn_name, identifier):
    if tau_settings["echo_level"] > 0:
        print "TAU SOLVER ImportData"
    data = CoSimIO.ImportData(conn_name, identifier)

    # identifier is the data-name in json
    if identifier == "Interface_disp":
        # Deform mesh
        if tau_mpi_rank() == 0:
            TauFunctions.ExecuteBeforeMeshDeformation(data, working_path, step, para_path_mod)
        Deform.run(read_primgrid=1, write_primgrid=1, read_deformation=0, field_io=1)
    else:
        raise Exception('TauSolver::ExportData::identifier "{}" not valid! Please use Interface_disp'.format(identifier))
    if tau_settings["echo_level"] > 0:
	print 'displacementKratos = ', data
        print "TAU SOLVER After ImportData"

def ExportData(conn_name, identifier):
    if tau_settings["echo_level"] > 0:
        print "TAU SOLVER ExportData"
    # identifier is the data-name in json
    if identifier == "Interface_force":
        data = TauFunctions.ComputeFluidForces(working_path, step)
	data *= 1000000
    else:
        raise Exception('TauSolver::ExportData::identifier "{}" not valid! Please use Interface_force'.format(identifier))

    CoSimIO.ExportData(conn_name, identifier, data)
    if tau_settings["echo_level"] > 0:
        print 'data = ', data
        print "TAU SOLVER After ExportData"

def ExportMesh(conn_name, identifier):
    if tau_settings["echo_level"] > 0:
        print "TAU SOLVER ExportMesh"
    # identifier is the data-name in json
    if identifier == "Fluid.Interface":
        nodal_coords, elem_connectivities, element_types = TauFunctions.GetFluidMesh(working_path, step, para_path_mod)
    else:
        raise Exception(
            'TauSolver::ExportMesh::identifier "{}" not valid! Please use Fluid.Interface'.format(identifier))

    CoSimIO.ExportMesh(conn_name, identifier, nodal_coords, elem_connectivities, element_types)
    if tau_settings["echo_level"] > 0:
        print "TAU SOLVER ExportMesh End"

connection_name = "TAU"

settings = {
    "echo_level" : "0",
    "print_timing" : "1",
    "communication_format" : "file"
}

CoSimIO.Connect(connection_name, settings)

CoSimIO.Register_AdvanceInTime(connection_name, AdvanceInTime)
CoSimIO.Register_InitializeSolutionStep(connection_name, InitializeSolutionStep)
CoSimIO.Register_SolveSolutionStep(connection_name, SolveSolutionStep)
CoSimIO.Register_FinalizeSolutionStep(connection_name, FinalizeSolutionStep)
CoSimIO.Register_ImportData(connection_name, ImportData)
CoSimIO.Register_ExportData(connection_name, ExportData)
CoSimIO.Register_ExportMesh(connection_name, ExportMesh)

# Run the coupled simulation, this returns after the entire CoSim is done
CoSimIO.Run(connection_name)

CoSimIO.Disconnect(connection_name)
tau("exit")
