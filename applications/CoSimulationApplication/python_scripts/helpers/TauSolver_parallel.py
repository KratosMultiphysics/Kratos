# -*- coding: utf-8 -*-
import shutil, sys, os, time, json
from datetime import datetime
import numpy as np
import CoSimIO
import subprocess
sys.path.append("/work/piquee/Softwares/TAU/TAU_2016.2/2016.2.0/bin/py_turb1eq/")
import PyPara, PyPrep, PySolv, PyDeform, PyCopyCluster, PyMotionExternalDelegate
from tau_python import tau_solver_unsteady_get_physical_time
from tau_python import *
rank = tau_mpi_rank()

print rank

with open('/work/piquee/MembraneWing/kratos_fsi_big/tau_settings.json') as json_file:
    tau_settings = json.load(json_file)
start_step = tau_settings["start_step"]
tau_path = tau_settings["tau_path"]
sys.path.append(tau_settings["kratos_path"])
sys.path.append(tau_path + "py_turb1eq/")
working_path = os.getcwd() + '/'
rotate = tau_settings["rotate"]
tau_parallel_sync()
NCores='3'
mpirun='/work/piquee/Softwares/openmpi_1.6.4/bin/mpirun'

# tau_functions can only be imported after appending kratos' path
sys.path.append("/work/piquee/Softwares/Kratos/applications/CoSimulationApplication/python_scripts/helpers/")
import tau_functions as TauFunctions

# Remove output files and deform mesh files from previous simulations
if (rank == 0):
    TauFunctions.RemoveFilesFromPreviousSimulations()
tau_parallel_sync()

# Definition of the parameter file
if rotate:
    para_path='airfoil_Structured_rotation.cntl'
else:
    para_path='airfoil_Structured.cntl'
tau_parallel_sync()
para_path_mod = para_path + ".mod"
shutil.copy(para_path, para_path_mod)
tau_parallel_sync()

# Initialize Tau python classes and auxiliary variable step
Para = PyPara.Parafile(para_path_mod)
Prep = PyPrep.Preprocessing(para_path_mod)


tau_time_step = float(Para.get_para_value('Unsteady physical time step size'))
tau_parallel_sync()

if rotate:
    # --- Prepare parameters for unsteady simulation ---
    Para.update({"Unsteady allow external control over progress in time (0/1)": 1,
                 "Unsteady enable external control over progress in time (0/1)": 1})

    # external output period (needed if output period in para_path not working properly)
    pitchDeg = 0  # mean pitch angle

    # --- Load external excitation files ---
    # --- thetaDeg -> pitch angle per timestep in deg
    # --- thetaRate -> pitch rate per timestep in deg/s
    thetaDeg = np.loadtxt(working_path + 'signal/APRBSDeg_membrane.dat')
    thetaRate = np.loadtxt(working_path + 'signal/APRBSRate_membrane.dat')

    MyMotionStringGenerator = TauFunctions.MotionStringGenerator(
        tau_time_step, pitchDeg, thetaDeg, thetaRate)

    # --- Prepare parameters for unsteady rotating simulation ---
    nodeNames = ["MEMBRANE"]

    # --- Instanciate required modules ---
    MyTauMotionDelegate = PyMotionExternalDelegate.MotionExternalDelegate(
        nodeNames)

    Solver = PySolv.Solver(para_path_mod, delegate=MyTauMotionDelegate)

    # --- Init motion stack ---
    motionString = MyMotionStringGenerator(0)
    MyTauMotionDelegate.UpdateMotion("MEMBRANE", motionString)
    MyTauMotionDelegate.PushMotion()
    TauFunctions.PrintBlockHeader("Inital Motionstring: %s" %(motionString))
else:
    Solver = PySolv.Solver(para_path_mod)

Deform = PyDeform.Deformation(para_path_mod)
step = start_step
Prep.run(write_dualgrid=True,free_primgrid=False)
Solver.init(verbose = 1, reset_steps = True, n_time_steps = 1)
#------------------------------------------------------------------
# Start Partitioning
#-------------------------------------------------------------------
#if (rank == 0):
#  print 'Partitioning start', str(datetime.now())	
#  print(mpirun + ' -np ' + '%s '%NCores + tau_path + 'ptau3d.subgrids ' + para_path_mod + ' out.partitioning use_mpi')
#  subprocess.call(mpirun + ' -np ' + '%s '%NCores + tau_path + 'ptau3d.subgrids ' + para_path_mod + ' out.partitioning use_mpi',shell=True)
#  subprocess.call('sleep 5' ,shell=True)
#  print 'Partitioning end'
#tau_parallel_sync()

#------------------------------------------------------------------
# End Partitioning
#-------------------------------------------------------------------

def AdvanceInTime(current_time):
    # Preprocessing needs to be done before getting the time and time step
    TauFunctions.PrintBlockHeader("Start Preprocessing at time %s" %(str(time)))
    print "rank bei advance in time / preprocessing: " + str(rank)
    #Prep.run(write_dualgrid=True,free_primgrid=False)
    tau_parallel_sync()
    print "rank bei advance in time / solving: " + str(rank)
    TauFunctions.PrintBlockHeader("Stop Preprocessing at time %s" %(str(time)))
    TauFunctions.PrintBlockHeader("Initialize Solver at time %s" %(str(time)))
   # Solver.init(verbose = 1, reset_steps = True, n_time_steps = 1)
    tau_parallel_sync() 

    if rotate:
        motionString = MyMotionStringGenerator(step - start_step)
        MyTauMotionDelegate.UpdateMotion("MEMBRANE", motionString)
        MyTauMotionDelegate.InitExchange()

    # Get current time and time step from tau
    tau_current_time = float(tau_solver_unsteady_get_physical_time())

    if tau_settings["echo_level"] > 0:
        print "TAU SOLVER AdvanceInTime"
        print 'tau_current_time = ', tau_current_time
        print 'tau_time_step = ', tau_time_step
    return tau_current_time + tau_time_step


def InitializeSolutionStep():
    if tau_settings["echo_level"] > 0:
        if (rank == 0):
            print("TAU SOLVER InitializeSolutionStep")
        tau_parallel_sync() 

def SolveSolutionStep():
    if tau_settings["echo_level"] > 0:
        if (rank == 0):
            print("TAU SOLVER SolveSolutionStep")
        tau_parallel_sync() 
      
    Solver.outer_loop()
    tau_parallel_sync() 
    Solver.output()
    tau_parallel_sync() 
    tau_plt_init_tecplot_params(para_path_mod)
    tau_solver_write_output_conditional()
    # Convert tau output to dat file using tau2plt
    TauFunctions.ConvertOutputToDat(working_path, tau_path, step, para_path_mod, start_step)

def FinalizeSolutionStep():
    if tau_settings["echo_level"] > 0:
        if (rank == 0):
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
        if (rank == 0):
            print "TAU SOLVER ImportData"
        tau_parallel_sync() 
    displacements = CoSimIO.ImportData(conn_name, identifier)

    # identifier is the data-name in json
    if identifier == "Interface_disp":
        TauFunctions.ExecuteBeforeMeshDeformation(displacements, step, para_path_mod, start_step)
        # Deform mesh
        Deform.run(read_primgrid=1, write_primgrid=1, read_deformation=0, field_io=1)
    else:
        raise Exception('TauSolver::ExportData::identifier "{}" not valid! Please use Interface_disp'.format(identifier))
    if tau_settings["echo_level"] > 0:
        if (rank == 0):
            print 'maximum_displacement_kratos = ', max(displacements)
            print "TAU SOLVER After ImportData"
        tau_parallel_sync() 

def ExportData(conn_name, identifier):
    if tau_settings["echo_level"] > 0:
        if (rank == 0):
            print "TAU SOLVER ExportData"
        tau_parallel_sync() 
    # identifier is the data-name in json
    if identifier == "Interface_force":
        if (rank == 0):
            forces = TauFunctions.ComputeFluidForces(working_path, step)
        tau_parallel_sync() 
    else:
        raise Exception('TauSolver::ExportData::identifier "{}" not valid! Please use Interface_force'.format(identifier))

    CoSimIO.ExportData(conn_name, identifier, forces)
    if tau_settings["echo_level"] > 0:
        if (rank == 0):
            print 'maximum_force_tau = ', max(forces)
            print "TAU SOLVER After ExportData"
        tau_parallel_sync() 

def ExportMesh(conn_name, identifier):
    if tau_settings["echo_level"] > 0:
        if (rank == 0):
            print "TAU SOLVER ExportMesh"
        tau_parallel_sync()
    # identifier is the data-name in json
    if identifier == "Fluid.Interface":
        if (rank == 0):
            nodal_coords, elem_connectivities, element_types = TauFunctions.GetFluidMesh(working_path, step, para_path_mod)
        tau_parallel_sync() 
    else:
        raise Exception(
            'TauSolver::ExportMesh::identifier "{}" not valid! Please use Fluid.Interface'.format(identifier))
    
    if (rank == 0):
        CoSimIO.ExportMesh(conn_name, identifier, nodal_coords, elem_connectivities, element_types)
    tau_parallel_sync() 
    if tau_settings["echo_level"] > 0:
        if (rank == 0):
            print "TAU SOLVER ExportMesh End"
        tau_parallel_sync() 

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
