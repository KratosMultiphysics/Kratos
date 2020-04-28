# -*- coding: utf-8 -*-
import shutil, sys, glob, os, time, subprocess, json
import numpy as np
import CoSimIO
import PyPara, PyPrep, PySolv, PyDeform, PySurfDeflect

with open('tau_settings.json') as json_file:
    tau_settings = json.load(json_file)

TAU_path = tau_settings["tau_path"]
sys.path.append(tau_settings["kratos_path"])
sys.path.append(TAU_path + "py_turb1eq/")
working_path = os.getcwd() + '/'
interface_file_path_pattern =  working_path + "Outputs/"
mesh_file_path_pattern = working_path + "Mesh/"

# tau_functions can only be imported after appending kratos' path
import tau_functions as tauFunctions

# Definition of the parameter file
para_path='airfoil_Structured.cntl'
para_path_mod = para_path + ".mod"
shutil.copy(para_path, para_path_mod)

# Initialize Tau python classes
Para = PyPara.Parafile(para_path_mod)
Prep = PyPrep.Preprocessing(para_path_mod)
Solver = PySolv.Solver(para_path_mod)
Deform = PyDeform.Deformation(para_path_mod)

# Initialize variables
this_step_out = 0
NodesNr = 0

def AdvanceInTime(current_time):
    # Preprocessing needs to be done before getting the time and time step
    tauFunctions.PrintBlockHeader("Start Preprocessing at time %s" %(str(time)))
    Prep.run(write_dualgrid=False,free_primgrid=False) 
    tauFunctions.PrintBlockHeader("Stop Preprocessing at time %s" %(str(time)))
    tauFunctions.PrintBlockHeader("Initialize Solver at time %s" %(str(time)))
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

def FinalizeSolutionStep():
    if tau_settings["echo_level"] > 0:
        print("TAU SOLVER FinalizeSolutionStep")
    global this_step_out
    tau_parallel_sync()
    Solver.finalize()
    tau_free_dualgrid()
    tau_free_prims()
    Para.free_parameters()
    this_step_out += 1

def ImportData(conn_name, identifier):
    if tau_settings["echo_level"] > 0:
        print "TAU SOLVER ImportData"
    data = CoSimIO.ImportData(conn_name, identifier)

    # TODO do sth with the data
    # identifier is the data-name in json
    if identifier == "Interface_disp":
        deformMesh(data)
    else:
        raise Exception('TauSolver::ExportData::identifier "{}" not valid! Please use Interface_disp'.format(identifier))
    if tau_settings["echo_level"] > 0:
        print "TAU SOLVER After ImportData"

def ExportData(conn_name, identifier):
    if tau_settings["echo_level"] > 0:
        print "TAU SOLVER ExportData"
    # identifier is the data-name in json
    if identifier == "Interface_force":
        interface_file_name_surface, interface_file_number_of_lines = tauFunctions.findSolutionAndConvert(interface_file_path_pattern, mesh_file_path_pattern, this_step_out, para_path_mod)
        data = caculatePressure(interface_file_name_surface, interface_file_number_of_lines, this_step_out)
    else:
        raise Exception('TauSolver::ExportData::identifier "{}" not valid! Please use Interface_force'.format(identifier))

    CoSimIO.ExportData(conn_name, identifier, data)
    if tau_settings["echo_level"] > 0:
        print "TAU SOLVER After ExportData"


def ExportMesh(conn_name, identifier):
    if tau_settings["echo_level"] > 0:
        print "TAU SOLVER ExportMesh"
        print "conn_name = ", conn_name
        print "identifier = ", identifier
    # identifier is the data-name in json
    if identifier == "Fluid.Interface":
        interface_file_name_surface, interface_file_number_of_lines = tauFunctions.findSolutionAndConvert(
            interface_file_path_pattern, mesh_file_path_pattern, this_step_out, para_path_mod)
        NodesNr, ElemsNr, X, Y, Z, CP, P, elemTable, liste_number = tauFunctions.readPressure(interface_file_name_surface + '.dat', interface_file_number_of_lines, 0, 20)
        nodal_coords, nodesID, elem_connectivities, element_types = tauFunctions.interfaceMeshFluid(NodesNr, ElemsNr, elemTable, X, Y, Z)
        # In vtk format element connectivities start from 0, not from 1
        elem_connectivities -= 1
        # nodal_coords, elem_connectivities, element_types = GetFluidMesh()
    else:
        raise Exception(
            'TauSolver::ExportMesh::identifier "{}" not valid! Please use Fluid.Interface'.format(identifier))

    CoSimIO.ExportMesh(conn_name, identifier, nodal_coords, elem_connectivities, element_types)
    if tau_settings["echo_level"] > 0:
        print "TAU SOLVER ExportMesh End"

# read the solution file name in '/Outputs' and calculate the pressure
def caculatePressure(interface_file_name_surface, interface_file_number_of_lines, this_step_out):
    global NodesNr
    NodesNr,ElemsNr,X,Y,Z,CP,P,elemTable,liste_number=tauFunctions.readPressure( interface_file_name_surface + '.dat', interface_file_number_of_lines, 0, 20)
    print 'NodesNr = ', NodesNr

    nodes,nodesID,elems,element_types=tauFunctions.interfaceMeshFluid(NodesNr,ElemsNr,elemTable,X,Y,Z)
  
    # calculating cp at the center of each interface element    
    pCell=tauFunctions.calcpCell(ElemsNr,P,X,elemTable)

    # calculating element area and normal vector
    area,normal = tauFunctions.calcAreaNormal(ElemsNr,elemTable,X,Y,Z,(this_step_out+1))

    # calculating the force vector
    forcesTauNP = tauFunctions.calcFluidForceVector(ElemsNr,elemTable,NodesNr,pCell,area,normal,(this_step_out+1))

    return forcesTauNP

#------------------------------------------------------------------
# Deformation
#------------------------------------------------------------------
def deformMesh(dispTau):
    global dispTauOld

    if tau_mpi_rank() == 0:
        print "deformationstart"
        interface_file_name_surface, interface_file_number_of_lines = tauFunctions.findSolutionAndConvert(interface_file_path_pattern, mesh_file_path_pattern, this_step_out, para_path_mod)

        NodesNr,ElemsNr,X,Y,Z,CP,P,elemTable,liste_number=tauFunctions.readPressure( interface_file_name_surface + '.dat', interface_file_number_of_lines, 0, 20)

        nodes,nodesID,elems,element_types=tauFunctions.interfaceMeshFluid(NodesNr,ElemsNr,elemTable,X,Y,Z)

        if(this_step_out==0):
            dispTauOld=np.zeros(3*NodesNr)
            dispTau_transpose = np.transpose(dispTau)
            print 'dispTau =', dispTau_transpose

        [ids,coordinates,globalID,coords]=tauFunctions.meshDeformation(NodesNr,nodes,dispTau,dispTauOld,0,para_path_mod)
        PySurfDeflect.write_test_surface_file('deformation_file',coords[:,0:2],coords[:,3:5])
        print "afterPySurfDeflect"


    Deform.run(read_primgrid=1, write_primgrid=1, read_deformation=0, field_io=1) 


    for i in xrange(0,3*NodesNr):
        dispTauOld[i]=dispTau[i]
    print "afterDeformation"

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