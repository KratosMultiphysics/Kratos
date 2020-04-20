# -*- coding: utf-8 -*-
import shutil, sys, glob, os
import CoSimIO
# import modules for TAU
import PyPara, PyPrep, PySolv, PyDataSet
from tau_python import tau_plt_init_tecplot_params
# TODO Find a better way of indicating this script's path
this_scripts_path = "/work/piquee/Softwares/Kratos/applications/CoSimulationApplication/python_scripts/helpers"
sys.path.append(this_scripts_path)
from tau_functions import findInterfaceFileName
from tau_functions import findInterfaceFileNumberOfLines

print " Hello Test"

print "Hello world test"

##### Set up and initialize TAU #####
# Definition of the parameter file
para_path='airfoil_Structured.cntl'
para_path_mod = para_path + ".mod"
shutil.copy(para_path, para_path_mod)
working_path = "/work/piquee/MembraneWing/run_tau_from_kratos"
interface_file_path_pattern =  working_path + "/Outputs/*.plt"

# Initialize Tau python classes
Para = PyPara.Parafile(para_path_mod)
Prep = PyPrep.Preprocessing(para_path_mod)
Solver = PySolv.Solver(para_path_mod)
DataSetList = PyDataSet.DataSetList()

Prep.run(write_dualgrid=0,free_primgrid=False)
Solver.init(verbose = 1, reset_steps = True, n_time_steps = 1)

surfaces = ["MEMBRANE"]
DS = PyDataSet.DataSet(dataset_identifier = "name", output_functions = "surface",dataset_type ="surface",surf_def  = "name",surf_zone_list  = surfaces)
for i in range(0, len(surfaces)):
    DS.define_output(output_name = surfaces[i],
                    output_period     = 1,
                    output_variables  = ["cp"],
                    output_gather     = 1)

DataSetList.store_dataset(DS)

##### CoSimulation #####
def AdvanceInTime(current_time):
    print "TAU SOLVER AdvanceInTime"
    ts_tau = 0.1
    return 100.0#current_time + ts_tau

def InitializeSolutionStep():
    print("TAU SOLVER InitializeSolutionStep")

def SolveSolutionStep():
    print("TAU SOLVER SolveSolutionStep")
    Solver.outer_loop()
    Solver.output()
    tau_plt_init_tecplot_params(para_path_mod)
    tau_solver_write_output_conditional()
    DataSetList.write_output()
    list_of_interface_file_paths = glob.glob(interface_file_path_pattern)
    print list_of_interface_file_paths
    interface_file_name = findInterfaceFileName(
        list_of_interface_file_paths, working_path, this_step_out=0)
    print interface_file_name
    interface_file_number_of_lines = findInterfaceFileNumberOfLines(interface_file_name)
    print interface_file_number_of_lines

def FinalizeSolutionStep():
    print("TAU SOLVER FinalizeSolutionStep")

def ImportData(conn_name, identifier):
    print "TAU SOLVER ImportData"
    #data = CoSimIO.ImportData(conn_name, identifier)
    print "TAU SOLVER After ImportData"

    # TODO do sth with the data
    # identifier is the data-name in json
    if identifier == "displacements":
        pass
    else:
        raise Exception

def ExportData(conn_name, identifier):
    print "TAU SOLVER ExportData"
    # identifier is the data-name in json
    if identifier == "forces":
        data = GetFluidForces()
    else:
        raise Exception

    CoSimIO.ExportData(conn_name, identifier, data)
    print "TAU SOLVER After ExportData"

def ExportMesh(conn_name, identifier):
    # identifier is the data-name in json
    if identifier == "wing_fsi_interface":
        nodal_coords, elem_connectivities, element_types = GetFluidMesh()
    else:
        raise Exception

    CoSimIO.ExportMesh(conn_name, identifier, nodal_coords, elem_connectivities, element_types)


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

# Run the coupled simulation
print "Before Run"
CoSimIO.Run(connection_name) #this returns after the entire CoSim is done

CoSimIO.Disconnect(connection_name)

DataSetList.free_data()
Solver.finalize()
Para.free_parameters()
tau("exit")
