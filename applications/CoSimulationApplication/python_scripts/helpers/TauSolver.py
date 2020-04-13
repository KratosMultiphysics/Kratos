# import modules for TAU
# import PyPara, PyPrep, PySolv
import shutil, sys

import CoSimIO

##### Set up and initialize TAU #####

Definition of the parameter file
para_path='airfoil_Structured.cntl'
para_path_mod = para_path + ".mod"
shutil.copy(para_path, para_path_mod)

# Init Tau python classes
Para = PyPara.Parafile(para_path_mod)
Prep = PyPrep.Preprocessing(para_path_mod)
Solver = PySolv.Solver(para_path_mod)

Prep.run(write_dualgrid=0,free_primgrid=False)
Solver.init(verbose = 1, reset_steps = True, n_time_steps = 1)


##### CoSimulation #####

def AdvanceInTime(current_time):
    ts_tau = 0.1
    return current_time + ts_tau

def InitializeSolutionStep():
    pass

def SolveSolutionStep():
    pass

def FinalizeSolutionStep():
    Solver.outer_loop()
    Solver.output()


def ImportData(conn_name, identifier):
    data = CoSimIO.ImportData(conn_name, identifier)

    # TODO do sth with the data
    # identifier is the data-name in json
    if identifier == "displacements":
        pass
    else:
        raise Exception

def ExportData(conn_name, identifier):
    # identifier is the data-name in json
    if identifier == "forces":
        data = GetFluidForces()
    else:
        raise Exception

    CoSimIO.ExportData(conn_name, identifier, data)

def ExportMesh(conn_name, identifier):
    # identifier is the data-name in json
    if identifier == "wing_fsi_interface":
        nodal_coords, elem_connectivities, element_types = GetFluidMesh()
    else:
        raise Exception

    CoSimIO.ExportMesh(conn_name, identifier, nodal_coords, elem_connectivities, element_types)


connection_name = "TAU"

CoSimIO.Connect(connection_name, "tau_so_sim_io_settings.txt") # TODO @Philipp in the future this can also be a dict

CoSimIO.Register_AdvanceInTime(connection_name, AdvanceInTime)
CoSimIO.Register_InitializeSolutionStep(connection_name, InitializeSolutionStep)
CoSimIO.Register_SolveSolutionStep(connection_name, SolveSolutionStep)
CoSimIO.Register_FinalizeSolutionStep(connection_name, FinalizeSolutionStep)

CoSimIO.Register_ImportData(connection_name, ImportData)
CoSimIO.Register_ExportData(connection_name, ExportData)
CoSimIO.Register_ExportMesh(connection_name, ExportMesh)

# Run the coupled simulation
CoSimIO.Run(connection_name) #this returns after the entire CoSim is done

CoSimIO.Disconnect(connection_name)

###############################################################################
######## ALTERNATIVE: WRITE COUPLING LOOP EXPLICITLY (NOT RECOMMENDED) ########
###############################################################################
'''
CoSimIO.Connect(connection_name, "tau_so_sim_io_settings.txt")

ExportMesh(connection_name, "wing_fsi_interface")

while current_time < end_time:
    current_time = AdvanceInTime(current_time)
    InitializeSolutionStep()

    while True:
        ImportData(connection_name, "displacements")
        SolveSolutionStep()
        ExportData(connection_name, "forces")
        if CoSimIO.IsConverged(): break

    FinalizeSolutionStep()

CoSimIO.Disconnect(connection_name)
'''


Solver.finalize()
Para.free_parameters()
tau("exit")
