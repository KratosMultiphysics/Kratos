# -*- coding: utf-8 -*-
import shutil, sys, glob, os, time, subprocess
import CoSimIO
# import modules for TAU
import PyPara, PyPrep, PySolv, PyDeform
#from tau_python import tau_plt_init_tecplot_params ####### import everything
from tau_python import *
# TODO Find a better way of indicating this script's path
this_scripts_path = "/work/piquee/Softwares/Kratos/applications/CoSimulationApplication/python_scripts/helpers/"
sys.path.append(this_scripts_path)
import tau_functions_testPiquee as tauFunctions

#-------------------------------------------------------------------------------
# Definitions
#-------------------------------------------------------------------------------
# Definition of the parameter file
para_path='airfoil_Structured.cntl'
para_path_mod = para_path + ".mod"
shutil.copy(para_path, para_path_mod)
working_path = "/work/piquee/MembraneWing/run_tau_from_kratos/"
interface_file_path_pattern =  working_path + "Outputs/"
mesh_file_path_pattern = working_path + "Mesh/"

# Definition of TAU path
TAU_path = "/work/piquee/Softwares/TAU/TAU_2016.2/2016.2.0/bin/"

#-------------------------------------------------------------------------------
# Init Tau python classes + get the informations necessary for the preprocessing 
#-------------------------------------------------------------------------------
Para = PyPara.Parafile(para_path_mod)
Prep = PyPrep.Preprocessing(para_path_mod)
Solver = PySolv.Solver(para_path_mod)
Deform = PyDeform.Deformation(para_path_mod)

# Definitions
# Primary grid filename
grid = Para.get_para_value("Primary grid filename") # 
# Surfaces where Cp has to be written 
surfaces = ["MEMBRANE"] # 
# get the name of the file where the deformation of the TAU Mesh are saved
deformfile = Para.get_para_value('RBF basis coordinates and deflections filename') # 

#------------------------------------------------------------------
# Convert the initial TAU Mesh in '.dat' to find the information necessary  ######### TO DO ######### maybe at first !!!!!!
#-------------------------------------------------------------------
tauFunctions.PrintBlockHeader("Initial TAU Mesh at time %d" %(str(time)))



##### CoSimulation #####
def AdvanceInTime(current_time):
    print "TAU SOLVER AdvanceInTime"
    ts_tau = 0.1
    return 100.0#current_time + ts_tau

#------------------------------------------------------------------
# Preprocessing            ############  ich glaube es muss in ein Schleife sein - für jede Schritte muss mann es machen #########
#-------------------------------------------------------------------
def InitializePreprocessingStep():
    tauFunctions.PrintBlockHeader("Start Preprocessing at time %d" %(str(time)))
    Prep.run(write_dualgrid=False,free_primgrid=False) 
    tauFunctions.PrintBlockHeader("Stop Preprocessing at time %d" %(str(time)))

#------------------------------------------------------------------
# Solving Initialize       ############  ich glaube es muss in ein Schleife sein - für jede Schritte muss mann es machen #########
#-------------------------------------------------------------------
def InitializeSolutionStep():
    tauFunctions.PrintBlockHeader("Initialize Solver at time %d" %(str(time)))
    Solver.init(verbose = 1, reset_steps = True, n_time_steps = 1)
    print("TAU SOLVER InitializeSolutionStep")

#------------------------------------------------------------------
# Solving Solve       ############  ich glaube es muss in ein Schleife sein - für jede Schritte muss mann es machen #########
#-------------------------------------------------------------------
def SolveSolutionStep(para_path_mod):
    print("TAU SOLVER SolveSolutionStep")
    Solver.outer_loop()
    Solver.output()
    tau_plt_init_tecplot_params(para_path_mod)
    tau_solver_write_output_conditional()

#------------------------------------------------------------------
# find the solution file name in '/Outputs' and convert in .dat ############  ich glaube es muss in ein Schleife sein - für jede Schritte muss mann es machen #########
#------------------------------------------------------------------
def findSolutionAndConvert(interface_file_path_pattern, mesh_file_path_pattern, this_step_out)
    list_of_interface_file_paths = glob.glob(interface_file_path_pattern + "*") 
    print "list_of_interface_file_path = ", list_of_interface_file_paths 

    interface_file_name = tauFunctions.findFileName(list_of_interface_file_paths, interface_file_path_pattern, "airfoilSol.pval.unsteady_i=",this_step_out) ###### hier ich denke this_step_out ist falsche
    print "interface_file_name = ", interface_file_name 

    list_of_meshes = glob.glob(mesh_file_path_pattern+ "*") 
    if this_step_out == 0:
        mesh_iteration = tauFunctions.findFileName0(list_of_meshes, mesh_file_path_pattern,'airfoil_Structured_scaliert.grid')
    else:
        mesh_iteration = tauFunctions.findFileName(list_of_meshes, mesh_file_path_pattern,'airfoil_Structured_scaliert.grid', this_step_out)
    print "mesh_iteration = ", mesh_iteration

    tauFunctions.PrintBlockHeader("Start Writting Solution Data at time %d" %(str(time)))
    subprocess.call('rm ' +  working_path + '/Tautoplt.cntl' ,shell=True)
    tauFunctions.readTautoplt(working_path + 'Tautoplt.cntl', working_path + 'Tautoplt_initial.cntl', mesh_iteration, interface_file_name,working_path + para_path_mod)
    subprocess.call(TAU_path + 'tau2plt ' + working_path + '/Tautoplt.cntl' ,shell=True)
    tauFunctions.PrintBlockHeader("Stop Writting Solution Data at time %d" %(str(time)))

    if interface_file_name + '.dat' not in glob.glob(interface_file_path_pattern + "*"):
        interface_file_name = interface_file_name[0:interface_file_name.find('+')]+ interface_file_name[interface_file_name.find('+')+1:len(interface_file_name)]
    
    if 'surface' not in interface_file_name + '.dat':
        interface_file_name_surface = interface_file_name[0:interface_file_name.find('.pval')]+ 'surface' + interface_file_name[interface_file_name.find('.pval')+1:len(interface_file_name)]
    else: 
        interface_file_name_surface = interface_file_name

    interface_file_number_of_lines = tauFunctions.findInterfaceFileNumberOfLines(interface_file_name_surface + '.dat')
    print 'interface_file_number_of_lines =', interface_file_number_of_lines

#------------------------------------------------------------------
# read the solution file name in '/Outputs' and calculate the pressure ############  ich glaube es muss in ein Schleife sein - für jede Schritte muss mann es machen #########
#------------------------------------------------------------------
def caculatePressure():
    NodesNr,ElemsNr,X,Y,Z,CP,P,elemTable_Sol,liste_number=tauFunctions.readPressure( interface_file_name_surface + '.dat', interface_file_number_of_lines, 0)
    elemTable = elemTable_Sol.astype(int)

    with open('xpNode','w') as f:
        for i in xrange(0,NodesNr_Sol):											
            f.write('%d\t%f\t%f\t%f\t%f\n'%(i,X[i],Y[i],Z[i],P[i]))	

    nodes,nodesID,elems,numNodesPerElem=tauFunctions.interfaceMeshFluid(NodesNr,ElemsNr,elemTable,X,Y,Z)
  
  #if (this_step_out==0):  
   #   TAUclient.setMesh('myMeshTau', NodesNr, ElemsNr, nodes, nodesID, numNodesPerElem, elems, 'TAUclient')
 #     

  # calculating cp at the center of each interface element    
 # pCell=Functions.calcpCell(ElemsNr,P,X,elemTable)
 # print 'calcCpCell'
  # calculating element area and normal vector
 # area,normal = Functions.calcAreaNormal(ElemsNr,elemTable,X,Y,Z,fluidIter*(this_step_out+1))
 # print 'calcAreaNormal'
  # calculating the force vector
  #forcesTauNP = Functions.calcFluidForceVector(ElemsNr,elemTable,NodesNr,pCell,area,normal,fluidIter*(this_step_out+1))



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
CoSimIO.Register_InitializePreprocessingStep()
CoSimIO.Register_InitializeSolutionStep()
#########   CoSimIO.Register_InitializeSolutionStep(connection_name, InitializeSolutionStep) ######### welche
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
