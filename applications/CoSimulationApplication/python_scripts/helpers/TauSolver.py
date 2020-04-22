# -*- coding: utf-8 -*-
import shutil, sys, glob, os, time, subprocess
import CoSimIO
# import modules for TAU
import PyPara, PyPrep, PySolv, PyDeform
#from tau_python import tau_plt_init_tecplot_params ####### import everything
sys.path.append("/work/piquee/Softwares/TAU/TAU_2016.2/2016.2.0/bin/py_turb1eq")
import tau_python 
from tau_python import tau_msg
# TODO Find a better way of indicating this script's path
this_scripts_path = "/work/piquee/Softwares/Kratos/applications/CoSimulationApplication/python_scripts/helpers/"
sys.path.append(this_scripts_path)

import tau_functions as tauFunctions

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
# Initialize TimeStep
this_step_out = 0 ### TO do  read from TAU
# Start Time
startTime = 0 ### TO do  read from TAU

#------------------------------------------------------------------
# Convert the initial TAU Mesh in '.dat' to find the information necessary  ######### TO DO ######### maybe at first !!!!!!
#-------------------------------------------------------------------
tauFunctions.PrintBlockHeader("Initial TAU Mesh at time %s" %(str(time)))



##### CoSimulation #####
def AdvanceInTime(current_time):
    print "TAU SOLVER AdvanceInTime"
    ts_tau = 0.005 ### TO do  read from TAU
    return 100.0#current_time + ts_tau

#------------------------------------------------------------------
#         ############  ich glaube es muss in ein Schleife sein - für jede Schritte muss mann es machen #########
#-------------------------------------------------------------------
#def InitializePreprocessingStep():


#------------------------------------------------------------------
#  Preprocessing und Solving Initialize       ############  ich glaube es muss in ein Schleife sein - für jede Schritte muss mann es machen #########
#-------------------------------------------------------------------
def InitializeSolutionStep():
    tauFunctions.PrintBlockHeader("Start Preprocessing at time %s" %(str(time)))
    Prep.run(write_dualgrid=False,free_primgrid=False) 
    tauFunctions.PrintBlockHeader("Stop Preprocessing at time %s" %(str(time)))
    tauFunctions.PrintBlockHeader("Initialize Solver at time %s" %(str(time)))
    Solver.init(verbose = 1, reset_steps = True, n_time_steps = 1)
    print("TAU SOLVER InitializeSolutionStep")

#------------------------------------------------------------------
# Solving Solve       ############  ich glaube es muss in ein Schleife sein - für jede Schritte muss mann es machen #########
#-------------------------------------------------------------------
def SolveSolutionStep():
    print("TAU SOLVER SolveSolutionStep")
    Solver.outer_loop()
    Solver.output()
    tau_plt_init_tecplot_params(para_path_mod)
    tau_solver_write_output_conditional()
    interface_file_name_surface, interface_file_number_of_lines = findSolutionAndConvert(interface_file_path_pattern, mesh_file_path_pattern, this_step_out)
    forcesTauNP = caculatePressure(interface_file_name_surface, interface_file_number_of_lines, this_step_out)

#------------------------------------------------------------------
# find the solution file name in '/Outputs' and convert in .dat ############  ich glaube es muss in ein Schleife sein - für jede Schritte muss mann es machen
#------------------------------------------------------------------
def findSolutionAndConvert(interface_file_path_pattern, mesh_file_path_pattern, this_step_out):
    list_of_interface_file_paths = glob.glob(interface_file_path_pattern + "*") 
    print "list_of_interface_file_path = ", list_of_interface_file_paths 

    interface_file_name = tauFunctions.findFileName(list_of_interface_file_paths, interface_file_path_pattern, "airfoilSol.pval.unsteady_i=",this_step_out+1) 
    print "interface_file_name = ", interface_file_name 

    list_of_meshes = glob.glob(mesh_file_path_pattern+ "*") 
    if this_step_out == 0:
        mesh_iteration = tauFunctions.findFileName0(list_of_meshes, mesh_file_path_pattern,'airfoil_Structured_scaliert.grid')
    else:
        mesh_iteration = tauFunctions.findFileName(list_of_meshes, mesh_file_path_pattern,'airfoil_Structured_scaliert.grid', this_step_out)
    print "mesh_iteration = ", mesh_iteration

    tauFunctions.PrintBlockHeader("Start Writting Solution Data at time %s" %(str(time)))
    subprocess.call('rm ' +  working_path + '/Tautoplt.cntl' ,shell=True)
    tauFunctions.readTautoplt(working_path + 'Tautoplt.cntl', working_path + 'Tautoplt_initial.cntl', mesh_iteration, interface_file_name,working_path + para_path_mod)
    subprocess.call(TAU_path + 'tau2plt ' + working_path + '/Tautoplt.cntl' ,shell=True)
    tauFunctions.PrintBlockHeader("Stop Writting Solution Data at time %s" %(str(time)))

    if interface_file_name + '.dat' not in glob.glob(interface_file_path_pattern + "*"):
        interface_file_name = interface_file_name[0:interface_file_name.find('+')]+ interface_file_name[interface_file_name.find('+')+1:len(interface_file_name)]
    
    if 'surface' not in interface_file_name + '.dat':
        interface_file_name_surface = interface_file_name[0:interface_file_name.find('.pval')]+ '.surface.' + interface_file_name[interface_file_name.find('.pval')+1:len(interface_file_name)]
    else: 
        interface_file_name_surface = interface_file_name

    interface_file_number_of_lines = tauFunctions.findInterfaceFileNumberOfLines(interface_file_name_surface + '.dat')
    print 'interface_file_number_of_lines =', interface_file_number_of_lines

    return interface_file_name_surface, interface_file_number_of_lines

#------------------------------------------------------------------
# read the solution file name in '/Outputs' and calculate the pressure ## ich glaube es muss in ein Schleife sein - für jede Schritte muss mann es machen ##
#------------------------------------------------------------------
def caculatePressure(interface_file_name_surface, interface_file_number_of_lines, this_step_out):
    NodesNr,ElemsNr,X,Y,Z,CP,P,elemTable_Sol,liste_number=tauFunctions.readPressure( interface_file_name_surface + '.dat', interface_file_number_of_lines, 0, 20)

    elemTable = elemTable_Sol.astype(int)

    with open('xpNode','w') as f:
        for i in xrange(0,NodesNr):											
            f.write('%d\t%f\t%f\t%f\t%f\n'%(i,X[i],Y[i],Z[i],P[i]))	

    nodes,nodesID,elems,numNodesPerElem=tauFunctions.interfaceMeshFluid(NodesNr,ElemsNr,elemTable,X,Y,Z)
  
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
#def deformationTAU():

#    [ids,coordinates,globalID,coords]=Functions.meshDeformation(NodesNr,nodes,dispTau,dispTauOld,0)
 #   PySurfDeflect.write_test_surface_file('deformation_file',coords[:,0:2],coords[:,3:5])


#------------------------------------------------------------------
# 
#------------------------------------------------------------------
  #if (this_step_out==0):  
   #   TAUclient.setMesh('myMeshTau', NodesNr, ElemsNr, nodes, nodesID, numNodesPerElem, elems, 'TAUclient')
 #  


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

Solver.finalize()
Para.free_parameters()
tau("exit")
