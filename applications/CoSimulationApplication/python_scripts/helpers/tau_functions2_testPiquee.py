import shutil, sys, glob, os, time
#from tau_python import tau_plt_init_tecplot_params ####### import everything
from tau_python import *
tau_path = "/work/piquee/Softwares/TAU/TAU_2016.2/2016.2.0/bin/py_turb1eq"
sys.path.append(tau_path)
import tau_functions_testPiquee as tauFunctions
import PyPara, PyPrep, PySolv, PyDeform


#------------------------------------------------------------------
# find the solution file name in '/Outputs' and convert in .dat
#------------------------------------------------------------------
def findSolutionAndConvert(interface_file_path_pattern, mesh_file_path_pattern, this_step_out):
    list_of_interface_file_paths = glob.glob(interface_file_path_pattern + "*") ###### hier ich denke es war falsche... wir brauchen alle datei in Outputs sieht oben linie 21
    print "list_of_interface_file_path = ", list_of_interface_file_paths 

    interface_file_name = tauFunctions.findFileName(list_of_interface_file_paths, interface_file_path_pattern, "airfoilSol.pval.unsteady_i=",this_step_out) ###### hier ich denke this_step_out ist falsche
    print "interface_file_name = ", interface_file_name ######## first we need the mesh before looking for the number of lines

    list_of_meshes = glob.glob(mesh_file_path_pattern + "*") 
    if this_step_out == 0:
        mesh_iteration = tauFunctions.findFileName0(list_of_meshes, mesh_file_path_pattern,'airfoil_Structured_scaliert.grid')#### hier we need to separate .grid and .grid.def -> the this_step_out
    else:
        mesh_iteration = tauFunctions.findFileName(list_of_meshes, mesh_file_path_pattern,'airfoil_Structured_scaliert.grid', this_step_out)
    print "mesh_iteration = ", mesh_iteration

    tauFunctions.PrintBlockHeader("Start Writting Solution Data at time %s" %(str(time)))
    subprocess.call('rm ' +  working_path + '/Tautoplt.cntl' ,shell=True)
    tauFunctions.readTautoplt(working_path + '/Tautoplt.cntl', mesh_iteration, interface_file_name, para_path_mod)
    subprocess.call(TAU_path + '/tau2plt ' + working_path + '/Tautoplt.cntl' ,shell=True)
    tauFunctions.PrintBlockHeader("Stop Writting Solution Data at time %s" %(str(time)))

    if interface_file_name + '.dat' not in list_of_interface_file_paths:
        interface_file_name = interface_file_name[0:interface_file_name.find('+')]+ interface_file_name[interface_file_name.find('+')+1:len(interface_file_name)]
		
    interface_file_number_of_lines = tauFunctions.findInterfaceFileNumberOfLines(interface_file_name + '.dat')
    print 'interface_file_number_of_lines =', interface_file_number_of_lines
