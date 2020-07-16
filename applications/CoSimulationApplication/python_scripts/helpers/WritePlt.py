import os
import sys
import math
import re, glob, subprocess, time, os, warnings
###############
# Definitions
###############
cwd = os.getcwd()
sys.path.append(cwd)
working_path = cwd
n_step = 150
para_path_mod = working_path + '/airfoil_Structured.cntl'
start_step = 0
primary_grid_filename = working_path + '/Mesh/NeueGeo_WKA2D_MembraneSeparate_scaliert.grid'
tau_path = '/work/piquee/Softwares/TAU/taudir_repos.2019.13.08/bin'
rotate = 0

def WriteTautoplt(working_path, step, para_path_mod, start_step, primary_grid_filename):
    # Define Tautoplt.cntl file name and check if it already exists
    tautoplt_filename = working_path + '/Tautoplt.cntl'
    RemoveFileIfExists(tautoplt_filename)
    initial_tautoplt_filename = working_path + '/Tautoplt_initial.cntl'

    # Read and write simultaneously
    tautoplt_file_writing = open(tautoplt_filename, 'w')
    with open(initial_tautoplt_filename, 'r+') as tautoplt_file_reading:
        # Loop over lines
        for line in tautoplt_file_reading:
            # Check if line is from IO section and modify it
            line = ModifyFilesIOLines(line, working_path, step, para_path_mod, start_step, primary_grid_filename)
            tautoplt_file_writing.write(line)

        # Close files
        tautoplt_file_writing.close()
        tautoplt_file_reading.close()

    return tautoplt_filename

def ModifyFilesIOLines(line, working_path, step, para_path_mod, start_step, primary_grid_filename):
    if 'Primary grid filename:' in line:
        line = 'Primary grid filename:' + primary_grid_filename + ' \n'
    elif 'Boundary mapping filename:' in line:
        parameter_filename = para_path_mod
        line = 'Boundary mapping filename:' + parameter_filename + ' \n'
    elif 'Restart-data prefix:' in line:
        output_filename = FindOutputFilename(working_path, step)
        line = 'Restart-data prefix:' + output_filename + ' \n'

    return line

# Finds output filename
def FindOutputFilename(working_path, step):
    outputs_path = working_path + "/Outputs"
    CheckIfPathExists(outputs_path)
    if rotate:
        ouput_file_pattern = "/airfoilSol.pval.deform_i="
    else:
        ouput_file_pattern = "/airfoilSol.pval.unsteady_i="
    #CheckIfPathExists(FindFilename(outputs_path, ouput_file_pattern, step + 1))
    return FindFilename(outputs_path, ouput_file_pattern, step + 1)

# Looks for a file matching the given pattern within the given path
def FindFilename(path, pattern, step):
    files_list = glob.glob(path + "/*")
    for file in files_list:
        if file.startswith('%s' % path + '%s' % pattern + '%s' % step + '_' ):
            print('step = ', step)
            print('file1 = ', file)
            filename = file
            break

    if 'domain' in filename:
        position = filename.find('.domain_')
        filename = filename[0:position]
        print('file2 = ', filename)
    return filename
        
    #raise Exception('File: "{}" not found'.format(path + pattern + str(step)))

# Check if file exist and remove it, otherwise print a warning
def RemoveFileIfExists(path):
    if os.path.exists(path):
        os.remove(path)
    else:
        msg = 'The file ' + path + ' does not exist.'
        warnings.warn(msg)

# Check if the path exists
def CheckIfPathExists(path):
    if not os.path.exists(path):
        raise Exception('Path: "{}" not found'.format(path))


###############
# Loop to write all the outputs
###############

#command = "source /work/piquee/Softwares/setPaths/setTau2019.sh"
#subprocess.call(command, shell=True)
#del command

for iteration in range(5,n_step):
    print(iteration)
    # Write Tautoplt.cntl file
    tautoplt_filename = WriteTautoplt(working_path, iteration, para_path_mod, start_step, primary_grid_filename)
    
    # Execute tau2plt to convert output file into dat
    print("Start gathering Solution Data at time %s" % (str(time)))
    command = tau_path + '/gather ' + tautoplt_filename
    subprocess.call(command, shell=True)
    print(command)
    print("Start gathering Solution Data at time %s" % (str(time)))
    del command
    command = tau_path + '/tau2plt ' + tautoplt_filename
    print("Start Writting Solution Data at time %s" % (str(time)))
    subprocess.call(command, shell=True)
    print("Stop Writting Solution Data at time %s" % (str(time)))