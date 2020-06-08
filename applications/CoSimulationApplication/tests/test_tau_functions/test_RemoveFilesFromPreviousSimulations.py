import sys, os

tau_functions_path = os.path.join(os.path.dirname(__file__), '../../python_scripts/helpers')
sys.path.append(tau_functions_path)

import tau_functions as TauFunctions

# Create dummy outputs and mesh directory
TauFunctions.RemoveFilesFromPreviousSimulations()
path = os.getcwd() + '/'
if os.path.exists(path + "Outputs"):
    os.rmdir(path + "Outputs")
os.mkdir('Outputs')
os.mkdir('Mesh')

# Create dummy empty output files
for i in range(3):
    output_file_name = 'Outputs/test_file_' + str(i)
    open(output_file_name,'w').close()

    mesh_file_name = 'Mesh/airfoil_Structured_scaliert.grid.def' + str(i)
    open(mesh_file_name,'w').close()

TauFunctions.RemoveFilesFromPreviousSimulations()

# Check if files have been successfully removed
os.rmdir('Outputs')
os.rmdir('Mesh')