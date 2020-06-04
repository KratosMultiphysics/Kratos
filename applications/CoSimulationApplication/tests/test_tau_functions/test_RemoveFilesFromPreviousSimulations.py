import sys, os

parent_path = os.path.join(os.path.dirname(__file__), '../../python_scripts/helpers')
sys.path.append(parent_path)

import tau_functions as TauFunctions

TauFunctions.RemoveFilesFromPreviousSimulations()