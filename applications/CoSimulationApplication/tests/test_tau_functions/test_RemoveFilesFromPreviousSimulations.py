import sys, os

tau_functions_path = os.path.join(os.path.dirname(__file__), '../../python_scripts/helpers')
sys.path.append(tau_functions_path)

import tau_functions as TauFunctions

TauFunctions.RemoveFilesFromPreviousSimulations()