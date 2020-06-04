import sys, os
import numpy as np

tau_functions_path = os.path.join(os.path.dirname(__file__), '../../python_scripts/helpers')
sys.path.append(tau_functions_path)

import tau_functions as TauFunctions

# Define cell geometry
X = [0, 1, 1, 0]
Y = [0, 0, 1, 1]
Z = [0, 1, 1, 0]

node_ids = np.array([0, 1, 2, 3], dtype=int)

# Compute the normal vector
cell_area = TauFunctions.CalculateCellArea(X, Y, Z, node_ids)

# Define reference normal
reference_area = 1.4142135623730951

# Check
np.testing.assert_almost_equal(cell_area, reference_area, decimal=16)