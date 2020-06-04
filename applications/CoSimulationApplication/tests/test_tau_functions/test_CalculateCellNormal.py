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
cell_normal = TauFunctions.CalculateCellNormal(X, Y, Z, node_ids)

# Define reference normal
reference_normal = np.array([-0.7071067811865475, 0.0, 0.7071067811865475])

# Check
np.testing.assert_almost_equal(cell_normal, reference_normal, decimal=16)