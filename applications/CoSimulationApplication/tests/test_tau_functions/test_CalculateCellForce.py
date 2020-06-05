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

# Define nodal pressure
nodal_pressures = [1.4, 3.6, 5.7, 7.9]

# Compute the normal vector
cell_force = TauFunctions.CalculateCellForce(node_ids, nodal_pressures, X, Y, Z)

# Define reference normal
reference_cell_force = np.array([-4.65, 0.0, 4.65])

# Check
np.testing.assert_almost_equal(cell_force, reference_cell_force, decimal=16)