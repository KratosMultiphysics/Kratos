import sys, os
import numpy as np

tau_functions_path = os.path.join(os.path.dirname(__file__), '../../python_scripts/helpers')
sys.path.append(tau_functions_path)

import tau_functions as TauFunctions

# Define nodal pressure
nodal_pressures = [1.4, 3.6, 5.7, 7.9]

node_ids = np.array([0, 1, 2, 3], dtype=int)

# Compute the normal vector
cell_pressure = TauFunctions.CalculateCellPressure(nodal_pressures, node_ids)

# Define reference normal
reference_pressure = 4.65

# Check
np.testing.assert_almost_equal(cell_pressure, reference_pressure, decimal=16)