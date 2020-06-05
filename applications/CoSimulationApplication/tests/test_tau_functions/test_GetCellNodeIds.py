import sys, os
import numpy as np

tau_functions_path = os.path.join(os.path.dirname(__file__), '../../python_scripts/helpers')
sys.path.append(tau_functions_path)

import tau_functions as TauFunctions

# Define elem_connectivities
elem_connectivities = np.array([24, 35, 87, 94], dtype=int)
cell = 0

# Compute the normal vector
node_ids = TauFunctions.GetCellNodeIds(elem_connectivities, cell)

# Check
np.testing.assert_almost_equal(node_ids, elem_connectivities - 1, decimal=16)