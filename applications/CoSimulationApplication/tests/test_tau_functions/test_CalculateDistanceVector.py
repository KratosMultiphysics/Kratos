import sys, os
import numpy as np

tau_functions_path = os.path.join(os.path.dirname(__file__), '../../python_scripts/helpers')
sys.path.append(tau_functions_path)

import tau_functions as TauFunctions

# Populate two nodes with random coordinates
X = [1.6, 19.1]
Y = [-0.3, 3.6]
Z = [4.36, -6.78]

# Compute the distance vector
distance = TauFunctions.CalculateDistanceVector(X, Y, Z, 0, 1)

# Define reference distance
reference_distance = np.array([17.5, 3.9, -11.14])

# Check
np.testing.assert_almost_equal(reference_distance, distance, decimal=16)