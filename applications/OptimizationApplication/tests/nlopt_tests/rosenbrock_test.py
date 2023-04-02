import KratosMultiphysics.OptimizationApplication as KratosOA
import numpy as np
#import pybind11

# Define the Rosenbrock function
def rosenbrock(x):
    return (1 - x[0])**2 + 100 * (x[1] - x[0]**2)**2

# Define the gradient of the Rosenbrock function
def rosenbrock_gradient(x):
    return np.array([
        -2 * (1 - x[0]) - 400 * x[0] * (x[1] - x[0]**2),
        200 * (x[1] - x[0]**2)
    ])

# Create an instance of the NLOptOptimizer class
optimizer = KratosOA.NLOptOptimizer(KratosOA.nlopt_algorithm.NLOPT_LD_LBFGS, 2)
optimizer.set_objective_function(rosenbrock)
optimizer.set_gradient(rosenbrock_gradient)
optimizer.set_lower_bounds([-5, -5])
optimizer.set_upper_bounds([5, 5])
optimizer.set_initial_guess([0, 0])
optimizer.set_relative_tolerance(1e-6)
optimizer.set_max_iterations(1000)

# # Optimize the Rosenbrock function
#result = optimizer.optimize()
# print("Optimal solution:", result)
# print("Optimal value:", rosenbrock(result))