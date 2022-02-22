# Importing Kratos
import KratosMultiphysics

# Other imports
from importlib import import_module

def CreateSolverByParameters(model, solver_settings, parallelism):


    solver_type = solver_settings["solver_type"].GetString()

    try_import_custom_solver = False

    if solver_type == "neural_network" or solver_type == "Neural_network":
        solver_module_name = "neural_network_solver"
    elif solver_type == "pinn" or solver_type == "PINN":
        solver_module_name = "physics_informed_neural_network_solver"
    elif solver_type == "neural_network_mesh_moving" or solver_type == "Neural_network_mesh_moving":
        solver_module_name = "neural_network_mesh_moving_solver"
    else:
        print("the selected solver is not avaibale in the solvers, available solvers are neural_network, neural_network_mesh_moving and physics_informed_neural_network")

    kratos_module = "KratosMultiphysics.NeuralNetworkApplication"

    solver = import_module(kratos_module + "." + solver_module_name).CreateSolver(solver_settings, model)

    return solver


def CreateSolver(model, custom_settings):

    if not isinstance(model, KratosMultiphysics.Model):
        raise Exception("input is expected to be provided as a Kratos Model object")#

    if not isinstance(custom_settings, KratosMultiphysics.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    solver_settings = custom_settings
    parallelism = custom_settings["parallel_type"].GetString()

    return CreateSolverByParameters(model, solver_settings, parallelism)
