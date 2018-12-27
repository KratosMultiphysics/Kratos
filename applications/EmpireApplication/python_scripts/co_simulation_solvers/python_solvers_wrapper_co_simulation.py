from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the base class of the cosim-solvers
from co_simulation_solvers.co_simulation_base_solver import CoSimulationBaseSolver

def CreateSolver(cosim_solver_settings, level):
    """This function creates and returns the solvers used for CoSimulation
    The solver-module has to be on the PYTHONPATH
    Naming-Convention: The module-file has to end with "_solver"
    """
    if (type(cosim_solver_settings) != dict):
        raise Exception("Input is expected to be provided as a python dictionary")

    solver_module_name = cosim_solver_settings["solver_type"] + "_solver"

    solver_module = __import__(solver_module_name)
    solver = solver_module.CreateSolver(cosim_solver_settings, level+1)
    if not isinstance(solver, CoSimulationBaseSolver):
        err_msg  = 'The requested solver "' + solver_type
        err_msg += '" does not derive from "CoSimulationBaseSolver"!'
        raise Exception(err_msg)
    return solver