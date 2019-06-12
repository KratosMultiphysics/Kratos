from __future__ import print_function, absolute_import, division  # makes backward compatible with python 2.6 and 2.7

def CreateSolverInterface(model, settings, solver_name):
    """
    This function creates and returns the coupled solver.
    """
    co_sim_solver_type = settings["type"].GetString()
    module_full  = 'KratosMultiphysics.CoSimulationApplication.'+co_sim_solver_type
    co_sim_solver_module = __import__(module_full, fromlist=[co_sim_solver_type])
    return co_sim_solver_module.CreateSolver(model, settings, solver_name)
