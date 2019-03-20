from __future__ import print_function, absolute_import, division  # makes backward compatible with python 2.6 and 2.7

def CreateCoupledSolver(model, settings, solver_name):
    """
    This function creates and returns the coupled solver.
    """
    coupled_solver_type = settings["solver_type"]
    module_full  = 'KratosMultiphysics.CoSimulationApplication.coupled_solvers.'+coupled_solver_type
    module_full += '_coupled_solver'
    coupled_solver_module = __import__(module_full, fromlist=[coupled_solver_type])
    return coupled_solver_module.Create(model, settings, solver_name)
