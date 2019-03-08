from __future__ import print_function, absolute_import, division  # makes backward compatible with python 2.6 and 2.7


def CreateCoupledSolver(settings):
    """
    This function creates and returns the coupled solver.
    """
    coupled_solver_type = settings["coupled_solver_settings"]["solver_type"]
    module_name = coupled_solver_type.GetString()
    module_full = 'KratosMultiphysics.CoSimulationApplication.custom_coupled_solvers.'+module_name+'_coupled_solver'
    coupled_solver_module = __import__(module_full, fromlist=[module_name])
    return coupled_solver_module.Create(settings)
