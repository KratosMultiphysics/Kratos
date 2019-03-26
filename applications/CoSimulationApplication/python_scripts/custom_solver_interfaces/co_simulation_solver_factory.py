from __future__ import print_function, absolute_import, division  # makes backward compatible with python 2.6 and 2.7


def CreateSolverInterface(solver_name, settings):
    """
    This function creates and returns the solver interface."
    """
    solver_type = settings["solver_type"]
    module_name = solver_type.GetString()
    module_full = 'KratosMultiphysics.CoSimulationApplication.custom_solver_interfaces.'+module_name
    solver_module = __import__(module_full, fromlist=[module_name])
    return solver_module.Create(solver_name, settings)
