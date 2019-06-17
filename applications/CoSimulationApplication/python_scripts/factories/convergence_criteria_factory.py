from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

def CreateConvergenceCriteria(convergence_criteria_settings, solver_wrapper):
    """This function creates and returns the Convergence Criteria used for CoSimulation
    """
    module_name = convergence_criteria_settings["type"].GetString()
    module_full = "KratosMultiphysics.CoSimulationApplication.convergence_criteria."+module_name

    criteria_module = __import__(module_full,fromlist=[module_name])
    return criteria_module.Create(convergence_criteria_settings, solver_wrapper)
