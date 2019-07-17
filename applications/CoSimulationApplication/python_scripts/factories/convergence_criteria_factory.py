from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

from . import base_factory

def CreateConvergenceCriteria(convergence_criteria_settings, solver_wrapper):
    """This function creates and returns the Convergence Criteria used for CoSimulation"""
    return base_factory.Create(convergence_criteria_settings, [solver_wrapper], "KratosMultiphysics.CoSimulationApplication.convergence_criteria")
