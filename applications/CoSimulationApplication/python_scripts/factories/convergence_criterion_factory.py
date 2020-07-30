from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

from KratosMultiphysics.CoSimulationApplication.factories import base_factory

def CreateConvergenceCriterion(convergence_criterion_settings):
    """This function creates and returns the Convergence Criterion used for CoSimulation"""
    return base_factory.Create(convergence_criterion_settings, [], "KratosMultiphysics.CoSimulationApplication.convergence_criteria")
