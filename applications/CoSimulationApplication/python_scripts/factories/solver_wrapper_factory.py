from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

from KratosMultiphysics.CoSimulationApplication.factories import base_factory

def CreateSolverWrapper(settings, models, solver_name):
    """This function creates and returns the Wrapper for the Solver used for CoSimulation"""
    return base_factory.Create(settings, [models, solver_name], "KratosMultiphysics.CoSimulationApplication")
