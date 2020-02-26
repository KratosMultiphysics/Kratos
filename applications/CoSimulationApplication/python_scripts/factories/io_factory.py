from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

from KratosMultiphysics.CoSimulationApplication.factories import base_factory

def CreateIO(io_settings, model, solver_name, io_name):
    """This function creates and returns the IO used for CoSimulation"""
    return base_factory.Create(io_settings, [model, solver_name], "KratosMultiphysics.CoSimulationApplication.solver_wrappers", io_name)
