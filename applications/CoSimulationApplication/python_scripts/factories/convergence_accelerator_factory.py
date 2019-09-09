from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

from KratosMultiphysics.CoSimulationApplication.factories import base_factory

def CreateConvergenceAccelerator(convergence_accelerator_settings, solver_wrapper):
    """This function creates and returns the Convergence Accelerator used for CoSimulation"""
    return base_factory.Create(convergence_accelerator_settings, [solver_wrapper], "KratosMultiphysics.CoSimulationApplication.convergence_accelerators")
