from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import base class file
from KratosMultiphysics.ParticleMechanicsApplication.mpm_solver import MpmSolver


def CreateSolver(model, custom_settings):
    return MpmStaticSolver(model, custom_settings)


class MpmStaticSolver(MpmSolver):
    """The mpm solver for static problems.

    See mpm_solver.py for more information.
    """

    def __init__(self, model, custom_settings):
        # Set defaults and validate custom settings in the base class.
        # Construct the base solver.
        super(MpmStaticSolver, self).__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[MpmStaticSolver]:: ", "Construction finished")

    def _create_solution_scheme(self):
        return KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()