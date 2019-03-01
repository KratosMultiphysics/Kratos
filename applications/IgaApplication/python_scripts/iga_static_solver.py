from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("IgaApplication")

# Import applications
import KratosMultiphysics.IgaApplication as IgaApplication

# Import base class file
import iga_solver


def CreateSolver(model, custom_settings):
    return StaticIgaSolver(model, custom_settings)


class StaticIgaSolver(iga_solver.IgaSolver):
    """The iga static solver.

    This class creates the iga solvers for static analysis. It currently
    supports line search, linear, strategies.

    Public member variables:

    See iga_solver.py for more information.
    """
    def __init__(self, model, custom_settings):
        # Set defaults and validate custom settings in the base class.
        # Construct the base solver.
        super(StaticIgaSolver, self).__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[StaticIgaSolver]:: ", "Construction finished")

    def Initialize(self):
        KratosMultiphysics.Logger.PrintInfo("::[StaticIgaSolver]:: ", "Initializing ...")
        super(StaticIgaSolver, self).Initialize() # The mechanical solver is created here.
        KratosMultiphysics.Logger.PrintInfo("::[StaticIgaSolver]:: ", "Finished initialization.")

    def _create_solution_scheme(self):
        return KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()