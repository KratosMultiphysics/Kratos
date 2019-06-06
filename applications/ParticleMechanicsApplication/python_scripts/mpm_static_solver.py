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
    @classmethod
    def GetDefaultSettings(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "grid_model_import_settings": {
                "input_type"     : "mdpa",
                "input_filename" : "unknown_name_Grid"
            },
            "pressure_dofs"                      : false,
            "element_search_settings"            : {
                "max_number_of_results"          : 1000,
                "searching_tolerance"            : 1.0E-5
            },
        }""")
        this_defaults.AddMissingParameters(super(MpmSolver, cls).GetDefaultSettings())
        return this_defaults
    def _create_solution_scheme(self):
        return KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()