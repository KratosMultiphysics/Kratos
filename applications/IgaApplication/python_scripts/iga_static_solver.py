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
        self.print_on_rank_zero("::[StaticIgaSolver]:: ", "Construction finished")

    def Initialize(self):
        self.print_on_rank_zero("::[StaticIgaSolver]:: ", "Initializing ...")
        super(StaticIgaSolver, self).Initialize() # The mechanical solver is created here.
        self.print_on_rank_zero("::[StaticIgaSolver]:: ", "Finished initialization.")

    def Solve(self):
        super(StaticIgaSolver, self).Solve()

    def _create_solution_scheme(self):
        return KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()

    def _create_iga_solution_strategy(self):
        analysis_type = self.settings["analysis_type"].GetString()
        if analysis_type == "linear":
            iga_solution_strategy = self._create_linear_strategy()
        elif analysis_type == "non_linear":
            if(self.settings["line_search"].GetBool() == False):
                iga_solution_strategy = self._create_newton_raphson_strategy()
            else:
                iga_solution_strategy = self._create_line_search_strategy()
        else:
            err_msg =  "The requested analysis type \"" + analysis_type + "\" is not available!\n"
            err_msg += "Available options are: \"linear\", \"non_linear\""
            raise Exception(err_msg)
        return iga_solution_strategy