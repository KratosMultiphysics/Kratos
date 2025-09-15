# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_convergence_criteria import CoSimulationConvergenceCriteria

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
import KratosMultiphysics.CoSimulationApplication.colors as colors

# Other imports
import numpy as np
from numpy import linalg as la

def Create(settings):
    cs_tools.SettingsTypeCheck(settings)
    return RelativeNormInitialResidualConvergenceCriteria(settings)

class RelativeNormInitialResidualConvergenceCriteria(CoSimulationConvergenceCriteria):
    def __init__(self, settings):
        super().__init__(settings)

        self.abs_tolerance = self.settings["abs_tolerance"].GetDouble()
        self.rel_tolerance = self.settings["rel_tolerance"].GetDouble()
        self.label = self.settings["label"].GetString()

    def InitializeSolutionStep(self):
        self.initial_iteration = True

    def IsConverged(self, residual, current_data):
        abs_norm = la.norm(residual) / np.sqrt(residual.size)

        if self.initial_iteration:
            self.initial_iteration = False
            self.initial_norm = abs_norm

        rel_norm = abs_norm / self.initial_norm

        is_converged = bool(abs_norm < self.abs_tolerance or rel_norm < self.rel_tolerance)

        info_msg = ""

        if self.echo_level > 1:
            info_msg  = 'Convergence '

            if self.label != "":
                info_msg += 'for "{}": '.format(self.label)

            if is_converged:
                info_msg += colors.green("ACHIEVED")
            else:
                info_msg += colors.red("NOT ACHIEVED")

        if self.echo_level > 2:
            info_msg += '\n\t abs-norm = {:.2e} | abs-tol = {} || rel-norm = {:.2e} | rel-tol = {}'.format(abs_norm, self.abs_tolerance, rel_norm, self.rel_tolerance)

        if info_msg != "":
            cs_tools.cs_print_info(self._ClassName(), info_msg)

        return is_converged

    def Check(self):
        pass

    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "abs_tolerance" : 1e-5,
            "rel_tolerance" : 1e-5,
            "label"         : ""
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())
        return this_defaults

