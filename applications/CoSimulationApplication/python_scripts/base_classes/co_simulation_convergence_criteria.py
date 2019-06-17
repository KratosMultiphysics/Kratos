from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Other imports
import numpy as np
from numpy import linalg as la
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import classprint, bold, green, red
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

def CreateConvergenceCriteria(settings, solver_wrapper):
    cs_tools.SettingsTypeCheck(settings)
    return CoSimulationConvergenceCriteria(settings, solver_wrapper)

class CoSimulationConvergenceCriteria(object):
    def __init__(self, settings, solver_wrapper):
        self.settings = settings
        self.settings.RecursivelyValidateAndAssignDefaults(self._GetDefaultSettings())

        self.interface_data = solver_wrapper.GetInterfaceData(self.settings["data_name"].GetString())

        self.echo_level = self.settings["echo_level"].GetInt()
        self.abs_tolerance = self.settings["abs_tolerance"].GetDouble()
        self.rel_tolerance = self.settings["rel_tolerance"].GetDouble()


    def Initialize(self):
        pass

    def Finalize(self):
        pass

    def InitializeSolutionStep(self):
        pass

    def FinalizeSolutionStep(self):
        pass

    def InitializeCouplingIteration(self):
        # Saving the previous data (at beginning of iteration) for the computation of the residual
        self.old_data = self.interface_data.GetNumpyArray()

    def FinalizeCouplingIteration(self):
        pass

    def IsConverged(self):
        new_data = self.interface_data.GetNumpyArray()

        residual = new_data - self.old_data
        res_norm = la.norm(residual)
        norm_new_data = la.norm(new_data)
        if norm_new_data < 1e-15:
            norm_new_data = 1.0 # to avoid division by zero
        abs_norm = res_norm / np.sqrt(residual.size)
        rel_norm = res_norm / norm_new_data
        is_converged = abs_norm < self.abs_tolerance or rel_norm < self.rel_tolerance
        if self.echo_level > 1:
            info_msg  = 'Convergence for "'+bold(self.interface_data.variable.Name())+'": '
            if is_converged:
                info_msg += green("ACHIEVED")
            else:
                info_msg += red("NOT ACHIEVED")
            classprint(self._Name(), info_msg)
        if self.echo_level > 2:
            info_msg  = bold("abs_norm") + " = " + str(abs_norm) + " | "
            info_msg += bold("abs_tol")  + " = " + str(self.abs_tolerance) + " || "
            info_msg += bold("rel_norm") + " = " + str(rel_norm) + " | "
            info_msg += bold("rel_tol")  + " = " + str(self.rel_tolerance)
            classprint(self._Name(), info_msg)

        return is_converged

    def PrintInfo(self):
        classprint("Convergence Criteria", bold(self._Name()))

    def Check(self):
        print("ConvCrit does not implement Check yet!")

    def SetEchoLevel(self, level):
        self.echo_level = level

    def _Name(self):
        return self.__class__.__name__

    @classmethod
    def _GetDefaultSettings(cls):
        return cs_tools.cs_data_structure.Parameters("""{
            "type"       : "UNSPECIFIED",
            "solver"     : "UNSPECIFIED",
            "data_name"  : "UNSPECIFIED",
            "abs_tolerance" : 1e-5,
            "rel_tolerance" : 1e-5,
            "echo_level" : 0
        }""")
