from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Other imports
import numpy as np
from numpy import linalg as la
from co_simulation_tools import classprint, bold, green, red
import co_simulation_tools as cs_tools

def CreateConvergenceCriteria(settings, solvers, level):
    return CoSimulationConvergenceCriteria(settings, solvers, level)

class CoSimulationConvergenceCriteria(object):
    def __init__(self, settings, solvers, level):
        self.settings = settings
        self.solvers = solvers
        self.echo_level = 0
        self.lvl = level
        self.abs_tolerances = [data_entry["abs_tolerance"] for data_entry in self.settings["data_list"]]
        self.rel_tolerances = [data_entry["rel_tolerance"] for data_entry in self.settings["data_list"]]

        ## Here we preallocate the arrays that will be used to exchange data
        num_data = len(self.settings["data_list"])
        self.old_data = [np.array([]) for e in range(num_data)]
        self.new_data = np.array([])

    def Initialize(self):
        pass

    def Finalize(self):
        pass

    def InitializeSolutionStep(self):
        pass

    def FinalizeSolutionStep(self):
        pass

    def InitializeNonLinearIteration(self):
        # Saving the previous data (at beginning of iteration) for the computation of the residual
        for i, data_entry in enumerate(self.settings["data_list"]):
            solver = self.solvers[data_entry["solver"]]
            data_name = data_entry["data_name"]
            cs_tools.ImportArrayFromSolver(solver, data_name, self.old_data[i])

    def FinalizeNonLinearIteration(self):
        pass

    def IsConverged(self):
        convergence_list = []
        for i, data_entry in enumerate(self.settings["data_list"]):
            solver = self.solvers[data_entry["solver"]]
            data_name = data_entry["data_name"]
            cs_tools.ImportArrayFromSolver(solver, data_name, self.new_data)

            residual = self.new_data - self.old_data[i]
            res_norm = la.norm(residual)
            norm_new_data = la.norm(self.new_data)
            if norm_new_data < 1e-15:
                norm_new_data = 1.0 # to avoid division by zero
            abs_norm = res_norm / np.sqrt(residual.size)
            rel_norm = res_norm / norm_new_data
            convergence_list.append(abs_norm < self.abs_tolerances[i] or rel_norm < self.rel_tolerances[i])
            if self.echo_level > 1:
                info_msg  = 'Convergence for "'+bold(data_entry["data_name"])+'": '
                if convergence_list[i]:
                    info_msg += green("ACHIEVED")
                else:
                    info_msg += red("NOT ACHIEVED")
                classprint(self.lvl, self._Name(), info_msg)
            if self.echo_level > 2:
                info_msg  = bold("abs_norm")+" = " + str(abs_norm) + " | "
                info_msg += bold("abs_tol")+" = " + str(self.abs_tolerances[i])
                info_msg += " || "+bold("rel_norm")+" = " + str(rel_norm) + " | "
                info_msg += bold("rel_tol") +" = " + str(self.rel_tolerances[i])
                classprint(self.lvl, self._Name(), info_msg)

        return min(convergence_list) # return false if any of them did not converge!

    def PrintInfo(self):
        classprint(self.lvl, "Convergence Criteria", bold(self._Name()))

    def Check(self):
        print("ConvCrit does not implement Check yet!")

    def SetEchoLevel(self, level):
        self.echo_level = level

    def _Name(self):
        return self.__class__.__name__