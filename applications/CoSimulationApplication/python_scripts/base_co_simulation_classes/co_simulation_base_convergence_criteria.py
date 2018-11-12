from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Other imports
import co_simulation_tools as cs_tools
data_structure = cs_tools.cs_data_structure

def Create(settings, solver):
    return CoSimulationConvergenceCriteria(settings, solver)

class CoSimulationConvergenceCriteria(object):
    def __init__(self, settings, solver):
        self.settings = settings
        self.solver = solver
        self.echo_level = 0
        self.data_name = self.settings["data_name"].GetString()
        self.abs_tolerance = self.settings["abs_tolerance"].GetDouble()
        self.rel_tolerance = self.settings["rel_tolerance"].GetDouble()
        self.iteration = 0
        self.initial_residual_norm = 0.0
        self.data_variable = KratosMultiphysics.KratosGlobals.GetVariable(self.data_name)

        self.data_prev_iter = []
        self.data_current_iter = []

    def Initialize(self):
        pass

    def Finalize(self):
        pass

    def InitializeSolutionStep(self):
        self.iteration = 0

    def FinalizeSolutionStep(self):
        pass

    def InitializeNonLinearIteration(self):
        # storing the data for residual calculation
        self.data_prev_iter = self._GetDataAsList(self.solver, self.data_name)

    def FinalizeNonLinearIteration(self):
        self.data_current_iter = self._GetDataAsList(self.solver, self.data_name)
        self.iteration = self.iteration + 1

    def IsConverged(self):
        residual_list = self._CalculateResidual()
        abs_residual_norm = self._CalculateNorm(residual_list)
        if(self.iteration == 0):
            self.initial_residual_norm = abs_residual_norm

        rel_residual_norm = abs_residual_norm / self.initial_residual_norm

        if self.echo_level > 1:
            info_msg  = cs_tools.bcolors.HEADER + 'Convergence for "'+cs_tools.bcolors.BOLD + self.data_name+'": '
            if convergence_list[i]:
                info_msg += cs_tools.bcolors.GREEN+ "ACHIEVED"+
            else:
                info_msg += cs_tools.bcolors.FAIL + "NOT ACHIEVED"
            info_msg += cs_tools.bcolors.ENDC
            print(info_msg)
        if self.echo_level > 2:
            info_msg  = cs_tools.bcolors.HEADER + "abs_norm" + " = " + cs_tools.bcolors.BOLD + str(abs_residual_norm) + " | "
            info_msg += "abs_tol" + " = " + cs_tools.bcolors.BOLD + str(self.abs_tolerance)
            info_msg += " || " + "rel_norm" + " = " + cs_tools.bcolors.BOLD + str(rel_residual_norm) + " | "
            info_msg += "rel_tol = " + cs_tools.bcolors.BOLD + str(self.rel_tolerance)
            print(info_msg)

        return abs_residual_norm < self.abs_tolerance or rel_residual_norm < self.rel_tolerance

    def PrintInfo(self):
        classprint(self.lvl, "Convergence Criteria", bold(self._Name()))

    def Check(self):
        print("Check from Base Convergence Criteria")

    def SetEchoLevel(self, level):
        self.echo_level = level

    def _Name(self):
        return self.__class__.__name__

    def _GetDataAsList(self):
        data = []
        data_def = self.solver.data_list[self.data_name]
        data_mesh = self.solver.model[data_def["geometry_name"].GetString()]
        for node in data_mesh.Nodes:
            data.append(node.GetSolutionStepValue(self.data_variable,0))

        return data

    def _CalculateResidual(self):
        residual = []
        for i in range(len(self.data_prev_iter)):
            residual.append(self.data_current_iter[i] - self.data_prev_iter[i])

        return residual

    def _CalculateNorm(self, residual_list):
        norm = 0.0
        for residual in residual_list
            for component in
            norm += 0