from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Other imports
import co_simulation_tools as cs_tools
data_structure = cs_tools.cs_data_structure
import math

def Create(settings, solver):
    return CoSimulationConvergenceCriteria(settings, solver)

class CoSimulationConvergenceCriteria(object):
    ## __init__ : Initializer
    #
    #  @param self            The object pointer.
    #  @param settings: setting of the convergence criteria
    #  @param solver: solver to which convergence criteria is to be attached
    def __init__(self, settings, solver):
        default_settings = self._GetDefaultSettings()
        self.settings = settings
        self.settings.ValidateAndAssignDefaults(default_settings)
        self.solver = solver
        self.echo_level = 3
        self.data_name = self.settings["data_name"].GetString()
        self.abs_tolerance = self.settings["abs_tolerance"].GetDouble()
        self.rel_tolerance = self.settings["rel_tolerance"].GetDouble()
        self.iteration = 0
        self.initial_residual_norm = 0.0
        self.data_variable = data_structure.KratosGlobals.GetVariable(self.data_name)

        self.data_prev_iter = []
        self.data_current_iter = []

    def _GetDefaultSettings(self):
        default_setting = data_structure.Parameters("""
            {
                "solver"        : "",
                "data_name"     : "",
                "abs_tolerance" : 1e-5,
                "rel_tolerance" : 1e-5
            }
        """)
        return default_setting

    ## Initialize : Initialize function of the class
    #                   To be called once at the beginning of the simulation.
    #
    #  @param self            The object pointer.
    def Initialize(self):
        pass

    ## Finalize : Finalize function of the class
    #                   To be called once at the end of the simulation.
    #
    #  @param self            The object pointer.
    def Finalize(self):
        pass

    ## InitializeSolutionStep : InitializeSolutionStep function of the class.
    #                           To be called once at the beginning of the SolutionStep.
    #
    #  @param self            The object pointer.
    def InitializeSolutionStep(self):
        self.iteration = 0

    ## FinalizeSolutionStep : FinalizeSolutionStep function of the class.
    #                           To be called once at the end of the SolutionStep.
    #
    #  @param self            The object pointer.
    def FinalizeSolutionStep(self):
        pass

    ## InitializeNonLinearIteration : InitializeNonLinearIteration function of the class.
    #                           To be called once at the beginning of the non-linear coupling iteration.
    #
    #  @param self            The object pointer.
    def InitializeNonLinearIteration(self):
        # storing the data for residual calculation
        self.data_prev_iter = cs_tools.GetDataAsList(self.solver, self.data_name)

    ## FinalizeNonLinearIteration : FinalizeNonLinearIteration function of the class.
    #                           To be called once at the end of the non-linear coupling iteration.
    #
    #  @param self            The object pointer.
    def FinalizeNonLinearIteration(self):
        self.data_current_iter = cs_tools.GetDataAsList(self.solver, self.data_name)
        self.iteration = self.iteration + 1


    ## IsConverged : IsConverged function of the class.
    #                  To be called called when the convergence is to be enquired for this criteria
    #
    #  @param self            The object pointer.
    def IsConverged(self):
        residual_list = self._CalculateResidual()
        abs_residual_norm = cs_tools.CalculateNorm(residual_list)
        if(abs_residual_norm == 0):
            abs_residual_norm = 1.0
        if(self.iteration == 1):
            self.initial_residual_norm = abs_residual_norm

        rel_residual_norm = abs_residual_norm / self.initial_residual_norm

        is_converged = abs_residual_norm < self.abs_tolerance or rel_residual_norm < self.rel_tolerance
        if self.echo_level > 1:
            info_msg  = cs_tools.bcolors.HEADER + '\tConvergence for "'+cs_tools.bcolors.BOLD + self.data_name+'": '
            if is_converged:
                info_msg += cs_tools.bcolors.GREEN+ "ACHIEVED"
            else:
                info_msg += cs_tools.bcolors.FAIL + "NOT ACHIEVED"
            info_msg += cs_tools.bcolors.ENDC
            print(info_msg)
        if self.echo_level > 2:
            info_msg  = cs_tools.bcolors.HEADER + "\tabs_norm" + " = " + cs_tools.bcolors.BOLD + str(abs_residual_norm) + " | "
            info_msg += "abs_tol" + " = " + cs_tools.bcolors.BOLD + str(self.abs_tolerance)
            info_msg += " || " + "rel_norm" + " = " + cs_tools.bcolors.BOLD + str(rel_residual_norm) + " | "
            info_msg += "rel_tol = " + cs_tools.bcolors.BOLD + str(self.rel_tolerance)
            print(info_msg)

        return is_converged

    ## PrintInfo : PrintInfo function of the class.
    #                  Prints the information about the object
    #
    #  @param self            The object pointer.
    def PrintInfo(self):
        print("Convergence criteria for data : ", self.data_name, " from solver : ", self.solver.name)

    ## Check : Check function of the class.
    #
    #  @param self            The object pointer.
    def Check(self):
        print("Check from Base Convergence Criteria")


    ## SetEchoLevel : Function to set the echo level of this class.
    #                   Default is 3
    #
    #  @param self            The object pointer.
    def SetEchoLevel(self, level):
        self.echo_level = level


    ## _Name : Gets the name of the class
    #
    #  @param self            The object pointer.
    def _Name(self):
        return self.__class__.__name__

    ## _CalculateResidual : Calculates residual of the data specified in the settings
    #                       Numpy can be used in the variants of this class.
    #                       residual = data_in_current_iter - data_in_previous_iter
    #
    #  @param self            The object pointer
    def _CalculateResidual(self):
        residual = []
        for i in range(len(self.data_prev_iter)):
            residual.append(self.data_current_iter[i] - self.data_prev_iter[i])

        return residual
