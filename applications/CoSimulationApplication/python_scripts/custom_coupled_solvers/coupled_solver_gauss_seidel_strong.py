from KratosMultiphysics.CoSimulationApplication.base_classes.base_coupled_solver import CoSimulationBaseCoupledSolver
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools


def Create(custom_settings):
    return GaussSeidelIterativeStrongCoupledSolver(custom_settings)


class GaussSeidelIterativeStrongCoupledSolver(CoSimulationBaseCoupledSolver):
    def __init__(self, custom_settings):
        super(GaussSeidelIterativeStrongCoupledSolver, self).__init__(custom_settings)
        if not self.number_of_participants == 2:
            raise Exception(cs_tools.bcolors.FAIL + "Exactly two solvers have to be specified for the " + self.__class__.__name__ + "!")

        # Importing the Participant modules
        self.participants_setting_dict = self.full_settings["coupled_solver_settings"]["participants"]
        self.participating_solver_names = []

        for p in range(0,self.number_of_participants) :
            self.participating_solver_names.append(self.participants_setting_dict[p]['name'])

        # Comment how the settings are specified has to be consistent!
        # Making the convergence accelerator for this strategy
        self.convergence_accelerators_list = self._CreateConvergenceAccelerators(self.settings["convergence_accelerators"])

        # Creating the convergence criterion
        self.convergence_criteria_list = self._CreateConvergenceCriteria(self.settings["convergence_criteria_settings"]["data_list"])

    def Initialize(self):
        super(GaussSeidelIterativeStrongCoupledSolver, self).Initialize()
        for conv_accelerator in self.convergence_accelerators_list:
            conv_accelerator.Initialize()
        for conv_criteria in self.convergence_criteria_list:
            conv_criteria.Initialize()

    def Finalize(self):
        super(GaussSeidelIterativeStrongCoupledSolver, self).Finalize()
        for conv_accelerator in self.convergence_accelerators_list:
            conv_accelerator.Finalize()
        for conv_criteria in self.convergence_criteria_list:
            conv_criteria.Finalize()

    def InitializeSolutionStep(self):
        super(GaussSeidelIterativeStrongCoupledSolver, self).InitializeSolutionStep()
        for accelerator in self.convergence_accelerators_list:
            accelerator.InitializeSolutionStep()
        for conv_criteria in self.convergence_criteria_list:
            conv_criteria.InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        super(GaussSeidelIterativeStrongCoupledSolver, self).FinalizeSolutionStep()
        for accelerator in self.convergence_accelerators_list:
            accelerator.FinalizeSolutionStep()
        for conv_criteria in self.convergence_criteria_list:
            conv_criteria.FinalizeSolutionStep()

    def SolveSolutionStep(self):
        if self.coupling_started:
            for iteration in range(self.num_coupling_iterations):
                if self.echo_level > 0:
                    cs_tools.PrintInfo("\t"+ cs_tools.bcolors.HEADER + str(self._Name()) ,
                                        cs_tools.bcolors.MEGENTA + "Coupling iteration: ", cs_tools.bcolors.BOLD + str(iteration+1) +
                                        " / " + cs_tools.bcolors.BLUE + str(self.num_coupling_iterations) + cs_tools.bcolors.ENDC)

                for accelerator in self.convergence_accelerators_list:
                    accelerator.InitializeNonLinearIteration()
                for conv_criteria in self.convergence_criteria_list:
                    conv_criteria.InitializeNonLinearIteration()

                for solver_name, solver in self.participating_solvers.items():
                    self._SynchronizeInputData(solver_name)
                    cs_tools.PrintInfo("\t"+cs_tools.bcolors.GREEN + cs_tools.bcolors.BOLD + "SolveSolutionStep for Solver", solver_name + cs_tools.bcolors.ENDC)
                    solver.SolveSolutionStep()
                    self._SynchronizeOutputData(solver_name)

                for conv_criteria in self.convergence_criteria_list:
                    conv_criteria.FinalizeNonLinearIteration()

                is_converged = True #Comment I think this would be suitable for list-comprehension
                for conv_criteria in self.convergence_criteria_list:
                    is_converged = is_converged and conv_criteria.IsConverged()
                if is_converged:
                    if self.echo_level > 0:
                        cs_tools.PrintInfo(cs_tools.bcolors.GREEN + "### CONVERGENCE WAS ACHIEVED ###" + cs_tools.bcolors.ENDC )
                    break

                if iteration+1 >= self.num_coupling_iterations and self.echo_level > 0:
                    cs_tools.PrintWarning("\t"+cs_tools.bcolors.FAIL + "### CONVERGENCE NOT ACHIEVED IN STRONG COUPLING ITERATIONS ###" + cs_tools.bcolors.ENDC)

                for accelerator in self.convergence_accelerators_list:
                    accelerator.FinalizeNonLinearIteration()
        else:
            for solver_name, solver in self.participating_solvers.items():
                cs_tools.PrintInfo("\t"+cs_tools.bcolors.GREEN + cs_tools.bcolors.BOLD + "SolveSolutionStep for Solver", solver_name + cs_tools.bcolors.ENDC)
                solver.SolveSolutionStep()

    def _Name(self):
        return self.settings['name'].GetString()


    def PrintInfo(self):
        super(GaussSeidelIterativeStrongCoupledSolver, self).PrintInfo()
        for accelerator in self.convergence_accelerators_list:
            accelerator.PrintInfo()
        for accelerator in self.convergence_accelerators:
            accelerator.PrintInfo()
