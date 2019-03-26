from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_base_coupled_solver import CoSimulationBaseCoupledSolver
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools


def Create(custom_settings):
    return GaussSeidelLooseCoupledSolver(custom_settings)


class GaussSeidelLooseCoupledSolver(CoSimulationBaseCoupledSolver):
    def __init__(self, custom_settings):
        super(GaussSeidelLooseCoupledSolver, self).__init__(custom_settings)
        if not self.number_of_participants == 2:
            raise Exception(cs_tools.bcolors.FAIL + "Exactly two solvers have to be specified for the " + self.__class__.__name__ + "!")

        # Importing the Participant modules
        self.participants_setting_dict = self.full_settings["coupled_solver_settings"]["participants"]
        self.participating_solver_names = []

        for i, participant_settings in enumerate(self.participants_setting_dict):
            self.participating_solver_names.append(participant_settings['name'])

    def Initialize(self):
        super(GaussSeidelLooseCoupledSolver, self).Initialize()

    def Finalize(self):
        super(GaussSeidelLooseCoupledSolver, self).Finalize()

    def InitializeSolutionStep(self):
        super(GaussSeidelLooseCoupledSolver, self).InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        super(GaussSeidelLooseCoupledSolver, self).FinalizeSolutionStep()

    def SolveSolutionStep(self):
        if self.coupling_started:
            for solver_name, solver in self.participating_solvers.items():
                self._SynchronizeInputData(solver_name)
                cs_tools.PrintInfo("\t"+cs_tools.bcolors.GREEN + cs_tools.bcolors.BOLD + "SolveSolutionStep for Solver", solver_name + cs_tools.bcolors.ENDC)
                solver.SolveSolutionStep()
                self._SynchronizeOutputData(solver_name)

        else:
            for solver_name, solver in self.participating_solvers.items():
                cs_tools.PrintInfo("\t"+cs_tools.bcolors.GREEN + cs_tools.bcolors.BOLD + "SolveSolutionStep for Solver", solver_name + cs_tools.bcolors.ENDC)
                solver.SolveSolutionStep()

    def _Name(self):
        return self.settings['name'].GetString()

    def PrintInfo(self):
        super(GaussSeidelLooseCoupledSolver, self).PrintInfo()
