# co simulation imports
from co_simulation_base_coupled_solver import CoSimulationBaseCoupledSolver
import co_simulation_tools as tools

# Other imports
import os

def Create(custom_settings):
    return GaussSeidelIterativeStrongCouplingSolver(custom_settings)

class GaussSeidelIterativeStrongCouplingSolver(CoSimulationBaseCoupledSolver):
    def __init__(self, custom_settings):
        default_settings = {}
        super(GaussSeidelIterativeStrongCouplingSolver, self).__init__(custom_settings)
        default_settings["convergence_accelerator_settings"] = dict    #MANDATORY
        default_settings["convergence_criteria_settings"] = dict    #MANDATORY
        self.settings = tools.ValidateAndAssignInputParameters(default_settings, self.settings, False)
        self.number_of_participants = len( self.settings['participants'] )

        if not self.number_of_participants == 2:
            raise Exception("Exactly two solvers have to be specified for the " + self.__class__.__name__ + "!")

        ### Importing the Participant modules
        self.participants_setting_dict = self.full_settings["coupled_solver_settings"]["participants"]
        self.participating_solver_names = []

        for p in range(0,self.number_of_participants) :
            self.participating_solver_names.append(self.participants_setting_dict[p]['name'])

        solvers = tools.GetSolvers(self.full_settings['solvers'])

        ### Making the convergence accelerator for this strategy
        self.convergence_accelerator = None

        ### Creating the convergence criterion
        #self.convergence_criteria = CoSimApp.CoSimulationBaseConvergenceCriterion(self.settings['residual_relative_tolerance'].GetDouble(), self.settings['residual_relative_tolerance'].GetDouble())

    def Initialize(self):
        super(GaussSeidelIterativeStrongCouplingSolver, self).Initialize()
        #self.convergence_accelerator.Initialize()
        #self.convergence_criteria.Initialize()

    def Finalize(self):
        super(GaussSeidelIterativeStrongCouplingSolver, self).Finalize()
        #self.convergence_accelerator.Finalize()
        #self.convergence_criteria.Finalize()


    def InitializeSolutionStep(self):
        super(GaussSeidelIterativeStrongCouplingSolver, self).InitializeSolutionStep()
        #self.convergence_accelerator.InitializeSolutionStep()
        #self.convergence_criteria.InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        super(GaussSeidelIterativeStrongCouplingSolver, self).InitializeSolutionStep()
        #self.convergence_accelerator.InitializeSolutionStep()
        #self.convergence_criteria.InitializeSolutionStep()

    def SolveSolutionStep(self):
        pass

        """
        for iteration in range(self.num_coupling_iterations):
            if self.echo_level > 0:
                couplingsolverprint(self.lvl, self._Name(),
                                    cyan("Coupling iteration:"), bold(str(k+1)+" / " + str(self.num_coupling_iterations)))

            self.convergence_accelerator.InitializeNonLinearIteration()
            self.convergence_criteria.InitializeNonLinearIteration()

            for solver_name, solver in self.participating_solvers.items():
                self._SynchronizeInputData(solver, solver_name)
                solver.SolveSolutionStep()
                self._SynchronizeOutputData(solver, solver_name)

            self.convergence_accelerator.FinalizeNonLinearIteration()
            self.convergence_criteria.FinalizeNonLinearIteration()

            if self.convergence_criteria.IsConverged():
                if self.echo_level > 0:
                    couplingsolverprint(self.lvl, self._Name(), green("### CONVERGENCE WAS ACHIEVED ###"))
                break
            else:
                self.convergence_accelerator.ComputeUpdate()

            if iteration+1 >= self.num_coupling_iterations and self.echo_level > 0:
                couplingsolverprint(self.lvl, self._Name(), red("XXX CONVERGENCE WAS NOT ACHIEVED XXX"))
        """

