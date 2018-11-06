# co simulation imports
from CoSimulationApplication import *
from base_co_simulation_classes.co_simulation_base_coupled_solver import CoSimulationBaseCoupledSolver
import co_simulation_tools as tools

# Other imports
import os

def Create(custom_settings):
    return GaussSeidelIterativeStrongCouplingSolver(custom_settings)

class GaussSeidelIterativeStrongCouplingSolver(CoSimulationBaseCoupledSolver):
    def __init__(self, custom_settings):
        default_settings = {}
        super(GaussSeidelIterativeStrongCouplingSolver, self).__init__(custom_settings)
        default_settings["convergence_accelerators"] = list    #MANDATORY
        default_settings["convergence_criteria_settings"] = dict    #MANDATORY
        self.settings = tools.ValidateAndAssignInputParameters(default_settings, self.settings, False)
        self.number_of_participants = len( self.settings['participants'] )
        self.max_num_coupling_iterations = self.settings['max_coupling_iterations']

        if not self.number_of_participants == 2:
            raise Exception(tools.bcolors.FAIL + "Exactly two solvers have to be specified for the " + self.__class__.__name__ + "!")

        ### Importing the Participant modules
        self.participants_setting_dict = self.full_settings["coupled_solver_settings"]["participants"]
        self.participating_solver_names = []

        for p in range(0,self.number_of_participants) :
            self.participating_solver_names.append(self.participants_setting_dict[p]['name'])

        ### Making the convergence accelerator for this strategy
        self.convergence_accelerators = self._GetConvergenceAccelerators(self.settings["convergence_accelerators"])

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
        for iteration in range(self.max_num_coupling_iterations):
            if self.echo_level > 0:
                print("\t"+ tools.bcolors.HEADER + str(self._Name()) + " : "+
                                    tools.bcolors.MEGENTA + "Coupling iteration: ", tools.bcolors.BOLD + str(iteration+1) +
                                    " / " + tools.bcolors.BLUE + str(self.max_num_coupling_iterations) + tools.bcolors.ENDC)

            for accelerator in self.convergence_accelerators:
                accelerator.InitializeNonLinearIteration()
            #self.convergence_criteria.InitializeNonLinearIteration()

            for solver_name, solver in self.participating_solvers.items():
                self._SynchronizeInputData(solver_name)
                solver.SolveSolutionStep()
                self._SynchronizeOutputData(solver_name)

            for accelerator in self.convergence_accelerators:
                accelerator.FinalizeNonLinearIteration()
            #self.convergence_criteria.FinalizeNonLinearIteration()

            """if self.convergence_criteria.IsConverged():
                if self.echo_level > 0:
                    print(tools.bcolors.GREEN + "### CONVERGENCE WAS ACHIEVED ###" + tools.bcolors.ENDC )
                break
            else:
                self.convergence_accelerator.ComputeUpdate()"""

            if iteration+1 >= self.max_num_coupling_iterations and self.echo_level > 0:
                print("\t"+tools.bcolors.FAIL + "### CONVERGENCE NOT ACHIEVED IN STRONG COUPLING ITERATIONS ###" + tools.bcolors.ENDC)

    def _Name(self):
        return self.settings['name']


    def PrintInfo(self):
        super(GaussSeidelIterativeStrongCouplingSolver, self).PrintInfo()
        for accelerator in self.convergence_accelerators:
            accelerator.PrintInfo()
