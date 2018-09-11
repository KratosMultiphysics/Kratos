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

        print(self.participating_solver_names)
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
        max_iter = self.settings["co_simulation_solver_settings"]['max_iteration_per_step'].GetInt()
        iter = 0
        while(iter < max_iter):
        #while(iter < maxIter and not self.conv_criterion.IsConverged(self.appOne.GetModelPart(), DISPLACEMENT)):
            print('\t############## ')
            print('\tCoupling iteration :: ', iter)
            for i in range (0,len(self.participating_solver_names)): #prevent systems to be solved in random order
                solver_name = self.participating_solver_names[i]
                solver = self.participating_solvers[solver_name]
                print('\tSolving :: ', solver_name)
                self.__SynchronizeInputData(solver)
                solver.SolveTimeStep()
                self.__SynchronizeOutputData(solver)

            iter = iter + 1

    ###############################
    def __SynchronizeInputData(self, solver):
        input_data_list = self.solver_cosim_details[solver.name]['input_data_list'] #.Name()
        for i in range(0, input_data_list.size()):
            from_solver  = self.participating_solvers[ input_data_list[i]['from_solver'].GetString() ]
            self.participating_solvers[solver.name].ImportData(input_data_list[i]['data_name'].GetString(), from_solver)

    def __SynchronizeOutputData(self, solver):
        output_data_list = self.solver_cosim_details[solver.name]['output_data_list']   #.Name()
        for i in range(0, output_data_list.size()):
            to_solver  = self.participating_solvers[ output_data_list[i]['to_solver'].GetString() ]
            self.participating_solvers[solver.name].ExportData(output_data_list[i]['data_name'].GetString(), to_solver)


