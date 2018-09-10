# co simulation imports
import tools

# Other imports
import os

def Create(custom_settings):
    return GaussSeidelIterativeStrongCouplingSolver(custom_settings)

class GaussSeidelIterativeStrongCouplingSolver(CoSimApp.CoSimulationBaseCouplingStrategy):
    def __init__(self, custom_settings):

        ##settings string in json format
        # default values for all available settings
        # for mandatory settings, the type is defined
        defaultSettings = {}
        defaultSettings["echo_level"] = 1
        defaultSettings["convergence_acceleration"] = dict    #MANDATORY
        defaultSettings["convergence_criteria"] = dict    #MANDATORY
        defaultSettings["max_iteration_per_step"] = 10
        defaultSettings["participants"] = list
        self.settings = ValidateAndAssignInputParameters(defaultSettings, custom_settings)

        self.number_of_participants = ( self.settings['participants'] ).size()
        print('Number of participants :: ', self.number_of_participants)
        if(self.number_of_participants <= 1):
            raise RuntimeError('Number of participants in a Co-Simulation strategy cannot be less than 2. \nPlease check the input ... !')
        if(self.number_of_participants > 2):
            raise RuntimeError('Number of participants in a Co-Simulation strategy cannot be more than 2. \nPlease check the input ... !')

        ### Importing the Participant modules
        self.participants = self.settings["co_simulation_solver_settings"]["participants"]
        self.participating_solver_names = []

        for p in range(0,self.number_of_participants) :
            self.participating_solver_names.append(self.participants[p]['name'].GetString())

        solvers = tools.GetSolvers(self.settings['solvers'])

        ### obtaining the participants
        self.participating_solvers = tools.Extract(solvers, self.participating_solver_names)

        ### obtaining participants co sim details
        self.solver_cosim_details = tools.GetSolverCoSimulationDetails(self.participants)

        ### Making the convergence accelerator for this strategy
        self.convAccelerator = None

        ### Creating the convergence criterion
        #self.conv_criterion = CoSimApp.CoSimulationBaseConvergenceCriterion(self.settings['residual_relative_tolerance'].GetDouble(), self.settings['residual_relative_tolerance'].GetDouble())

    def Initialize(self):
        for solver_name, solver in self.participating_solvers.items():
            print('Initializing solver :: ', solver_name)
            solver.Initialize()

    def Finalize(self):
        for solver_name, solver in self.participating_solvers.items():
            print('Finalizing solver :: ', solver_name)
            solver.Finalize()


    def InitializeTimeStep(self):
        for solver_name, solver in self.participating_solvers.items():
            print('InitializeSolutionStep for solver :: ', solver_name)
            solver.InitializeTimeStep()

    def FinalizeTimeStep(self):
        for solver_name, solver in self.participating_solvers.items():
            print('FinalizeTimeStep for solver :: ', solver_name)
            solver.FinalizeTimeStep()

    def ImportData(self, DataName, FromClient):
        pass
    def ImportMesh(self, MeshName, FromClient):
        pass

    def ExportData(self, DataName, ToClient):
        pass
    def ExportMesh(self, MeshName, ToClient):
        pass

    def MakeDataAvailable(self, DataName, ToClient):
        pass
    def MakeMeshAvailable(self, MeshName, ToClient):
        pass


    def SolveTimeStep(self):
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


