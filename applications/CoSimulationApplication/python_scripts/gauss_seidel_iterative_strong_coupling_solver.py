import KratosMultiphysics
# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("CoSimulationApplication")
import KratosMultiphysics.CoSimulationApplication as CoSimApp

# Other imports
import os

def CreateSolver(custom_settings):
    return GaussSeidelIterativeStrongCouplingSolver(custom_settings)

class GaussSeidelIterativeStrongCouplingSolver(CoSimApp.CoSimulationBaseCouplingStrategy):
    def __init__(self, custom_settings):

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "type"                        : "gauss_seidel_iterative_strong_coupling_solver",
            "echo_level"                         : 1,
            "convergence_acceleration"           : {
                "type" : "constant",
                "settings" : {
                    "factor":0.1
                }
            },
            "residual_relative_tolerance"        : 0.001,
            "residual_absolute_tolerance"        : 1e-6,
            "max_iteration_per_step"             : 10,
            "participants"                       : []          
        }""")

        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)
        numberOfParticipants = self.settings['participants'].size()
        print('Number of participants :: ', numberOfParticipants)
        if(numberOfParticipants <= 1):
            raise RuntimeError('Number of participants in a Co-Simulation strategy cannot be less than 2. \nPlease check the input ... !')     
        if(numberOfParticipants > 2):
            raise RuntimeError('Number of participants in a Co-Simulation strategy cannot be more than 2. \nPlease check the input ... !')

        ### Importing the Participant modules
        solverOneModule = __import__(self.settings['participants'][0]['type'].GetString())
        solverTwoModule = __import__(self.settings['participants'][1]['type'].GetString())

        ### Creating the participants
        self.appOne = solverOneModule.CreateSolver(self.settings['participants'][0])
        self.appTwo = solverTwoModule.CreateSolver(self.settings['participants'][1])

        ### Making the convergence accelerator for this strategy
        self.convAccelerator = None

        ### Creating the convergence criterion
        self.conv_criterion = CoSimApp.CoSimulationBaseConvergenceCriterion(self.settings['residual_relative_tolerance'].GetDouble(), self.settings['residual_relative_tolerance'].GetDouble())

        for i in range(self.settings['participants'].size()):            
            print('Participant type :: ', self.settings['participants'][i]['type'].GetString())   

        #super(IterativeStrongCouplingSolver, self).__init__(self.appOne, self.appTwo, convAccelerator)
    def InitializeSolutionStep(self):
        self.appOne.InitializeSolutionStep()
        self.appTwo.InitializeSolutionStep()
        
        self.appOne.SynchronizeInputData()
        self.appTwo.SynchronizeInputData()
    def Predict(self):
        print("GaussSeidelIterativeStrongCouplingSolver predicting ... !!")
    def Initialize(self):
        self.appOne.Initialize()
        self.appTwo.Initialize()
    def Clear(self):
        pass
    def IsConverged(self):
        pass
    def CalculateOutputData(self):
        pass
    def FinalizeSolutionStep(self):
        self.appOne.SynchronizeOutputData()
        self.appTwo.SynchronizeOutputData()
    def SolveSolutionStep(self):
        pass
    def SetEchoLevel(self):
        pass
    def GetEchoLevel(self):
        pass
    def GetResidualNorm(self):
        pass
    def SynchronizeInputData(self):
        pass
    def SynchronizeOutputData(self):
        pass    
    def Solve(self):
        maxIter = self.settings['max_iteration_per_step'].GetInt()
        iter = 0
        while(iter < maxIter):
        #while(iter < maxIter and not self.conv_criterion.IsConverged(self.appOne.GetModelPart(), DISPLACEMENT)):
            print('\t############## ')
            print('\tCoupling iteration :: ', iter)
            self.appOne.SynchronizeInputData()
            self.appOne.Solve()
            self.appOne.SynchronizeOutputData()
            self.appTwo.SynchronizeInputData()
            self.appTwo.Solve()
            self.appTwo.SynchronizeOutputData()
            
            iter = iter + 1
        