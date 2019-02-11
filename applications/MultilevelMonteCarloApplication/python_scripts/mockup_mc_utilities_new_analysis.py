class MonteCarlo(object):
    '''The base class for the MonteCarlo-classes'''
    def __init__(self,custom_settings):
        self.analysis = None

        default_settings = KratosMultiphysics.Parameters()
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)
        self.convergence = False
        self.current_number_levels = 0

        self.QoI = StatisticalVariable(self.current_number_levels) # <--- here define qoi class inside MC class


    '''
    function defining the Kratos specific application analysis stage of the problem
    '''
    def SetAnalysis(self,application_analysis_stage):
        self.analysis = application_analysis_stage

    '''
    function getting the Kratos specific application analysis stage of the problem previously defined
    '''
    def GetAnalysis(self):
        if (self.analysis is not None):
            return self.analysis
        else:
            print("Please provide a Kratos specific application analysis stage for the current problem the avaiables are: None")



    def CheckConvergence(self,level):


    def AddResults(self,simulation_results):
        '''simulation_results = [MultilevelMonteCarloResults class, level (integer type, not compss.future)]'''
        level = 0
        QoI_value = AddResultsAux_Task(simulation_results,level)
        '''update values of QoI'''
        self.QoI.values[level] = np.append(self.QoI.values[level],QoI_value)

    def InitializeMCPhase(self):

    def FinalizeMCPhase(self):

    def ScreeningInfoInitializeMCPhase(self):

    def ScreeningInfoFinalizeMCPhase(self):

    def SetConvergenceCriteria(self,convergence_string_name):
