import co_simulation_tools as cs_tools
data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return PredictorLinear(parameters)


class LinearPredictor(object):
    def __init__( self, settings, solver ):
        self.data_prev_iter = []
        self.data_current_iter = []
        self.ResidualStorage = deque( maxlen = 2 ) # 0 is the latest , 1 is the previous




    def InitializeSolutionStep(self):
        self.iteration = 0
        self.alpha_old = self.initial_alpha
        self.data_prev_iter = self.data_current_iter

    def FinalizeSolutionStep(self):
        pass

    ## InitializeNonLinearIteration : Function initializes the non linear iteration (coupling iteration)
    #                               Called at the beginning of the nonlinear iteration (coupling iteration)
    #
    def InitializeNonLinearIteration(self):
        pass

    ## FinalizeNonLinearIteration : Function finalizes the non linear iteration (coupling iteration)
    #                               Called at the end of the nonlinear iteration (coupling iteration)
    #
    def FinalizeNonLinearIteration(self):
        pass

    ## PrintInfo : Function to print the information of the convergence accelerator
    #
    def PrintInfo(self):
        cs_tools.PrintInfo(cs_tools.bcolors.HEADER + "This is an object of linear predictor." +cs_tools.bcolors.ENDC )

    ## Check : Function to Check the setup of the convergence accelerator
    #
    def Check(self):
        pass

    ## _Name : Function to get the name of the convergence accelerator
    #
    def _Name(self):
        return "linear_predictor"
