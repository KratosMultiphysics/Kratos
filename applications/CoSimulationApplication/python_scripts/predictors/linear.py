from KratosMultiphysics.CoSimulationApplication.co_simulation_component import CoSimulationComponent
import co_simulation_tools as cs_tools
data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return PredictorLinear(parameters)


class LinearPredictor(CoSimulationComponent):
    def __init__(self, parameters):
        self.data_prev_iter = []
        self.data_current_iter = []
        self.ResidualStorage = deque( maxlen = 2 ) # 0 is the latest , 1 is the previous




    def InitializeSolutionStep(self):
        self.iteration = 0
        self.alpha_old = self.initial_alpha
        self.data_prev_iter = self.data_current_iter
