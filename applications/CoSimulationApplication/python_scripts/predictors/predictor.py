from KratosMultiphysics.CoSimulationApplication.co_simulation_component import CoSimulationComponent
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return Predictor(parameters)


# Class Predictor: Base class for extrapolation based on the last two (linear), three (quadratic) or four (cubic) time steps,
# assuming constant time step size.
class Predictor(CoSimulationComponent):
    def __init__(self, _unused):
        super().__init__()

        self.updated = False
        self.dataprev = None
        self.order = None

    def Initialize(self, x):
        super().Initialize()

        self.dataprev = [x.GetNumpyArray()]

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

        self.updated = False

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        if not self.updated:
            raise Exception("Not updated")

    def Linear(self, x_in):
        x = x_in.deepcopy()
        if not self.updated:
            if len(self.dataprev) == 1:
                y = self.dataprev[0]
            else:
                y = 2 * self.dataprev[0] - self.dataprev[1]
            x.SetNumpyArray(y)
            return x
        else:
            raise Exception("Already updated")

    def Quadratic(self, x_in):
        x = x_in.deepcopy()
        if not self.updated:
            if len(self.dataprev) < 3:
                raise Exception("Not sufficient information for quadratic extrapolation")
            y = 3.0 * self.dataprev[0] - 3.0 * self.dataprev[1] + 1.0 * self.dataprev[2]
            x.SetNumpyArray(y)
            return x
        else:
            raise Exception("Already updated")

    def Cubic(self, x_in):
        x = x_in.deepcopy()
        if not self.updated:
            if len(self.dataprev) < 4:
                raise Exception("Not sufficient information for cubic extrapolation")
            y = 4.0 * self.dataprev[0] - 6.0 * self.dataprev[1] + 4.0 * self.dataprev[2] - 1.0 * self.dataprev[3]
            x.SetNumpyArray(y)
            return x
        else:
            raise Exception("Already updated")

    def Predict(self, x):
        pass

    def Update(self, x):
        if not self.updated:
            self.dataprev = [x.GetNumpyArray()] + self.dataprev
            if len(self.dataprev) > self.order + 1:
                self.dataprev.pop()
            self.updated = True
        else:
            raise Exception("Already updated")
