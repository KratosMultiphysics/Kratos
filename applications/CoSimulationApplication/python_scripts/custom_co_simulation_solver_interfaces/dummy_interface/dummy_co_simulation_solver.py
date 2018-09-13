# co simulation imports
from co_simulation_base_solver import CoSimulationBaseSolver
import co_simulation_tools as tools


def Create(settings):
    return DummyCoSimulationSolver(settings)


class DummyCoSimulationSolver(CoSimulationBaseSolver):
    def __init__(self, custom_settings):
        super(DummyCoSimulationSolver, self).__init__(custom_settings)
        self.name = custom_settings["name"]
        pass

    def PrintInfo(self):
        print("This is an object of DummyCoSimulationSolver with name : ", self.name)

    def AdvanceInTime(self, current_time):
        return current_time + self.GetDeltaTime()

    def GetDeltaTime(self):
        return 0.1
