# co simulation imports
from co_simulation_base_solver import CoSimulationBaseSolver
import co_simulation_tools as tools


def Create(settings):
    return DummyCoSimulationSolver(settings)


class DummyCoSimulationSolver(CoSimulationBaseSolver):
    def __init__(self, custom_settings):
        super(DummyCoSimulationSolver, self).__init__(custom_settings)
        self.name = custom_settings["name"]

    def Initialize(self):
        self.InitializeIO()

    def PrintInfo(self):
        print( self.data_list.keys() )
        print(tools.bcolors.BLUE+"This is an object of DummyCoSimulationSolver with name : ", self.name + tools.bcolors.ENDC)

    def AdvanceInTime(self, current_time):
        return current_time + self.GetDeltaTime()

    def GetDeltaTime(self):
        return 0.2

    def Check(self):
        print(tools.bcolors.GREEN+"Check from dummy co simulation solver : CHECKED !"+tools.bcolors.ENDC)

    def _GetIOName(self):
        return "dummy"