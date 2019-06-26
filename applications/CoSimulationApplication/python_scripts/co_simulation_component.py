from __future__ import print_function, absolute_import, division  # makes backward compatible with python 2.6 and 2.7

import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools


class CoSimulationComponent(object):
    def __init__(self):
        self.initialized = False
        self.initialized_solution_step = False

    def Initialize(self):
        if self.initialized:
            raise Exception("Already initialized")
        else:
            self.initialized = True

    def Finalize(self):
        if self.initialized_solution_step:
            raise Exception("Solution step not finalized")
        if self.initialized:
            self.initialized = False
        else:
            raise Exception("Not initialized")

    def InitializeSolutionStep(self):
        if self.initialized:
            if self.initialized_solution_step:
                raise Exception("Already solution step initialized")
            else:
                self.initialized_solution_step = True
        else:
            raise Exception("Not initialized")

    def FinalizeSolutionStep(self):
        if self.initialized:
            if self.initialized_solution_step:
                self.initialized_solution_step = False
            else:
                raise Exception("Solution step not initialized")
        else:
            raise Exception("Not initialized")

    def OutputSolutionStep(self):
        pass

    def Check(self):
        pass

    def PrintInfo(self):
        cs_tools.PrintInfo("The component ", self.__class__.__name__)