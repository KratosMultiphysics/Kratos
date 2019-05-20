from __future__ import print_function, absolute_import, division  # makes backward compatible with python 2.6 and 2.7

import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools


class CoSimulationComponent(object):
    def __init__(self):
        pass

    def Initialize(self):
        pass

    def Finalize(self):
        pass

    def InitializeSolutionStep(self):
        pass

    def Predict(self):
        pass

    def SolveSolutionStep(self):
        pass

    def FinalizeSolutionStep(self):
        pass

    def OutputSolutionStep(self):
        pass

    def Check(self):
        pass

    def PrintInfo(self):
        cs_tools.PrintInfo("The component ", self.__class__.__name__)