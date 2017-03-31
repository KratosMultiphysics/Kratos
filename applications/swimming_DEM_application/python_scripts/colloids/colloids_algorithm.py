from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
import swimming_DEM_algorithm
import swimming_DEM_procedures as SDP
import math
BaseAlgorithm = swimming_DEM_algorithm.Algorithm

class Algorithm(BaseAlgorithm):
    def __init__(self, pp):
        BaseAlgorithm.__init__(self, pp)
        self.cation_concentration_counter = self.GetCationConcentrationCounter()

    def SetBetaParamters(self):
        BaseAlgorithm.SetBetaParamters(self)
        self.pp.CFD_DEM.alpha = 1000.0
        self.pp.IntegrationScheme = 'TerminalVelocityScheme'

    def PerformInitialDEMStepOperations(self, time = None):
        if self.cation_concentration_counter.Tick():
            self.pp.concentration = self.pp.final_concentration + (self.pp.initial_concentration - self.pp.final_concentration) * math.exp(- self.pp.CFD_DEM.alpha * time)
            for node in self.spheres_model_part.Nodes:
                node.SetSolutionStepValue(CATION_CONCENTRATION, self.pp.concentration)
            if self.cation_concentration_counter.SuperTick(1000):
                print()
                print('Current cation concentration: ', str(self.pp.concentration))
                print()

    def GetCationConcentrationCounter(self):
        return SDP.Counter(self.pp.cation_concentration_frequence, 1)

    def DoSolveDEMVariable(self):
        self.pp.do_solve_dem = False
