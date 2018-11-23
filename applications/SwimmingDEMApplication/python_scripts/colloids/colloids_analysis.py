from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
import swimming_DEM_procedures as SDP
import math
import swimming_DEM_analysis
BaseAnalysis = swimming_DEM_analysis.SwimmingDEMAnalysis

class ColloidsAnalysis(BaseAnalysis):
    def __init__(self, varying_parameters = dict()):
        BaseAnalysis.__init__(self, varying_parameters)
        self.cation_concentration_counter = self.GetCationConcentrationCounter()

    def SetBetaParameters(self):
        BaseAnalysis.SetBetaParameters(self)
        self.pp.CFD_DEM.alpha = 0.01
        self.pp["TranslationalIntegrationScheme"].SetString('TerminalVelocityScheme')
        self.pp.CFD_DEM.AddEmptyValue("basset_force_type").SetInt(0)
        self.pp.CFD_DEM.PostCationConcentration = True
        self.pp.initial_concentration = 10
        self.pp.final_concentration = 0.01
        self.pp.fluid_speed = 1e-6
        self.pp.cation_concentration_frequence = 1
        self.pp.CFD_DEM.drag_force_type = 9

    def PerformInitialDEMStepOperations(self, time = None):
        if self.cation_concentration_counter.Tick():
            self.pp.concentration = self.pp.final_concentration + (self.pp.initial_concentration - self.pp.final_concentration) * math.exp(- self.pp.CFD_DEM.alpha * time / self.pp.CFD_DEM["MaxTimeStep"].GetDouble())
            for node in self.spheres_model_part.Nodes:
                node.SetSolutionStepValue(CATION_CONCENTRATION, self.pp.concentration)
            if self.cation_concentration_counter.SuperTick(1000):
                print()
                print('Current cation concentration: ', str(self.pp.concentration))
                print()

    def GetCationConcentrationCounter(self):
        return SDP.Counter(self.pp.cation_concentration_frequence, 1)

    def DoSolveDEMVariable(self):
        self.do_solve_dem = False
