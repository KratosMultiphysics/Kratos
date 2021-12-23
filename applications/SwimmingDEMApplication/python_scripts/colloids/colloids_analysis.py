import KratosMultiphysics as Kratos
import KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_procedures as SDP
import math
import KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_analysis as swimming_DEM_analysis
BaseAnalysis = swimming_DEM_analysis.SwimmingDEMAnalysis

class ColloidsAnalysis(BaseAnalysis):
    def __init__(self, varying_parameters = dict()):
        BaseAnalysis.__init__(self, varying_parameters)
        self.cation_concentration_counter = self.GetCationConcentrationCounter()

    def SetBetaParameters(self):
        BaseAnalysis.SetBetaParameters(self)
        Add = self.project_parameters.AddEmptyValue
        Add('PostCationConcentration').SetBool(True)
        Add('alpha').SetDouble(0.01)
        Add('initial_concentration').SetDouble(10)
        Add('final_concentration').SetDouble(10)
        Add('fluid_speed').SetDouble(10)
        Add('cation_concentration_frequence').SetDouble(10)
        self.project_parameters["custom_dem"]["translational_integration_scheme"].SetString('TerminalVelocityScheme')
        self.project_parameters["custom_dem"]["do_solve_dem"].SetBool(False)

    def PerformInitialDEMStepOperations(self, time = None):
        if self.cation_concentration_counter.Tick():
            self.concentration = (self.project_parameters['final_concentration'].GetDouble()
                                  + (self.project_parameters['initial_concentration'].GetDouble()
                                  - self.project_parameters['final_concentration'].GetDouble())
                                  * math.exp(- self.project_parameters['alpha'].GetDouble() * time
                                             / self.project_parameters["time_stepping"]["time_step"].GetDouble()))

            for node in self.spheres_model_part.Nodes:
                node.SetSolutionStepValue(Kratos.CATION_CONCENTRATION, self.concentration)

            if self.cation_concentration_counter.SuperTick(1000):
                Kratos.Logger.PrintInfo("SwimmingDEM", "")
                Kratos.Logger.PrintInfo("SwimmingDEM", 'Current cation concentration: ', str(self.concentration))
                Kratos.Logger.PrintInfo("SwimmingDEM", "")

    def GetCationConcentrationCounter(self):
        return SDP.Counter(self.project_parameters['cation_concentration_frequence'].GetInt(), 1)
