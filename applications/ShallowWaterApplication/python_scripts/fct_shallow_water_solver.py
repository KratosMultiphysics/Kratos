# importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

## Import base class file
from KratosMultiphysics.ShallowWaterApplication.stabilized_shallow_water_solver import StabilizedShallowWaterSolver

def CreateSolver(model, custom_settings):
    return FCTShallowWaterSolver(model, custom_settings)

class FCTShallowWaterSolver(StabilizedShallowWaterSolver):

    def Initialize(self):
        super().Initialize()
        self.main_model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 2)
        self.fct_utility = SW.AlgebraicFluxCorrectionUtility(self.main_model_part, self.settings["afc_parameters"])

    def InitializeSolutionStep(self):
        if self._TimeBufferIsInitialized():
            self.fct_utility.ExecuteInitializeLowOrderStep()
            self.__execute_before_low_order_step()
            super().InitializeSolutionStep()

    def SolveSolutionStep(self):
        if self._TimeBufferIsInitialized():
            is_converged = self.solver.SolveSolutionStep()
            self.fct_utility.ExecuteFinalizeLowOrderStep()
            self.fct_utility.ExecuteInitializeHighOrderStep()
            self.__execute_before_high_order_step()
            is_converged = self.solver.SolveSolutionStep()
            return is_converged
        else:
            return True

    def FinalizeSolutionStep(self):
        if self._TimeBufferIsInitialized():
            super().FinalizeSolutionStep()
            self.fct_utility.ExecuteFinalizeHighOrderStep()

    @classmethod
    def GetDefaultParameters(cls):
        afc_parameters = KM.Parameters("""{
            "name"              : "algebraic_flux_correction_utility",
            "rebuild_level"     : 0,
            "limiting_variable" : "VARIABLE_NAME"
        }""")
        default_parameters = super().GetDefaultParameters()
        default_parameters.AddValue("afc_parameters", afc_parameters)
        return default_parameters

    def __execute_before_low_order_step(self):
        dry_height = self.main_model_part.ProcessInfo.GetValue(SW.DRY_HEIGHT)
        SW.ShallowWaterUtilities().ResetDryDomain(self.main_model_part, dry_height)
        SW.ShallowWaterUtilities().IdentifyWetDomain(self.main_model_part, KM.ACTIVE, dry_height)

    def __execute_before_high_order_step(self):
        dry_height = 10 * self.main_model_part.ProcessInfo.GetValue(SW.DRY_HEIGHT)
        SW.ShallowWaterUtilities().IdentifyWetDomain(self.main_model_part, KM.ACTIVE, dry_height)
