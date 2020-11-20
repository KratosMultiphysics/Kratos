# importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

## Import base class file
from KratosMultiphysics.ShallowWaterApplication.stabilized_shallow_water_solver import StabilizedShallowWaterSolver

def CreateSolver(model, custom_settings):
    return FCTShallowWaterSolver(model, custom_settings)

class FCTShallowWaterSolver(StabilizedShallowWaterSolver):

    def AddVariables(self):
        super().AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(KM.INTERNAL_ENERGY)

    def Initialize(self):
        super().Initialize()
        self.main_model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 2)
        self.fct_utility = SW.AlgebraicFluxCorrectionUtility(self.main_model_part, self.settings["afc_parameters"])

    def InitializeSolutionStep(self):
        if self._TimeBufferIsInitialized():
            self.__execute_before_low_order_step()
            super().InitializeSolutionStep()

    def SolveSolutionStep(self):
        if self._TimeBufferIsInitialized():
            is_converged = self._GetSolutionStrategy().SolveSolutionStep()
            self.__execute_after_low_order_step()
            self.__execute_before_high_order_step()
            is_converged = self._GetSolutionStrategy().SolveSolutionStep()
            return is_converged
        else:
            return True

    def FinalizeSolutionStep(self):
        if self._TimeBufferIsInitialized():
            super().FinalizeSolutionStep()
            self.__execute_after_high_order_step()

    @classmethod
    def GetDefaultParameters(cls):
        afc_parameters = KM.Parameters("""{
            "name"               : "algebraic_flux_correction_utility",
            "limiting_variables" : ["VARIABLE_NAME"],
            "maximum_iterations" : 1,
            "rebuild_level"      : 0
        }""")
        default_parameters = super().GetDefaultParameters()
        default_parameters.AddValue("afc_parameters", afc_parameters)
        return default_parameters

    def __execute_before_low_order_step(self):
        self.main_model_part.ProcessInfo.SetValue(SW.LUMPED_MASS_FACTOR, 1.0)
        self.main_model_part.ProcessInfo.SetValue(SW.IS_MONOTONIC_CALCULATION, True)
        dry_height = self.main_model_part.ProcessInfo.GetValue(SW.DRY_HEIGHT)
        SW.ShallowWaterUtilities().ResetDryDomain(self.main_model_part, dry_height)
        SW.ShallowWaterUtilities().IdentifyWetDomain(self.main_model_part, KM.ACTIVE, dry_height)
        self.fct_utility.InitializeCorrection()

    def __execute_after_low_order_step(self):
        self.fct_utility.GetLowOrderValues()

    def __execute_before_high_order_step(self):
        self.main_model_part.ProcessInfo.SetValue(SW.LUMPED_MASS_FACTOR, 0.0)
        self.main_model_part.ProcessInfo.SetValue(SW.IS_MONOTONIC_CALCULATION, False)
        dry_height = 10 * self.main_model_part.ProcessInfo.GetValue(SW.DRY_HEIGHT)
        SW.ShallowWaterUtilities().IdentifyWetDomain(self.main_model_part, KM.ACTIVE, dry_height)

    def __execute_after_high_order_step(self):
        self.fct_utility.GetHighOrderValues()
        self.fct_utility.ApplyCorrection()
