# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# Import base class file
from KratosMultiphysics.FluidDynamicsApplication.navier_stokes_monolithic_solver import NavierStokesMonolithicSolver


def CreateSolver(model, custom_settings):
    return NavierStokesMonolithicIgaSolver(model, custom_settings)


class NavierStokesMonolithicIgaSolver(NavierStokesMonolithicSolver):

    @classmethod
    def GetDefaultParameters(cls):
        default_settings = KratosMultiphysics.Parameters("""
        {
            "skip_entities_replace_and_check": false
        }""")
        default_settings.AddMissingParameters(super().GetDefaultParameters())
        return default_settings

    def __init__(self, model, custom_settings):
        super().__init__(model, custom_settings)
        self._skip_entities_replace_and_check = self.settings["skip_entities_replace_and_check"].GetBool()

    def PrepareModelPart(self):
        if not self.is_restarted():
            materials_imported = self._SetPhysicalProperties()
            if not materials_imported:
                KratosMultiphysics.Logger.PrintWarning(
                    self.__class__.__name__,
                    "Material properties have not been imported. Check 'material_import_settings' in your ProjectParameters.json.")
            if self._enforce_element_and_conditions_replacement:
                self._ReplaceElementsAndConditions()
            self._SetAndFillBuffer()

        if not self._skip_entities_replace_and_check:
            self._ExecuteCheckAndPrepare()

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Model reading finished.")

    def Predict(self):
        scheme_type = self.settings["time_scheme"].GetString()
        if scheme_type in ("bdf2", "bdf2_higher_order_vms"):
            KratosMultiphysics.TimeDiscretization.BDF(2).ComputeAndSaveBDFCoefficients(
                self.GetComputingModelPart().ProcessInfo)
        self._GetSolutionStrategy().Predict()

    def GetComputingModelPart(self):
        if self._skip_entities_replace_and_check:
            return self.main_model_part
        return super().GetComputingModelPart()

    def _CreateScheme(self):
        if self.settings["time_scheme"].GetString() == "bdf2_higher_order_vms":
            return KratosCFD.BDF2HigherOrderVMSScheme()
        return super()._CreateScheme()

    def _SetTimeSchemeBufferSize(self):
        if self.settings["time_scheme"].GetString() == "bdf2_higher_order_vms":
            self.min_buffer_size = 3
        else:
            super()._SetTimeSchemeBufferSize()
