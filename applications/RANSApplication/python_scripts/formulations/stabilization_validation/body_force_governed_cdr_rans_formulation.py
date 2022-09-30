# import kratos
import KratosMultiphysics as Kratos

# import RANS
import KratosMultiphysics.RANSApplication as KratosRANS

# import formulation interface
from KratosMultiphysics.RANSApplication.formulations.scalar_rans_formulation import ScalarRansFormulation

class BodyForceGovernedCDRRansFormulation(ScalarRansFormulation):
    def __init__(self, model_part, settings, deprecated_settings):
        default_settings = Kratos.Parameters(r'''
        {
            "formulation_name": "body_force_governed_cdr",
            "stabilization_method": "residual_based_flux_corrected",
            "body_force_governed_cdr_solver_settings": {}
        }''')

        settings.ValidateAndAssignDefaults(default_settings)

        self.stabilization_method = settings["stabilization_method"].GetString()
        self.SetStabilizationMethod(self.stabilization_method)

        super().__init__(
            model_part,
            settings["body_force_governed_cdr_solver_settings"],
            deprecated_settings)

    def GetSolvingVariable(self):
        return KratosRANS.VELOCITY_POTENTIAL

    def GetElementNamePrefix(self):
        return "RansBodyForceGovernedCDR"

    def GetConditionNamePrefix(self):
        return ""

