# import kratos
import KratosMultiphysics as Kratos

# import RANS
import KratosMultiphysics.RANSApplication as KratosRANS

# import formulation interface
from KratosMultiphysics.RANSApplication.formulations.scalar_rans_formulation import ScalarRansFormulation

class CircularConvectionRansFormulation(ScalarRansFormulation):
    def __init__(self, model_part, settings):
        default_settings = Kratos.Parameters(r'''
        {
            "formulation_name": "circular_convection",
            "stabilization_method": "residual_based_flux_corrected",
            "circular_convection_solver_settings": {}
        }''')

        settings.ValidateAndAssignDefaults(default_settings)

        self.stabilization_method = settings["stabilization_method"].GetString()
        self.SetStabilizationMethod(self.stabilization_method)

        super().__init__(
            model_part,
            settings["circular_convection_solver_settings"])

    def GetSolvingVariable(self):
        return KratosRANS.VELOCITY_POTENTIAL

    def GetElementNamePrefix(self):
        return "RansCircularConvection"

    def GetConditionNamePrefix(self):
        return ""

