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

    def AddVariables(self):
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.NORMAL)

        super().AddVariables()

    def GetSolvingVariable(self):
        return KratosRANS.VELOCITY_POTENTIAL

    def GetElementNamePrefix(self):
        return "RansCircularConvection"

    def GetConditionNamePrefix(self):
        return ""

    def SetConstants(self, settings):
        defaults = Kratos.Parameters('''{
            "circular_convection_constants": {
                "is_clock_wise_rotation": true,
                "rotation_center": [0.0, 0.0, 0.0]
            },
            "stabilization_constants":{
                "dynamic_tau"                       : 0.0,
                "upwind_operator_coefficient"       : 1.2,
                "positivity_preserving_coefficient" : 1.2
            }
        }''')

        settings.RecursivelyValidateAndAssignDefaults(defaults)

        process_info = self.GetBaseModelPart().ProcessInfo

        # stabilization parameters
        constants = settings["circular_convection_constants"]
        process_info.SetValue(KratosRANS.CIRCULAR_CONVECTION_ROTATION_CLOCKWISE, constants["is_clock_wise_rotation"].GetBool())
        process_info.SetValue(KratosRANS.CIRCULAR_CONVECTION_ROTATION_CENTER, constants["rotation_center"].GetVector())

        constants = settings["stabilization_constants"]
        if (self.is_steady_simulation):
            self.GetBaseModelPart().ProcessInfo.SetValue(Kratos.DYNAMIC_TAU, 0.0)
        else:
            self.GetBaseModelPart().ProcessInfo.SetValue(Kratos.DYNAMIC_TAU, constants["dynamic_tau"].GetDouble())

        # stabilization parameters
        process_info.SetValue(KratosRANS.RANS_STABILIZATION_DISCRETE_UPWIND_OPERATOR_COEFFICIENT, constants["upwind_operator_coefficient"].GetDouble())
        process_info.SetValue(KratosRANS.RANS_STABILIZATION_DIAGONAL_POSITIVITY_PRESERVING_COEFFICIENT, constants["positivity_preserving_coefficient"].GetDouble())
