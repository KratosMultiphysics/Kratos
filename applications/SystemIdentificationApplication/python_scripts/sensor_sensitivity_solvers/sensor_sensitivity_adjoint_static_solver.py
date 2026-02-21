import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_adjoint_static_solver import StructuralMechanicsAdjointStaticSolver
from KratosMultiphysics.SystemIdentificationApplication.sensor_sensitivity_solvers.response_sensitivity_interface import ResponseSensitivityInterface

class SensorSensitivityAdjointStaticSolver(StructuralMechanicsAdjointStaticSolver, ResponseSensitivityInterface):
    def __init__(self, model: Kratos.Model, custom_settings: Kratos.Parameters):
        super().__init__(model, custom_settings)
        self.sensitivity_builder_initialized = False

        default_response_function_settings = Kratos.Parameters("""{
            "perturbation_size"      : 1e-8,
            "adapt_perturbation_size": true
        }""")
        self.settings["response_function_settings"].ValidateAndAssignDefaults(default_response_function_settings)

    def Initialize(self) -> None:
        model_part = self.GetComputingModelPart()

        response_function_settings: Kratos.Parameters = self.settings["response_function_settings"]
        model_part.ProcessInfo[KratosSI.PERTURBATION_SIZE] = response_function_settings["perturbation_size"].GetDouble()
        model_part.ProcessInfo[KratosSI.ADAPT_PERTURBATION_SIZE] = response_function_settings["adapt_perturbation_size"].GetBool()

        # create sensitivity builder
        self.sensitivity_builder = Kratos.SensitivityBuilder(self.settings["sensitivity_settings"], self.main_model_part, None)

        # skip calling Initialize of StructuralMechanicsAdjointStaticSolver and call the its base class Initialize
        super(StructuralMechanicsAdjointStaticSolver, self).Initialize()

        Kratos.Logger.PrintInfo(self.__class__.__name__, "Finished initialization.")

    def SetResponseFunction(self, response_function: Kratos.AdjointResponseFunction) -> None:
        self.response_function = response_function
        self.sensitivity_builder.SetResponseFunction(response_function)
        self._GetScheme().SetResponseFunction(response_function)
        if not self.sensitivity_builder_initialized:
            self.sensitivity_builder.Initialize()
            self.sensitivity_builder_initialized = True

    def GetResponseFunction(self) -> Kratos.AdjointResponseFunction:
        return self._GetScheme().GetResponseFunction()

    def GetSensitivityModelPart(self) -> Kratos.ModelPart:
        return self.main_model_part.GetSubModelPart(self.settings["sensitivity_settings"]["sensitivity_model_part_name"].GetString())

    def SolveSolutionStep(self) -> bool:
        return super(StructuralMechanicsAdjointStaticSolver, self).SolveSolutionStep()

    def _CreateScheme(self):
        return Kratos.ResidualBasedAdjointStaticScheme(None)