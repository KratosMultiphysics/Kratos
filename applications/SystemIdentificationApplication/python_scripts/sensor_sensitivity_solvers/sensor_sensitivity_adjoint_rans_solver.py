import KratosMultiphysics as Kratos
import KratosMultiphysics.RANSApplication as KratosRANS
from KratosMultiphysics.RANSApplication.formulations.utilities import CalculateNormalsOnConditions
from KratosMultiphysics.RANSApplication.adjoint_rans_solver import AdjointRANSSolver
from KratosMultiphysics.SystemIdentificationApplication.sensor_sensitivity_solvers.response_sensitivity_interface import ResponseSensitivityInterface

class SensorSensitivityAdjointRANSSolver(AdjointRANSSolver, ResponseSensitivityInterface):
    def __init__(self, model, custom_settings):
        super().__init__(model, custom_settings)
        self.sensitivity_builder_initialized = False

    def AddVariables(self):
        super().AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(KratosRANS.VELOCITY_SENSITIVITY)

    def SetResponseFunction(self, response_function: Kratos.AdjointResponseFunction) -> None:
        self.response_function = response_function
        self.GetSensitivityBuilder().SetResponseFunction(response_function)
        self._GetScheme().SetResponseFunction(response_function)
        if not self.sensitivity_builder_initialized:
            self.GetSensitivityBuilder().Initialize()
            self.sensitivity_builder_initialized = True

    def Initialize(self):
        self.communicator = None

        domain_size = self.main_model_part.ProcessInfo[Kratos.DOMAIN_SIZE]
        CalculateNormalsOnConditions(self.main_model_part)
        Kratos.NormalCalculationUtils().CalculateNormalShapeDerivativesOnSimplex(self.main_model_part.Conditions, domain_size)

        # Construct and set the solution strategy
        solution_strategy = self._GetSolutionStrategy()
        solution_strategy.SetEchoLevel(self.settings["echo_level"].GetInt())

        # Initialize the strategy and adjoint utilities
        solution_strategy.Initialize()

        # create sensitivity builder
        print(self.settings)
        self.sensitivity_builder = Kratos.SensitivityBuilder(self.adjoint_settings["sensitivity_settings"], self.main_model_part, None)

        Kratos.Logger.PrintInfo(self.__class__.__name__, "Solver initialization finished.")

    def GetResponseFunction(self) -> Kratos.AdjointResponseFunction:
        if not hasattr(self, 'response_function'):
            return None
        return self.response_function

    def GetSensitivityModelPart(self) -> Kratos.ModelPart:
        return self.main_model_part.GetSubModelPart(self.settings["sensitivity_settings"]["sensitivity_model_part_name"].GetString())

    def SolveSolutionStep(self) -> bool:
        return super().SolveSolutionStep()
