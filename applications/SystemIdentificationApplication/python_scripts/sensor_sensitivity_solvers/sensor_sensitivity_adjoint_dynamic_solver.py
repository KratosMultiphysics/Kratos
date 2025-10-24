import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_adjoint_dynamic_solver import StructuralMechanicsAdjointDynamicSolver
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.SystemIdentificationApplication.utilities.expression_utils import ExpressionUnionType
from KratosMultiphysics.SystemIdentificationApplication.sensor_sensitivity_solvers.response_sensitivity_interface import ResponseSensitivityInterface

class SensorSensitivityAdjointDynamicSolver(StructuralMechanicsAdjointDynamicSolver, ResponseSensitivityInterface):
    def __init__(self, model: Kratos.Model, custom_settings: Kratos.Parameters):
        super().__init__(model, custom_settings)
        self.sensitivity_builder_initialized = False

        default_response_function_settings = Kratos.Parameters("""{
            "perturbation_size"      : 1e-8,
            "adapt_perturbation_size": true,
            "bossak_alpha": -0.3
        }""")
        self.settings["response_function_settings"].ValidateAndAssignDefaults(default_response_function_settings)

    def Initialize(self) -> None:
        model_part = self.GetComputingModelPart()

        response_function_settings: Kratos.Parameters = self.settings["response_function_settings"]
        model_part.ProcessInfo[KratosSI.PERTURBATION_SIZE] = response_function_settings["perturbation_size"].GetDouble()
        model_part.ProcessInfo[KratosSI.ADAPT_PERTURBATION_SIZE] = response_function_settings["adapt_perturbation_size"].GetBool()

        # create sensitivity builder
        print("constructing a SensitivityBuilder")
        self.sensitivity_builder = Kratos.SensitivityBuilder(self.settings["sensitivity_settings"], self.main_model_part, None)

        # skip calling Initialize of StructuralMechanicsAdjointDynamicSolver and call the its base class Initialize
        super(StructuralMechanicsAdjointDynamicSolver, self).Initialize()

        # set the strategy to only build LHS once, and build RHS for all other solves
        self._GetSolutionStrategy().SetRebuildLevel(1)

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

    def Predict(self):
        pass 

    def SolveSolutionStep(self) -> bool:
        return super(StructuralMechanicsAdjointDynamicSolver, self).SolveSolutionStep()

    def GetSensitivities(self) -> 'dict[SupportedSensitivityFieldVariableTypes, ExpressionUnionType]':
        sensitivity_settings = self.settings["sensitivity_settings"]
        sensitivities_dict: 'dict[SupportedSensitivityFieldVariableTypes, ExpressionUnionType]' = {}
        self.__ReadSensitivities(sensitivities_dict, sensitivity_settings["nodal_solution_step_sensitivity_variables"].GetStringArray(), Kratos.Expression.NodalExpression, lambda x, y: Kratos.Expression.VariableExpressionIO.Read(x, y, True))
        self.__ReadSensitivities(sensitivities_dict, sensitivity_settings["condition_data_value_sensitivity_variables"].GetStringArray(), Kratos.Expression.ConditionExpression, lambda x, y: Kratos.Expression.VariableExpressionIO.Read(x, y))
        self.__ReadSensitivities(sensitivities_dict, sensitivity_settings["element_data_value_sensitivity_variables"].GetStringArray(), Kratos.Expression.ElementExpression, lambda x, y: Kratos.Expression.VariableExpressionIO.Read(x, y))
        return sensitivities_dict

    def __ReadSensitivities(self, sensitivity_dict: 'dict[SupportedSensitivityFieldVariableTypes, ExpressionUnionType]', list_of_variable_names: 'list[str]', expression_type, read_lambda) -> None:
        for sensitivity_variable_name in list_of_variable_names:
 
            sensitivity_variable_name = Kratos.SensitivityUtilities.GetSensitivityVariableName(Kratos.KratosGlobals.GetVariable(sensitivity_variable_name))
            sensitivity_variable = Kratos.KratosGlobals.GetVariable(sensitivity_variable_name)
            exp = expression_type(self.GetSensitivityModelPart())
            read_lambda(exp, sensitivity_variable)
            sensitivity_dict[sensitivity_variable] = exp.Clone()

    def _CreateScheme(self):
        settings = Kratos.Parameters("""{
            "name"         : "adjoint_bossak",
            "scheme_type"  : "bossak"
        }""")
        settings.AddDouble("alpha_bossak", self.settings["response_function_settings"]["bossak_alpha"].GetDouble())
        return Kratos.ResidualBasedAdjointBossakScheme(settings, None)
    
    # def _ComputeDeltaTime(self):
    #     """This function returns the delta time
    #     Note that the analysis does not support the automatic time stepping
    #     Also note that as it requires to go backwards in time, here we return minus the user time step
    #     """
    #     # if self.settings["time_stepping"]["automatic_time_step"].GetBool():
    #     #     raise Exception("Automatic time stepping is not supported by adjoint fluid solver.")

    #     delta_time = self.settings["time_stepping"]["time_step"].GetDouble()
    #     if delta_time > 0.0:
    #         # Expected positive time_step is reversed to advance backwards in time
    #         delta_time *= -1.0
    #     else:
    #         # TODO: Remove this check after the backwards compatibility period
    #         # In order to keep backwards compatibility, we throw a warning if a negative time_step is provided
    #         Kratos.Logger.PrintWarning("Setting a negative 'time_step' is no longer needed by 'SensorSensitivityAdjointDynamicSolver'.")

    #     return delta_time
    
    # def AdvanceInTime(self, current_time):
    #     dt = self._ComputeDeltaTime()
    #     print("dt is ", dt)
    #     print("current time is ", current_time)
    #     new_time = current_time + dt
    #     print("new time is ", new_time)
    #     self.main_model_part.ProcessInfo[Kratos.STEP] -= 1
    #     print("current step is ", self.main_model_part.ProcessInfo[Kratos.STEP])
    #     self.main_model_part.CloneTimeStep(new_time)
    #     return new_time