import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI

from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import WorkFolderScope
from KratosMultiphysics.RANSApplication.adjoint_rans_solver import AdjointRANSSolver
from KratosMultiphysics.RANSApplication.adjoint_rans_analysis import AdjointRANSAnalysis
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import CreateSensors
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.SystemIdentificationApplication.sensor_sensitivity_solvers.response_sensitivity_analysis import ResponseSensitivityAnalysis

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
    return SystemIdentificationRANSAnalysis(model, parameters["settings"], optimization_problem)

class SensorSensitivityAdjointRANSSolver(AdjointRANSSolver):
    def __init__(self, model, custom_settings, sensor_group_name: str, p_coefficient: float, list_of_sensors: 'list[KratosSI.Sensors.Sensor]'):
        self.sensor_group_name = sensor_group_name
        self.p_coefficient = p_coefficient
        self.list_of_sensors = list_of_sensors
        super().__init__(model, custom_settings)

    def Initialize(self):
        measurement_residual = self.GetResponseFunction()
        domain_model_part = self.GetComputingModelPart()
        sensor_model_part = self.model.CreateModelPart(self.sensor_group_name)
        temp_list_of_sensors = CreateSensors(sensor_model_part, domain_model_part, [sensor.GetSensorParameters() for sensor in self.list_of_sensors])

        for i, sensor in enumerate(temp_list_of_sensors):
            sensor.GetNode().SetValue(KratosSI.SENSOR_MEASURED_VALUE, self.list_of_sensors[i].GetNode().GetValue(KratosSI.SENSOR_MEASURED_VALUE))
            measurement_residual.AddSensor(sensor)

        super().Initialize()
        measurement_residual.CalculateValue(domain_model_part)

    def _CreateResponseFunction(self) -> Kratos.AdjointResponseFunction:
        return KratosSI.Responses.MeasurementResidualResponseFunction(self.p_coefficient)

    def SolveSolutionStep(self):
        super().SolveSolutionStep()

class SensorSensitivityAdjointRansAnalysis(AdjointRANSAnalysis):
    def __init__(self, model, parameters, sensor_group_name: str, p_coefficient: float, list_of_sensors: 'list[KratosSI.Sensors.Sensor]'):
        self.sensor_group_name = sensor_group_name
        self.p_coefficient = p_coefficient
        self.list_of_sensors = list_of_sensors
        super().__init__(model, parameters)

    def _CreateSolver(self):
        return SensorSensitivityAdjointRANSSolver(self.model, self.project_parameters["solver_settings"], self.sensor_group_name, self.p_coefficient, self.list_of_sensors)

class SystemIdentificationRANSAnalysis(ResponseSensitivityAnalysis):
    def __init__(self, model: Kratos.Model, settings: Kratos.Parameters, optimization_problem: OptimizationProblem, p_coefficient: float):
        default_settings = Kratos.Parameters("""{
            "execution_policy_name"               : "",
            "adjoint_project_parameters_file_name": ""
        }""")

        self.model = model
        self.optimization_problem = optimization_problem
        self.p_coefficient = p_coefficient

        settings.ValidateAndAssignDefaults(default_settings)

        self.execution_policy = self.optimization_problem.GetExecutionPolicy(settings["execution_policy_name"].GetString())
        self.adjoint_project_parameters_file_name = settings["adjoint_project_parameters_file_name"].GetString()

        self.analysis = None
        self.sensor_group_name = ""
        self.list_of_sensors: 'list[KratosSI.Sensors.Sensor]' = []

    def Initialize(self):
        pass

    def Finalize(self):
        pass

    def SetSensors(self, sensor_group_name: str, list_of_sensors: 'list[KratosSI.Sensors.Sensor]') -> None:
        self.sensor_group_name = sensor_group_name
        self.list_of_sensors = list_of_sensors

    def CalculateGradient(self, _: Kratos.AdjointResponseFunction) -> None:
        with WorkFolderScope(self.execution_policy.GetPath()):
            model = Kratos.Model()
            with open(self.adjoint_project_parameters_file_name, "r") as file_input:
                parameters = Kratos.Parameters(file_input.read())

            self.analysis = SensorSensitivityAdjointRansAnalysis(model, parameters, self.sensor_group_name, self.p_coefficient, self.list_of_sensors)
            self.analysis.Run()

    def GetGradient(self, sensitivity_variable: SupportedSensitivityFieldVariableTypes, gradient_tensor_adaptor: Kratos.TensorAdaptors.DoubleTensorAdaptor) -> None:
        sub_model_part_name = self.analysis.project_parameters["solver_settings"]["sensitivity_settings"]["sensitivity_model_part_name"].GetString()
        model_part = self.analysis._GetSolver().main_model_part.GetSubModelPart(sub_model_part_name)

        if isinstance(gradient_tensor_adaptor.GetContainer(), Kratos.NodesArray):
            temp_ta = Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor(model_part.Nodes, sensitivity_variable)
            temp_ta.CollectData()
            gradient_tensor_adaptor.data[:] = temp_ta.data
        elif isinstance(gradient_tensor_adaptor.GetContainer(), Kratos.ConditionsArray):
            temp_ta = Kratos.TensorAdaptors.VariableTensorAdaptor(model_part.Conditions, sensitivity_variable)
            temp_ta.CollectData()
            gradient_tensor_adaptor.data[:] = temp_ta.data
        elif isinstance(gradient_tensor_adaptor.GetContainer(), Kratos.ElementsArray):
            temp_ta = Kratos.TensorAdaptors.VariableTensorAdaptor(model_part.Elements, sensitivity_variable)
            temp_ta.CollectData()
            gradient_tensor_adaptor.data[:] = temp_ta.data
        else:
            raise RuntimeError(f"Unsupported container type in the tensor adaptor [ tensor adaptor = {gradient_tensor_adaptor} ].")

