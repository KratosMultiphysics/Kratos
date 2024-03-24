import typing

import KratosMultiphysics as Kratos
import KratosMultiphysics.DigitalTwinApplication as KratosDT

from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.DigitalTwinApplication.sensor_sensitivity_solvers.sensor_sensitivity_adjoint_static_solver import SensorSensitivityAdjointStaticSolver
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import GetSensors
from KratosMultiphysics.DigitalTwinApplication.utilities.expression_utils import GetContainerExpression
from KratosMultiphysics.DigitalTwinApplication.utilities.expression_utils import ExpressionDataLocation

class SensorSensitivityStaticAnalysis(AnalysisStage):
    def __init__(self, model: Kratos.Model, project_parameters: Kratos.Parameters):
        super().__init__(model, project_parameters)
        self.listof_sensors: 'list[KratosDT.Sensors.Sensor]' = {}

    def Initialize(self):
        super().Initialize()

        default_sensor_settings = Kratos.Parameters("""{
            "perturbation_size"            : 1e-8,
            "adapt_perturbation_size"      : true,
            "force_calculate_sensitivities": true,
            "list_of_sensors"              : []
        }""")

        sensor_settings = self.project_parameters["sensor_settings"]
        sensor_settings.ValidateAndAssignDefaults(default_sensor_settings)

        model_part: Kratos.ModelPart = self._GetSolver().GetComputingModelPart()
        model_part.ProcessInfo[KratosDT.PERTURBATION_SIZE] = sensor_settings["perturbation_size"].GetDouble()
        model_part.ProcessInfo[KratosDT.ADAPT_PERTURBATION_SIZE] = sensor_settings["adapt_perturbation_size"].GetBool()
        self.listof_sensors = GetSensors(model_part, sensor_settings["list_of_sensors"].values())
        self.force_calculate_sensitivities = sensor_settings["force_calculate_sensitivities"].GetBool()

    def _CreateSolver(self) -> SensorSensitivityAdjointStaticSolver:
        return SensorSensitivityAdjointStaticSolver(self.model, self.project_parameters["solver_settings"])

    def _GetSimulationName(self) -> str:
        return "::[SensorAnalysis]:: "

    def RunSolutionLoop(self):
        # clear sensor data containers
        for sensor in self.GetListOfSensors():
            sensor.ClearNodalExpressions()
            sensor.ClearConditionExpressions()
            sensor.ClearElementExpressions()

        computing_model_part: Kratos.ModelPart = self._GetSolver().GetComputingModelPart()

        for sensor in self.listof_sensors:
            computing_model_part.ProcessInfo[Kratos.STEP] += 1
            computing_model_part.ProcessInfo[KratosDT.SENSOR_NAME] = sensor.GetName()

            self._GetSolver().SetSensor(sensor)
            self.InitializeSolutionStep()
            sensor.InitializeSolutionStep()

            self._GetSolver().SolveSolutionStep()

            # calling calculate value to set the sensor value
            sensor.SetSensorValue(sensor.CalculateValue(computing_model_part))

            sensor.FinalizeSolutionStep()
            self.FinalizeSolutionStep()
            self.UpdateSensorSensitivities(sensor)

            self.OutputSolutionStep()

    def GetListOfSensors(self) -> 'list[KratosDT.Sensors.Sensor]':
        return self.listof_sensors

    def PrintAnalysisStageProgressInformation(self):
        process_info: Kratos.ProcessInfo = self._GetSolver().GetComputingModelPart().ProcessInfo
        Kratos.Logger.PrintInfo(self._GetSimulationName(), f"STEP: {process_info[Kratos.STEP]} - {process_info[KratosDT.SENSOR_NAME]}")

    def UpdateSensorSensitivities(self, sensor: KratosDT.Sensors.Sensor) -> None:
        sensitivity_model_part: Kratos.ModelPart = self._GetSolver().GetSensitivityModelPart()
        sensitivity_variables: 'dict[ExpressionDataLocation, list[typing.Union[Kratos.DoubleVariable, Kratos.Array1DVariable3]]]' = self._GetSolver().GetSensitivtyVariables()
        for data_location, variables in sensitivity_variables.items():
            for variable in variables:
                # first read the variable
                expression = GetContainerExpression(sensitivity_model_part, data_location, variable)
                if isinstance(expression, Kratos.Expression.NodalExpression):
                    sensor.AddNodalExpression(variable.Name(), expression)
                elif isinstance(expression, Kratos.Expression.ConditionExpression):
                    sensor.AddConditionExpression(variable.Name(), expression)
                elif isinstance(expression, Kratos.Expression.ElementExpression):
                    sensor.AddElementExpression(variable.Name(), expression)
                else:
                    raise RuntimeError("Unsupported expression type.")

    def GetSensitivityModelPart(self) -> Kratos.ModelPart:
        return self._GetSolver().GetSensitivityModelPart()

if __name__ == "__main__":
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default ProjectParameters (read from "ProjectParameters.json"):\n'
        err_msg += '    "python3 structural_mechanics_analysis.py"\n'
        err_msg += '- With custom ProjectParameters:\n'
        err_msg += '    "python3 structural_mechanics_analysis.py CustomProjectParameters.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        project_parameters_file_name = argv[1]
    else: # using default name
        project_parameters_file_name = "ProjectParameters.json"

    with open(project_parameters_file_name,'r') as parameter_file:
        parameters = Kratos.Parameters(parameter_file.read())

    model = Kratos.Model()
    simulation = SensorSensitivityStaticAnalysis(model, parameters)
    simulation.Run()
