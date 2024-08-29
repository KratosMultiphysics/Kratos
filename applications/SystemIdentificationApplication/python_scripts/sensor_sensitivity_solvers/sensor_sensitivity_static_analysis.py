import typing

import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI

from KratosMultiphysics.SystemIdentificationApplication.sensor_sensitivity_solvers.sensor_sensitivity_analysis import SensorSensitivityAnalysis
from KratosMultiphysics.SystemIdentificationApplication.sensor_sensitivity_solvers.sensor_sensitivity_adjoint_static_solver import SensorSensitivityAdjointStaticSolver
from KratosMultiphysics.SystemIdentificationApplication.utilities.expression_utils import GetContainerExpression
from KratosMultiphysics.SystemIdentificationApplication.utilities.expression_utils import ExpressionDataLocation
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes

class SensorSensitivityStaticAnalysis(SensorSensitivityAnalysis):
    def __init__(self, model: Kratos.Model, project_parameters: Kratos.Parameters):
        super().__init__(model, project_parameters)
        self.list_of_sensors: 'list[KratosSI.Sensors.Sensor]' = []

    def Initialize(self):
        super().Initialize()

        default_sensor_settings = Kratos.Parameters("""{
            "perturbation_size"      : 1e-8,
            "adapt_perturbation_size": true
        }""")

        sensor_settings = self.project_parameters["sensor_settings"]
        sensor_settings.ValidateAndAssignDefaults(default_sensor_settings)

        model_part: Kratos.ModelPart = self._GetSolver().GetComputingModelPart()
        model_part.ProcessInfo[KratosSI.PERTURBATION_SIZE] = sensor_settings["perturbation_size"].GetDouble()
        model_part.ProcessInfo[KratosSI.ADAPT_PERTURBATION_SIZE] = sensor_settings["adapt_perturbation_size"].GetBool()

    def _CreateSolver(self) -> SensorSensitivityAdjointStaticSolver:
        return SensorSensitivityAdjointStaticSolver(self.model, self.project_parameters["solver_settings"])

    def _GetSimulationName(self) -> str:
        return "SensorSensitivityStaticAnalysis"

    def CalculateGradient(self, sensor: KratosSI.Sensors.Sensor) -> dict[str, ContainerExpressionTypes]:
        computing_model_part: Kratos.ModelPart = self._GetSolver().GetComputingModelPart()

        computing_model_part.ProcessInfo[Kratos.STEP] = sensor.GetNode().Id
        computing_model_part.ProcessInfo[KratosSI.SENSOR_NAME] = sensor.GetName()

        self._GetSolver().SetSensor(sensor)
        self.InitializeSolutionStep()
        sensor.InitializeSolutionStep()

        self._GetSolver().SolveSolutionStep()

        # calling calculate value to set the sensor value
        sensor.SetSensorValue(sensor.CalculateValue(computing_model_part))

        sensor.FinalizeSolutionStep()
        self.FinalizeSolutionStep()
        self.OutputSolutionStep()

        gradients: 'dict[str, ContainerExpressionTypes]' = {}

        sensitivity_model_part: Kratos.ModelPart = self._GetSolver().GetSensitivityModelPart()
        sensitivity_variables: 'dict[ExpressionDataLocation, list[typing.Union[Kratos.DoubleVariable, Kratos.Array1DVariable3]]]' = self._GetSolver().GetSensitivityVariables()
        for data_location, variables in sensitivity_variables.items():
            for variable in variables:
                # first read the variable
                gradients[variable.Name] = GetContainerExpression(sensitivity_model_part, data_location, variable)

        return gradients

    def PrintAnalysisStageProgressInformation(self):
        process_info: Kratos.ProcessInfo = self._GetSolver().GetComputingModelPart().ProcessInfo
        Kratos.Logger.PrintInfo(self._GetSimulationName(), f"Sensor Id = {process_info[Kratos.STEP]:05d}, Sensor name = {process_info[KratosSI.SENSOR_NAME]}")

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