import typing
from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.DigitalTwinApplication as KratosDT
import KratosMultiphysics.HDF5Application as KratosHDF5

from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.DigitalTwinApplication.sensor_sensitivity_solvers.sensor_sensitivity_adjoint_static_solver import SensorSensitivityAdjointStaticSolver
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import GetSensors
from KratosMultiphysics.DigitalTwinApplication.utilities.expression_utils import GetContainerExpression
from KratosMultiphysics.DigitalTwinApplication.utilities.expression_utils import GetContainerExpressionType
from KratosMultiphysics.HDF5Application.core.file_io import CreateHDF5File

class SensorSensitivityStaticAnalysis(AnalysisStage):
    def __init__(self, model: Kratos.Model, project_parameters: Kratos.Parameters):
        super().__init__(model, project_parameters)
        self.listof_sensors: 'list[KratosDT.Sensors.Sensor]' = {}

    def Initialize(self):
        super().Initialize()

        default_sensor_settings = Kratos.Parameters("""{
            "perturbation_size"      : 1e-8,
            "adapt_perturbation_size": true,
            "list_of_sensors"        : []
        }""")

        sensor_settings = self.project_parameters["sensor_settings"]
        sensor_settings.ValidateAndAssignDefaults(default_sensor_settings)

        model_part: Kratos.ModelPart = self._GetSolver().GetComputingModelPart()
        model_part.ProcessInfo[KratosDT.PERTURBATION_SIZE] = sensor_settings["perturbation_size"].GetDouble()
        model_part.ProcessInfo[KratosDT.ADAPT_PERTURBATION_SIZE] = sensor_settings["adapt_perturbation_size"].GetBool()
        self.listof_sensors = GetSensors(model_part, sensor_settings["list_of_sensors"].values())

    def _CreateSolver(self) -> SensorSensitivityAdjointStaticSolver:
        return SensorSensitivityAdjointStaticSolver(self.model, self.project_parameters["solver_settings"])

    def _GetSimulationName(self) -> str:
        return "::[SensorSpecificationAnalysis]:: "

    def RunSolutionLoop(self):
        # clear sensor data containers
        for sensor in self.GetListOfSensors():
            sensor.ClearNodalExpressions()
            sensor.ClearConditionExpressions()
            sensor.ClearElementExpressions()

        computing_model_part: Kratos.ModelPart = self._GetSolver().GetComputingModelPart()

        output_path = Path("sensor_sensitivities.h5")
        hdf5_parameters = Kratos.Parameters("""
                {
                    "file_name": "sensor_sensitivities.h5",
                    "file_access_mode": "truncate"
                }""")
        if output_path.is_file:
            hdf5_parameters["file_access_mode"].SetString("read_write")

        h5_file = CreateHDF5File(computing_model_part, hdf5_parameters)
        for sensor in self.listof_sensors:
            computing_model_part.ProcessInfo[Kratos.STEP] += 1
            computing_model_part.ProcessInfo[KratosDT.SENSOR_NAME] = sensor.GetName()
            computing_model_part.ProcessInfo[KratosDT.SENSOR_LOCATION] = sensor.GetLocation()

            h5_path = f"/SensitivityData/{sensor.__class__.__name__}/{sensor.GetName()}"
            sensitivity_variables: 'dict[Kratos.Globals.DataLocation, list[typing.Union[Kratos.DoubleVariable, Kratos.Array1DVariable3]]]' = self._GetSolver().GetSensitivtyVariables()
            if not h5_file.HasPath(h5_path):
                # recursively creat group
                h5_path_name_data = h5_path.split("/")[1:]
                for i, _ in enumerate(h5_path_name_data):
                    current_long_path = "/" + "/".join(h5_path_name_data[0:i+1])
                    if not h5_file.HasPath(current_long_path):
                        h5_file.CreateGroup(current_long_path)

                self._GetSolver().SetSensor(sensor)
                self.InitializeSolutionStep()
                sensor.InitializeSolutionStep()

                self._GetSolver().SolveSolutionStep()

                # calling calculate value to set the sensor value
                sensor.SetSensorValue(sensor.CalculateValue(computing_model_part))

                sensor.FinalizeSolutionStep()
                self.FinalizeSolutionStep()
                self.UpdateSensorSensitivities(sensor)

                # now save everything for reloading
                attribs = Kratos.Parameters("""{
                    "value": 0.0
                }""")
                attribs["value"].SetDouble(sensor.GetSensorValue())
                h5_file.WriteAttribute(h5_path, attribs)

                for data_location, variables in sensitivity_variables.items():
                    h5_data_params = Kratos.Parameters("""{
                        "prefix": ""
                    }""")
                    h5_data_params["prefix"].SetString(f"{h5_path}/{data_location.name}/")
                    h5_file_data = KratosHDF5.ExpressionIO(h5_data_params, h5_file)
                    if data_location in [Kratos.Globals.DataLocation.NodeHistorical, Kratos.Globals.DataLocation.NodeNonHistorical]:
                        for variable in variables:
                            h5_file_data.Write(f"{variable.Name()}", sensor.GetNodalExpression(f"{variable.Name()}").Clone())
                    elif data_location == Kratos.Globals.DataLocation.Condition:
                        for variable in variables:
                            h5_file_data.Write(f"{variable.Name()}", sensor.GetConditionExpression(f"{variable.Name()}").Clone())
                    elif data_location == Kratos.Globals.DataLocation.Element:
                        for variable in variables:
                            h5_file_data.Write(f"{variable.Name()}", sensor.GetElementExpression(f"{variable.Name()}").Clone())
                    else:
                        raise RuntimeError(f"Unsupported data location type = {data_location}.")

                self.OutputSolutionStep()
            else:
                attribs = h5_file.ReadAttribute(h5_path)
                sensor.SetSensorValue(attribs["value"].GetDouble())
                sensitivity_model_part: Kratos.ModelPart = self._GetSolver().GetSensitivityModelPart()
                for data_location, variables in sensitivity_variables.items():
                    h5_data_params = Kratos.Parameters("""{
                        "prefix": ""
                    }""")
                    h5_data_params["prefix"].SetString(f"{h5_path}/{data_location.name}/")
                    h5_file_data = KratosHDF5.ExpressionIO(h5_data_params, h5_file)

                    exp_type = GetContainerExpressionType(data_location)
                    for variable in variables:
                        expression = exp_type(sensitivity_model_part)
                        h5_file_data.Read(f"{variable.Name()}", expression)
                        if isinstance(expression, Kratos.Expression.NodalExpression):
                            sensor.AddNodalExpression(variable.Name(), expression)
                        elif isinstance(expression, Kratos.Expression.ConditionExpression):
                            sensor.AddConditionExpression(variable.Name(), expression)
                        elif isinstance(expression, Kratos.Expression.ElementExpression):
                            sensor.AddElementExpression(variable.Name(), expression)
                        else:
                            raise RuntimeError("Unsupported expression type.")

        h5_file.Close()

    def GetListOfSensors(self) -> 'list[KratosDT.Sensors.Sensor]':
        return self.listof_sensors

    def PrintAnalysisStageProgressInformation(self):
        process_info: Kratos.ProcessInfo = self._GetSolver().GetComputingModelPart().ProcessInfo
        Kratos.Logger.PrintInfo(self._GetSimulationName(), f"STEP: {process_info[Kratos.STEP]} - {process_info[KratosDT.SENSOR_NAME]}")

    def UpdateSensorSensitivities(self, sensor: KratosDT.Sensors.Sensor) -> None:
        sensitivity_model_part: Kratos.ModelPart = self._GetSolver().GetSensitivityModelPart()
        sensitivity_variables: 'dict[Kratos.Globals.DataLocation, list[typing.Union[Kratos.DoubleVariable, Kratos.Array1DVariable3]]]' = self._GetSolver().GetSensitivtyVariables()
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
