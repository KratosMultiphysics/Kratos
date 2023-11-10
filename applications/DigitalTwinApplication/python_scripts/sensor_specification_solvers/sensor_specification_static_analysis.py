import typing
from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.DigitalTwinApplication as KratosDT
import KratosMultiphysics.HDF5Application as KratosHDF5

from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.DigitalTwinApplication.sensor_specification_solvers.sensor_specification_adjoint_static_solver import SensorSpecificationAdjointStaticSolver
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_specification_utils import GetSpecifications
from KratosMultiphysics.DigitalTwinApplication.utilities.expression_utils import GetContainerExpression
from KratosMultiphysics.DigitalTwinApplication.utilities.expression_utils import GetContainerExpressionType
from KratosMultiphysics.HDF5Application.core.file_io import CreateHDF5File

class SensorSpecificationStaticAnalysis(AnalysisStage):
    def __init__(self, model: Kratos.Model, project_parameters: Kratos.Parameters):
        super().__init__(model, project_parameters)
        self.dict_of_sensor_specifications: 'dict[KratosDT.Sensors.AdjointSensor, list[KratosDT.Sensors.SensorSpecification]]' = {}

    def Initialize(self):
        super().Initialize()

        default_sensor_settings = Kratos.Parameters("""{
            "sensor_type"           : "PLEASE_SPECIFY_SENSOR_TYPE",
            "base_settings"         : {},
            "list_of_specifications": []
        }""")

        # create the list of specifications
        sensor_setting: Kratos.Parameters
        for sensor_setting in self.project_parameters["sensor_settings"].values():
            sensor_setting.ValidateAndAssignDefaults(default_sensor_settings)

            sensor_type = sensor_setting["sensor_type"].GetString()

            if sensor_type == "adjoint_displacement_sensor":
                adjoint_sensor = KratosDT.Sensors.AdjointDisplacementSensor(self.model, sensor_setting["base_settings"])
                self.dict_of_sensor_specifications[adjoint_sensor] = GetSpecifications(self._GetSolver().GetComputingModelPart(), sensor_setting["list_of_specifications"].values())
            else:
                raise RuntimeError(f"Unsupported sensor_type = \"{sensor_type}\".")

    def _CreateSolver(self) -> SensorSpecificationAdjointStaticSolver:
        return SensorSpecificationAdjointStaticSolver(self.model, self.project_parameters["solver_settings"])

    def _GetSimulationName(self) -> str:
        return "::[SensorSpecificationAnalysis]:: "

    def RunSolutionLoop(self):
        # clear specification data containers
        for specification in self.GetListOfSpecifications():
            specification.ClearNodalExpressions()
            specification.ClearConditionExpressions()
            specification.ClearElementExpressions()

        computing_model_part: Kratos.ModelPart = self._GetSolver().GetComputingModelPart()

        output_path = Path("sensor_specification_sensitivities.h5")
        hdf5_parameters = Kratos.Parameters("""
                {
                    "file_name": "sensor_specification_sensitivities.h5",
                    "file_access_mode": "truncate"
                }""")
        if output_path.is_file:
            hdf5_parameters["file_access_mode"].SetString("read_write")

        h5_file = CreateHDF5File(computing_model_part, hdf5_parameters)
        for adjoint_sensor, list_of_specifications in self.dict_of_sensor_specifications.items():
            self._GetSolver().SetSensor(adjoint_sensor)

            adjoint_sensor.Initialize()

            for specification in list_of_specifications:
                computing_model_part.ProcessInfo[Kratos.STEP] += 1
                computing_model_part.ProcessInfo[KratosDT.SENSOR_NAME] = specification.GetType()
                computing_model_part.ProcessInfo[KratosDT.SENSOR_ID] = specification.Id
                computing_model_part.ProcessInfo[KratosDT.SENSOR_LOCATION] = specification.GetLocation()

                h5_path = f"/SensitivityData/{adjoint_sensor.__class__.__name__}/{specification.GetType()}/{specification.Id}"
                sensitivity_variables: 'dict[Kratos.Globals.DataLocation, list[typing.Union[Kratos.DoubleVariable, Kratos.Array1DVariable3]]]' = self._GetSolver().GetSensitivtyVariables()
                if not h5_file.HasPath(h5_path):
                    # recursively creat group
                    h5_path_name_data = h5_path.split("/")[1:]
                    for i, _ in enumerate(h5_path_name_data):
                        current_long_path = "/" + "/".join(h5_path_name_data[0:i+1])
                        if not h5_file.HasPath(current_long_path):
                            h5_file.CreateGroup(current_long_path)

                    adjoint_sensor.SetSensorSpecification(specification)

                    self.InitializeSolutionStep()
                    adjoint_sensor.InitializeSolutionStep()

                    self._GetSolver().SolveSolutionStep()
                    specification.SetSensorValue(adjoint_sensor.CalculateValue(computing_model_part))


                    adjoint_sensor.FinalizeSolutionStep()
                    self.FinalizeSolutionStep()
                    self.UpdateSpecification(specification)

                    # now save everything for reloading
                    attribs = Kratos.Parameters("""{
                        "value": 0.0
                    }""")
                    attribs["value"].SetDouble(specification.GetSensorValue())
                    h5_file.WriteAttribute(h5_path, attribs)

                    for data_location, variables in sensitivity_variables.items():
                        h5_data_params = Kratos.Parameters("""{
                            "prefix": ""
                        }""")
                        h5_data_params["prefix"].SetString(f"{h5_path}/{data_location.name}/")
                        h5_file_data = KratosHDF5.ExpressionIO(h5_data_params, h5_file)
                        if data_location in [Kratos.Globals.DataLocation.NodeHistorical, Kratos.Globals.DataLocation.NodeNonHistorical]:
                            for variable in variables:
                                h5_file_data.Write(f"{variable.Name()}", specification.GetNodalExpression(f"{variable.Name()}").Clone())
                        elif data_location == Kratos.Globals.DataLocation.Condition:
                            for variable in variables:
                                h5_file_data.Write(f"{variable.Name()}", specification.GetConditionExpression(f"{variable.Name()}").Clone())
                        elif data_location == Kratos.Globals.DataLocation.Element:
                            for variable in variables:
                                h5_file_data.Write(f"{variable.Name()}", specification.GetElementExpression(f"{variable.Name()}").Clone())
                        else:
                            raise RuntimeError(f"Unsupported data location type = {data_location}.")

                    self.OutputSolutionStep()
                else:
                    attribs = h5_file.ReadAttribute(h5_path)
                    specification.SetSensorValue(attribs["value"].GetDouble())
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
                                specification.AddNodalExpression(variable.Name(), expression)
                            elif isinstance(expression, Kratos.Expression.ConditionExpression):
                                specification.AddConditionExpression(variable.Name(), expression)
                            elif isinstance(expression, Kratos.Expression.ElementExpression):
                                specification.AddElementExpression(variable.Name(), expression)
                            else:
                                raise RuntimeError("Unsupported expression type.")

        h5_file.Close()

    def GetListOfSpecifications(self) -> 'list[KratosDT.Sensors.SensorSpecification]':
        result: 'list[KratosDT.Sensors.SensorSpecification]' = []
        for v in self.dict_of_sensor_specifications.values():
            result.extend(v)
        return result

    def PrintAnalysisStageProgressInformation(self):
        process_info: Kratos.ProcessInfo = self._GetSolver().GetComputingModelPart().ProcessInfo
        Kratos.Logger.PrintInfo(self._GetSimulationName(), f"STEP: {process_info[Kratos.STEP]} - {process_info[KratosDT.SENSOR_NAME]} - {process_info[KratosDT.SENSOR_ID]}")

    def UpdateSpecification(self, specification: KratosDT.Sensors.SensorSpecification) -> None:
        sensitivity_model_part: Kratos.ModelPart = self._GetSolver().GetSensitivityModelPart()
        sensitivity_variables: 'dict[Kratos.Globals.DataLocation, list[typing.Union[Kratos.DoubleVariable, Kratos.Array1DVariable3]]]' = self._GetSolver().GetSensitivtyVariables()
        for data_location, variables in sensitivity_variables.items():
            for variable in variables:
                # first read the variable
                expression = GetContainerExpression(sensitivity_model_part, data_location, variable)
                if isinstance(expression, Kratos.Expression.NodalExpression):
                    specification.AddNodalExpression(variable.Name(), expression)
                elif isinstance(expression, Kratos.Expression.ConditionExpression):
                    specification.AddConditionExpression(variable.Name(), expression)
                elif isinstance(expression, Kratos.Expression.ElementExpression):
                    specification.AddElementExpression(variable.Name(), expression)
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
    simulation = SensorSpecificationAnalysis(model, parameters)
    simulation.Run()
