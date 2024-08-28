import h5py

import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.model_part_controllers.model_part_controller import ModelPartController
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import GetSensors

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ModelPartController:
    if not parameters.Has("settings"):
        raise RuntimeError(f"SensorDataInputController instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return SensorDataInputController(model, parameters["settings"], optimization_problem)

class SensorDataInputController(Kratos.Process):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        default_settings = Kratos.Parameters("""{
            "h5_file_name"   : "",
            "data_field_name": ""
        }""")

        self.model = model
        self.optimization_problem = optimization_problem
        parameters.ValidateAndAssignDefaults(default_settings)

        self.h5_file_name = parameters["h5_file_name"].GetString()
        self.data_field_name = parameters["data_field_name"].GetString()

    def ImportModelPart(self) -> None:
        list_of_sensors = GetSensors(self.optimization_problem)

        with h5py.File(self.h5_file_name, "r") as h5_file:
            for sensor in list_of_sensors:
                sensor_id = sensor.GetNode().Id
                sensor_name = sensor.GetName()
                current_sensor_data_field_name = self.data_field_name.replace("<SENSOR_NAME>", sensor_name)
                current_sensor_data_field_name = current_sensor_data_field_name.replace("<SENSOR_ID>", sensor_id)

                dataset = h5_file.get(current_sensor_data_field_name)
                if isinstance(dataset, h5py.Dataset):
                    raise RuntimeError(f"The sensor data not found or not a dataset at \"{current_sensor_data_field_name}\" [ sensor name = \"{sensor_name}\", sensor id = \"{sensor_id}\" ].")

                expression_name = current_sensor_data_field_name.split("/")[-1]
                container_type = dataset.attrs["__container_type"]
                model_part = self.model[dataset.attrs["__model_part_name"]]
                if container_type == "NODES":
                    expression = Kratos.Expression.NodalExpression(model_part)
                elif container_type == "CONDITIONS":
                    expression = Kratos.Expression.ConditionExpression(model_part)
                elif container_type == "ELEMENTS":
                    expression = Kratos.Expression.ElementExpression(model_part)
                else:
                    raise RuntimeError(f"Unsupported container type = \"{container_type}\" requested for dataset at \"{current_sensor_data_field_name}\".")

                Kratos.Expression.CArrayExpressionIO.Read(expression, dataset)
                sensor.AddContainerExpression(expression_name, expression)

    def GetModelPart(self) -> Kratos.ModelPart:
        for sensor in GetSensors(self.optimization_problem):
            return list(sensor.GetContainerExpressionsMap().values())[0].GetModelPart()