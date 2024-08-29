import h5py

import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.model_part_controllers.model_part_controller import ModelPartController
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import GetSensors

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ModelPartController:
    if not parameters.Has("settings"):
        raise RuntimeError(f"SensorDataOutputController instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return SensorDataOutputController(model, parameters["settings"], optimization_problem)

class SensorDataOutputController(Kratos.OutputProcess):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        Kratos.OutputProcess.__init__(self)

        default_settings = Kratos.Parameters("""{
            "h5_file_name"     : "",
            "data_field_prefix": ""
        }""")

        self.model = model
        self.optimization_problem = optimization_problem
        parameters.ValidateAndAssignDefaults(default_settings)

        self.h5_file_name = parameters["h5_file_name"].GetString()
        self.data_field_prefix = parameters["data_field_prefix"].GetString()

    def IsOutputStep(self) -> bool:
        return True

    def PrintOutput(self) -> None:
        list_of_sensors = GetSensors(self.optimization_problem)

        with h5py.File(self.h5_file_name, "a") as h5_file:
            for sensor in list_of_sensors:
                sensor_id = sensor.GetNode().Id
                sensor_name = sensor.GetName()
                current_sensor_group = self.data_field_prefix.replace("<SENSOR_NAME>", sensor_name)
                current_sensor_group = current_sensor_group.replace("<SENSOR_ID>", str(sensor_id))

                for expression_name, expression in sensor.GetContainerExpressionsMap().items():
                    h5_file[f"/{current_sensor_group}/{expression_name}"] = expression.Evaluate()
                    h5_file[f"/{current_sensor_group}/{expression_name}"].attrs["__model_part_name"] = expression.GetModelPart().FullName()
                    h5_file[f"/{current_sensor_group}/{expression_name}"].attrs["__container_type"] = expression.__class__.__name__
