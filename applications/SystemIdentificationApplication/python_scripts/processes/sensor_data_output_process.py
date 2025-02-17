import json
import pathlib
import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import GetSensors

def Factory(_: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> Kratos.Process:
    if not parameters.Has("settings"):
        raise RuntimeError(f"SensorDataOutputProcess instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return SensorDataOutputProcess(parameters["settings"], optimization_problem)


class SensorDataOutputProcess(Kratos.OutputProcess):
    def __init__(self, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        Kratos.OutputProcess.__init__(self)

        default_parameters = Kratos.Parameters("""{
            "sensor_group_name"          : "",
            "sensor_activation_threshold": 0.9,
            "output_file_name"           : "",
            "output_interval"            : 1
        }""")

        parameters.ValidateAndAssignDefaults(default_parameters)
        self.sensor_group_name = parameters["sensor_group_name"].GetString()
        self.output_file_name = parameters["output_file_name"].GetString()
        self.output_interval = parameters["output_interval"].GetInt()
        self.sensor_activation_threshold = parameters["sensor_activation_threshold"].GetDouble()

        self.optimization_problem = optimization_problem

        self.last_output_step = self.optimization_problem.GetStep()

    def IsOutputStep(self):
        if self.optimization_problem.GetStep() - self.last_output_step >= self.output_interval:
            self.last_output_step = self.optimization_problem.GetStep()
            self.PrintOutput()

    def PrintOutput(self) -> None:
        list_of_sensors = GetSensors(ComponentDataView(self.sensor_group_name, self.optimization_problem))

        current_file_name = self.output_file_name.replace("<STEP>", str(self.optimization_problem.GetStep()))
        pathlib.Path(current_file_name).parent.mkdir(exist_ok=True, parents=True)

        json_sensors = {"list_of_sensors": []}
        with open(current_file_name, "w") as file_output:
            for sensor in list_of_sensors:
                if sensor.GetNode().GetValue(KratosSI.SENSOR_STATUS) > self.sensor_activation_threshold:
                    json_params = json.loads(sensor.GetSensorParameters().WriteJsonString())
                    json_sensors["list_of_sensors"].append(json_params)
            file_output.write(json.dumps(json_sensors, indent=4))

    def ExecuteFinalize(self):
        self.PrintOutput()

