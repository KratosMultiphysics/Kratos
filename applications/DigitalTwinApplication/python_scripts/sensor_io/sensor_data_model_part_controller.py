import KratosMultiphysics as Kratos
import KratosMultiphysics.DigitalTwinApplication as KratosDT
from KratosMultiphysics.OptimizationApplication.model_part_controllers.model_part_controller import ModelPartController
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_io import OpenSensorFile
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import GetSensors
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ModelPartController:
    if not parameters.Has("settings"):
        raise RuntimeError(f"SensorDataModelPartController instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return SensorDataModelPartController(model, parameters["settings"], optimization_problem)

class SensorDataModelPartController(ModelPartController):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        default_settings = Kratos.Parameters("""{
            "target_model_part_name": "",
            "sensor_data_file_name" : "",
            "sensor_data_prefix"    : "",
            "sensor_model_part_name": "",
            "list_of_sensors"       : []
        }""")

        self.model = model
        self.optimization_problem = optimization_problem
        self.parameters = parameters
        self.parameters.ValidateAndAssignDefaults(default_settings)

        self.target_model_part_name = parameters["target_model_part_name"].GetString()
        self.sensor_data_file_name = parameters["sensor_data_file_name"].GetString()
        self.prefix = parameters["sensor_data_prefix"].GetString()

        sensor_model_part_name = parameters["sensor_model_part_name"].GetString()
        if model.HasModelPart(sensor_model_part_name):
            raise RuntimeError(f"The sensor model part name \"{sensor_model_part_name}\" already exists.")
        self.sensor_model_part = model.CreateModelPart(sensor_model_part_name)

    def ImportModelPart(self) -> None:
        self.targe_model_part = self.model[self.target_model_part_name]
        self.list_of_sensors = GetSensors(self.targe_model_part, self.parameters["list_of_sensors"].values())

        # add the list of sensors to optimization problem
        sensors = ComponentDataView("sensor", self.optimization_problem)
        sensors.GetUnBufferedData().SetValue("list_of_sensors", self.list_of_sensors)

        # create nodes for every sensor
        with OpenSensorFile(self.targe_model_part, self.sensor_data_file_name, self.prefix, "r") as sensor_io:
            for i, sensor in enumerate(self.list_of_sensors):
                sensor_io.Read(sensor)
                location = sensor.GetLocation()
                node: Kratos.Node = self.sensor_model_part.CreateNewNode(i + 1, location[0], location[1], location[2])
                # set data
                for var_name in sensor.GetDataVariableNames():
                    var = Kratos.KratosGlobals.GetVariable(var_name)
                    node.SetValue(var, sensor.GetValue(var))
                sensor.SetValue(KratosDT.SENSOR_ID, i + 1)

    def GetModelPart(self) -> Kratos.ModelPart:
        return self.sensor_model_part