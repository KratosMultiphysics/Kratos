import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosDT
from KratosMultiphysics.OptimizationApplication.model_part_controllers.model_part_controller import ModelPartController
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import GetSensors
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ModelPartController:
    if not parameters.Has("settings"):
        raise RuntimeError(f"SensorDataModelPartController instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return SensorModelPartController(model, parameters["settings"], optimization_problem)

class SensorModelPartController(ModelPartController):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        default_settings = Kratos.Parameters("""{
            "sensor_model_part_name": "",
            "target_model_part_name": "",
            "list_of_sensors"       : []
        }""")

        self.model = model
        self.optimization_problem = optimization_problem
        self.parameters = parameters
        self.parameters.ValidateAndAssignDefaults(default_settings)

        self.target_model_part_name = parameters["target_model_part_name"].GetString()

        sensor_model_part_name = parameters["sensor_model_part_name"].GetString()
        if model.HasModelPart(sensor_model_part_name):
            raise RuntimeError(f"The sensor model part name \"{sensor_model_part_name}\" already exists.")
        self.sensor_model_part = model.CreateModelPart(sensor_model_part_name)

    def ImportModelPart(self) -> None:
        self.targe_model_part = self.model[self.target_model_part_name]
        self.list_of_sensors = GetSensors(self.targe_model_part, self.parameters["list_of_sensors"].values())

        # add the list of sensors to optimization problem
        sensors = ComponentDataView("sensors", self.optimization_problem)
        sensors.SetDataBuffer(1)
        sensors.GetUnBufferedData().SetValue("list_of_sensors", self.list_of_sensors)

        for i, sensor in enumerate(self.list_of_sensors):
            location = sensor.GetLocation()
            node: Kratos.Node = self.sensor_model_part.CreateNewNode(i + 1, location[0], location[1], location[2])
            for var_name in sensor.GetDataVariableNames():
                var = Kratos.KratosGlobals.GetVariable(var_name)
                node.SetValue(var, sensor.GetValue(var))
            # sensor.SetValue(KratosDT.SENSOR_ID, i + 1)

    def GetModelPart(self) -> Kratos.ModelPart:
        return self.sensor_model_part