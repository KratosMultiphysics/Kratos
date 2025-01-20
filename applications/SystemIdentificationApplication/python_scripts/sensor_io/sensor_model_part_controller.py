import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.model_part_controllers.model_part_controller import ModelPartController
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import CreateSensors
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import SetSensors
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ModelPartController:
    if not parameters.Has("settings"):
        raise RuntimeError(f"SensorModelPartController instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return SensorModelPartController(model, parameters["settings"], optimization_problem)

class SensorModelPartController(ModelPartController):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        default_settings = Kratos.Parameters("""{
            "sensor_group_name"     : "",
            "domain_model_part_name": "",
            "list_of_sensors"       : []
        }""")

        self.model = model
        self.optimization_problem = optimization_problem
        self.parameters = parameters
        self.parameters.ValidateAndAssignDefaults(default_settings)

        self.domain_model_part_name = parameters["domain_model_part_name"].GetString()
        self.h5_file_name = parameters["h5_file_name"].GetString()
        self.data_field_name = parameters["data_field_name"].GetString()

        self.sensor_group_name = parameters["sensor_group_name"].GetString()

        if model.HasModelPart(self.sensor_group_name):
            raise RuntimeError(f"The sensor model part name \"{self.sensor_group_name}\" already exists. Please provide a different \"sensor_group_name\".")
        self.sensor_model_part = model.CreateModelPart(self.sensor_group_name)

    def ImportModelPart(self) -> None:
        self.domain_model_part = self.model[self.domain_model_part_name]
        self.list_of_sensors = CreateSensors(self.sensor_model_part, self.domain_model_part, self.parameters["list_of_sensors"].values())

        # add the list of sensors to optimization problem
        sensor_group_data = ComponentDataView(self.sensor_group_name, self.optimization_problem)
        sensor_group_data.SetDataBuffer(1)

        SetSensors(sensor_group_data, self.list_of_sensors)

    def GetModelPart(self) -> Kratos.ModelPart:
        return self.sensor_model_part


