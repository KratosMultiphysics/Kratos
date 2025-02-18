import h5py # TODO: Remove once HDF5Application is properly setup to function with Expressions
import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
from KratosMultiphysics.OptimizationApplication.model_part_controllers.model_part_controller import ModelPartController
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import GetSensors
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import AddMaskStatusController
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ModelPartController:
    if not parameters.Has("settings"):
        raise RuntimeError(f"SensorDataModelPartController instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return SensorModelPartController(model, parameters["settings"], optimization_problem)

class SensorModelPartController(ModelPartController):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        default_settings = Kratos.Parameters("""{
            "sensor_group_name": "",
            "sensor_mask_name" : "",
            "h5_file_name"     : "",
            "data_field_name"  : "",
            "echo_level"       : 0,
            "kd_tree_settings" : {
                "use_kd_tree"  : true,
                "leaf_max_size": 100
            }
        }""")

        self.model = model
        self.optimization_problem = optimization_problem
        parameters.ValidateAndAssignDefaults(default_settings)

        self.sensor_group_name = parameters["sensor_group_name"].GetString()
        self.sensor_mask_name = parameters["sensor_mask_name"].GetString()
        self.h5_file_name = parameters["h5_file_name"].GetString()
        self.data_field_name = parameters["data_field_name"].GetString()
        self.echo_level = parameters["echo_level"].GetInt()

        parameters["kd_tree_settings"].ValidateAndAssignDefaults(default_settings["kd_tree_settings"])

        self.use_kd_tree = parameters["kd_tree_settings"]["use_kd_tree"].GetBool()
        self.leaf_max_size = parameters["kd_tree_settings"]["leaf_max_size"].GetInt()

    def ImportModelPart(self) -> None:
        sensor_group_data = ComponentDataView(self.sensor_group_name, self.optimization_problem)
        list_of_sensors = GetSensors(sensor_group_data)

        expression_name = self.data_field_name.split("/")[-1]
        list_of_masks = []
        with h5py.File(self.h5_file_name, "r") as h5_file:
            for sensor in list_of_sensors:
                sensor_id = sensor.GetNode().Id
                sensor_name = sensor.GetName()
                current_sensor_data_field_name = self.data_field_name.replace("<SENSOR_NAME>", sensor_name)
                current_sensor_data_field_name = current_sensor_data_field_name.replace("<SENSOR_ID>", str(sensor_id))

                dataset = h5_file.get(current_sensor_data_field_name)
                if not isinstance(dataset, h5py.Dataset):
                    raise RuntimeError(f"The sensor data not found or not a dataset at \"{current_sensor_data_field_name}\" [ sensor name = \"{sensor_name}\", sensor id = \"{sensor_id}\" ]. Found data = \n{dataset}")

                container_type = dataset.attrs["__container_type"]
                model_part = self.model[dataset.attrs["__model_part_name"]]
                if container_type == "NodalExpression":
                    expression = Kratos.Expression.NodalExpression(model_part)
                elif container_type == "ConditionExpression":
                    expression = Kratos.Expression.ConditionExpression(model_part)
                elif container_type == "ElementExpression":
                    expression = Kratos.Expression.ElementExpression(model_part)
                else:
                    raise RuntimeError(f"Unsupported container type = \"{container_type}\" requested for dataset at \"{current_sensor_data_field_name}\".")

                Kratos.Expression.CArrayExpressionIO.Read(expression, dataset[:])
                sensor.AddContainerExpression(self.sensor_mask_name, expression)
                list_of_masks.append(expression.Clone())

        # now create the mask
        self.sensor_mask_status = KratosSI.SensorMaskStatus(self.model[self.sensor_group_name], list_of_masks, self.echo_level)
        AddMaskStatusController(sensor_group_data, self.sensor_mask_name, self.sensor_mask_status)

        if self.use_kd_tree:
            self.sensor_mask_status_kd_tree = KratosSI.SensorMaskStatusKDTree(self.sensor_mask_status, self.leaf_max_size, self.echo_level)
            AddMaskStatusController(sensor_group_data, self.sensor_mask_name, self.sensor_mask_status_kd_tree)

    def GetModelPart(self) -> Kratos.ModelPart:
        return self.model[self.sensor_group_name]


