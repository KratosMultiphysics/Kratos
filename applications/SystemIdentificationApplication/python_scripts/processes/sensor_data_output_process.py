import json
import pathlib
import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import GetSensors
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes

def Factory(_: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> Kratos.Process:
    if not parameters.Has("settings"):
        raise RuntimeError(f"SensorDataOutputProcess instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return SensorDataOutputProcess(parameters["settings"], optimization_problem)


class SensorDataOutputProcess(Kratos.OutputProcess):
    def __init__(self, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        Kratos.OutputProcess.__init__(self)

        default_parameters = Kratos.Parameters("""{
            "sensor_group_name"          : "",
            "sensor_mask_name"           : "",
            "sensor_activation_threshold": 0.9,
            "output_file_name"           : "",
            "output_interval"            : 1,
            "compute_clustering"         : true
        }""")

        parameters.ValidateAndAssignDefaults(default_parameters)
        self.sensor_group_name = parameters["sensor_group_name"].GetString()
        self.output_file_name = parameters["output_file_name"].GetString()
        self.output_interval = parameters["output_interval"].GetInt()
        self.compute_clustering = parameters["compute_clustering"].GetBool()
        self.sensor_mask_name = parameters["sensor_mask_name"].GetString()
        self.sensor_activation_threshold = parameters["sensor_activation_threshold"].GetDouble()

        self.optimization_problem = optimization_problem

        self.last_output_step = self.optimization_problem.GetStep()

    def ExecuteInitialize(self):
        if self.compute_clustering:
            component_data_view = ComponentDataView(self.sensor_group_name, self.optimization_problem)
            list_of_sensors = GetSensors(component_data_view)
            overall_cluster = list_of_sensors[0].GetContainerExpression(self.sensor_mask_name) * 0.0
            component_data_view.GetUnBufferedData().SetValue("clustering", overall_cluster.Clone(), overwrite=True)

    def IsOutputStep(self):
        if self.optimization_problem.GetStep() - self.last_output_step >= self.output_interval:
            self.last_output_step = self.optimization_problem.GetStep()
            self.PrintOutput()

    def PrintOutput(self) -> None:
        component_data_view = ComponentDataView(self.sensor_group_name, self.optimization_problem)
        list_of_sensors = GetSensors(component_data_view)

        current_file_name = self.output_file_name.replace("<STEP>", str(self.optimization_problem.GetStep()))
        pathlib.Path(current_file_name).parent.mkdir(exist_ok=True, parents=True)

        json_sensors = {"list_of_sensors": []}
        with open(current_file_name, "w") as file_output:
            for sensor in list_of_sensors:
                if sensor.GetNode().GetValue(KratosSI.SENSOR_STATUS) > self.sensor_activation_threshold:
                    json_params = json.loads(sensor.GetSensorParameters().WriteJsonString())
                    json_sensors["list_of_sensors"].append(json_params)
            file_output.write(json.dumps(json_sensors, indent=4))

        if self.compute_clustering:
            overall_cluster = list_of_sensors[0].GetContainerExpression(self.sensor_mask_name) * 0.0

            list_of_masks: 'list[ContainerExpressionTypes]' = []
            for sensor in list_of_sensors:
                if sensor.GetNode().GetValue(KratosSI.SENSOR_STATUS) > self.sensor_activation_threshold:
                    list_of_masks.append(sensor.GetContainerExpression(self.sensor_mask_name))

            cluster_info = KratosSI.MaskUtils.ClusterMasks(list_of_masks)

            cluster_id_exp_list: 'list[tuple[int, ContainerExpressionTypes]]' = []
            cluster_id = 1
            for cluster_indices_list, cluster_exp in cluster_info:
                if cluster_indices_list != []:
                    cluster_id_exp_list.append((cluster_id, cluster_exp.Clone()))
                    cluster_id += 1
                else:
                    cluster_id_exp_list.append((0, cluster_exp.Clone()))

            for cluster_id, cluster_exp in cluster_id_exp_list:
                overall_cluster += cluster_exp * cluster_id

            component_data_view.GetUnBufferedData().SetValue("clustering", overall_cluster.Clone(), overwrite=True)

    def ExecuteFinalize(self):
        self.PrintOutput()

