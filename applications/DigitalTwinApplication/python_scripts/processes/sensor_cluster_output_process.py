from pathlib import Path
import typing
import KratosMultiphysics as Kratos
import KratosMultiphysics.DigitalTwinApplication as KratosDT
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import PrintSensorListToJson

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
    return SensorClusterOutputProcess(model, parameters["settings"], optimization_problem)

class SensorClusterOutputProcess(Kratos.OutputProcess):
    def __init__(self, model: Kratos.Model, settings: Kratos.Parameters, optimization_problem: OptimizationProblem):
        Kratos.OutputProcess.__init__(self)

        default_settings = Kratos.Parameters("""{
            "model_part_name"        : "PLEASE_PROVIDE_A_MODEL_PART_NAME",
            "output_file_name_prefix": "",
            "mask_expression_name"   : "YOUNG_MODULUS_SENSITIVITY_mask",
            "sensor_status_threshold": 0.6,
            "output_every_step"      : false
        }""")

        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[settings["model_part_name"].GetString()]
        self.output_file_name_prefix = settings["output_file_name_prefix"].GetString()
        self.sensor_status_threshold = settings["sensor_status_threshold"].GetDouble()
        self.mask_expression_name = settings["mask_expression_name"].GetString()
        self.output_every_step = settings["output_every_step"].GetBool()

        self.optimization_problem = optimization_problem
        self.mask_model_part: 'typing.Optional[Kratos.ModelPart]' = None
        self.cluster_data: 'dict[int, tuple[list[int], Kratos.Expression.ElementExpression]]' = {}

    def ExecuteInitialize(self) -> None:
        sensor_data = ComponentDataView("sensor", self.optimization_problem)
        list_of_sensors: 'list[KratosDT.Sensors.Sensor]' = [sensor_data.GetUnBufferedData()["list_of_sensors"][i] for i, node in enumerate(self.model_part.Nodes)]
        if len(list_of_sensors) == 0:
            raise RuntimeError("No sensors found.")

        self.mask_model_part = list_of_sensors[0].GetElementExpression(self.mask_expression_name).GetModelPart()

    def IsOutputStep(self) -> bool:
        return self.output_every_step

    def PrintOutput(self) -> None:
        file_name_prefix = self.output_file_name_prefix.replace("<step>", str(self.optimization_problem.GetStep()))
        Path(file_name_prefix).parent.mkdir(exist_ok=True, parents=True)
        self.__Output(file_name_prefix)

    def ExecuteFinalize(self) -> None:
        file_name_prefix = self.output_file_name_prefix.replace("<step>", "final")
        Path(file_name_prefix).parent.mkdir(exist_ok=True, parents=True)
        self.__Output(file_name_prefix)

    def __Output(self, output_file_name_prefix: str) -> None:
        sensor_data = ComponentDataView("sensor", self.optimization_problem)
        list_of_sensors: 'list[KratosDT.Sensors.Sensor]' = [sensor_data.GetUnBufferedData()["list_of_sensors"][i] for i, node in enumerate(self.model_part.Nodes) if node.GetValue(KratosDT.SENSOR_STATUS) > self.sensor_status_threshold]
        cluster_data_list = KratosDT.MaskUtils.ClusterMasks([sensor.GetElementExpression(self.mask_expression_name) for sensor in list_of_sensors])

        # add or modify existing clusters
        list_of_unavailable_clusters = list(self.cluster_data.keys())
        for cluster_mask_indices, cluster_exp in cluster_data_list:
            sensor_ids = [list_of_sensors[cluster_mask_index].GetValue(KratosDT.SENSOR_ID) for cluster_mask_index in cluster_mask_indices]

            found_cluster = False
            for cluster_id, (existing_cluster_sensor_ids, exp) in self.cluster_data.items():
                if sensor_ids == existing_cluster_sensor_ids:
                    exp.SetExpression(cluster_exp.GetExpression())
                    found_cluster = True
                    del list_of_unavailable_clusters[list_of_unavailable_clusters.index(cluster_id)]
                    break

            if not found_cluster:
                if len(self.cluster_data) > 0:
                    new_cluster_id = max(self.cluster_data.keys()) + 1
                else:
                    new_cluster_id = 1
                self.cluster_data[new_cluster_id] = (sensor_ids, cluster_exp)

        # remove unavailable clusters
        for cluster_id in list_of_unavailable_clusters:
            del self.cluster_data[cluster_id]

        domain_size_exp = Kratos.Expression.ElementExpression(self.mask_model_part)
        Kratos.Expression.DomainSizeExpressionIO.Read(domain_size_exp)
        total_size = Kratos.Expression.Utils.Sum(domain_size_exp)

        cluster_ids: 'list[int]' = sorted(self.cluster_data.keys())

        with open(f"{output_file_name_prefix}.csv", "w") as file_output:
            file_output.write( "# Cluster output:\n")
            file_output.write(f"#          Number of sensors: {len(list_of_sensors)}\n")
            file_output.write(f"#          Number of clusters: {len(self.cluster_data.keys())}\n")
            file_output.write(f"# Cluster Id, Cluster size, Cluster size ratio\n")
            for cluster_id in cluster_ids:
                cluster_size = Kratos.Expression.Utils.InnerProduct(domain_size_exp, self.cluster_data[cluster_id][1])
                file_output.write(f"{cluster_id:12d}, {cluster_size:12.6e}, {cluster_size/total_size:18.6e}\n")
            file_output.write("# End of cluster data.")

        overall_cluster_exp = Kratos.Expression.ElementExpression(self.mask_model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(overall_cluster_exp, 0)
        for cluster_id, (_, cluster_exp) in self.cluster_data.items():
            overall_cluster_exp = Kratos.Expression.Utils.Collapse(overall_cluster_exp + cluster_exp * cluster_id)
        vtu_output = Kratos.VtuOutput(self.mask_model_part)
        vtu_output.AddContainerExpression("clusters", overall_cluster_exp)
        vtu_output.PrintOutput(output_file_name_prefix)






