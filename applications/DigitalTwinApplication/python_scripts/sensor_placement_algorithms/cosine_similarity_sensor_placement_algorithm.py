import typing
import scipy.cluster.hierarchy as sch
from pathlib import Path
import numpy as np
import json

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.DigitalTwinApplication as KratosDT
from KratosMultiphysics.DigitalTwinApplication.sensor_placement_algorithms.sensor_placement_algorithm import SensorPlacementAlgorithm
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import GetNormalizedSensorViews
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import GetCosineDistances
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import PrintSensorViewsListToCSV
from KratosMultiphysics.DigitalTwinApplication.utilities.data_utils import GetKratosValueToPythonValueConverter
from KratosMultiphysics.DigitalTwinApplication.utilities.data_utils import SensorViewUnionType
from KratosMultiphysics.DigitalTwinApplication.utilities.data_utils import SupportedVariableUnionType
from KratosMultiphysics.DigitalTwinApplication.utilities.data_utils import SupportedValueUnionType

class CosineSimilaritySensorPlacementAlgorithm(SensorPlacementAlgorithm):
    @classmethod
    def GetDefaultParameters(cls) -> Kratos.Parameters:
        return Kratos.Parameters("""{
            "type"                             : "cosine_similarity_sensor_placement",
            "output_to_vtu"                    : true,
            "output_to_csv"                    : true,
            "output_folder"                    : "sensor_placement/cluster_average",
            "clustering_method"                : "average",
            "sensor_coverage_percentage"       : 5.0,
            "max_clustering_iterations"        : 100,
            "best_sensor_identification_method": "highest_minimum",
            "filtering": {
                "filter_radius"             : 5.0,
                "filter_function_type"      : "linear",
                "fixed_model_part_name"     : "",
                "damping_function_type"     : "sigmoidal",
                "max_nodes_in_filter_radius": 1000
            }
        }""")

    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters) -> None:
        self.model = model
        self.parameters = parameters
        self.list_of_sensors:'list[KratosDT.Sensors.Sensor]' = []
        self.parameters.RecursivelyValidateAndAssignDefaults(self.GetDefaultParameters())
        self.is_vtu_output = self.parameters["output_to_vtu"].GetBool()
        self.is_csv_output = self.parameters["output_to_csv"].GetBool()
        self.best_sensor_identification_method = self.parameters["best_sensor_identification_method"].GetString()

    def Execute(self, list_of_sensors: 'list[KratosDT.Sensors.Sensor]') -> None:
        self.list_of_sensors = list_of_sensors

        if len(self.list_of_sensors) == 0:
            raise RuntimeError("No specifications are given.")

        self.SmoothenSensitivityFields()

        first_specification = self.list_of_sensors[0]

        for k in first_specification.GetNodalExpressionsMap().keys():
            self.ComputeSensorPlacement(f"nodal/{k}", [KratosDT.Sensors.NodalSensorView(sensor, k) for sensor in list_of_sensors])

        for k in first_specification.GetConditionExpressionsMap().keys():
            self.ComputeSensorPlacement(f"condition/{k}", [KratosDT.Sensors.ConditionSensorView(sensor, k) for sensor in list_of_sensors])

        for k in first_specification.GetElementExpressionsMap().keys():
            self.ComputeSensorPlacement(f"element/{k}", [KratosDT.Sensors.ElementSensorView(sensor, k) for sensor in list_of_sensors])

    def ComputeSensorPlacement(self, name: str, list_of_sensor_views: 'list[SensorViewUnionType]'):
        normalized_sensor_views = GetNormalizedSensorViews(list_of_sensor_views)

        if len(normalized_sensor_views) == 0:
            raise RuntimeError("No sensors with non-zero vectors found.")

        dummy_cexp = normalized_sensor_views[0].GetContainerExpression()
        vtu_output = Kratos.VtuOutput(dummy_cexp.GetModelPart())
        total_number_of_entities = len(dummy_cexp.GetContainer())

        # compute best sensor for each entity
        entity_best_sensors_list:'list[SensorViewUnionType]' = []
        list_of_entity_ids: 'list[int]' = []
        for entity_index, entity in enumerate(dummy_cexp.GetContainer()):
            entity_best_spec = max(normalized_sensor_views, key=lambda x: x.GetContainerExpression().Evaluate()[entity_index])
            entity_best_sensors_list.append(entity_best_spec)
            list_of_entity_ids.append(entity.Id)

        unique_normalized_list = list(set(entity_best_sensors_list))

        # get the entity domain sizes
        entity_domain_size_exp = unique_normalized_list[0].GetContainerExpression().Clone()
        Kratos.Expression.EntityDomainSizeExpressionIO.Read(entity_domain_size_exp)
        total_domain_size = KratosOA.ExpressionUtils.Sum(entity_domain_size_exp)
        entity_domain_size_np_exp = entity_domain_size_exp.Evaluate()
        maximum_allowed_domain_size_per_cluster = total_domain_size * self.parameters["sensor_coverage_percentage"].GetDouble() / 100.0

        sensor_cluster_dict = {1: unique_normalized_list}
        entity_index_cluster_dict = {1: range(total_number_of_entities)}

        output_path = Path(self.parameters["output_folder"].GetString()) / f"sensors/{name}"

        list_of_cluster_ids_cannot_be_divided: 'list[int]' = []

        clustering_iteration = 0
        list_of_clusters_to_break: 'list[int]' = [1]
        list_of_best_sensor_views: 'list[SensorViewUnionType]' = []
        while clustering_iteration < self.parameters["max_clustering_iterations"].GetInt() and set(list_of_cluster_ids_cannot_be_divided) != set(list_of_clusters_to_break) and len(list_of_clusters_to_break) > 0:
            Kratos.Logger.PrintInfo(self.__class__.__name__,f"Clustering iteration {clustering_iteration}")
            clustering_iteration += 1

            for cluster_id_to_break in list_of_clusters_to_break:
                if not cluster_id_to_break in list_of_cluster_ids_cannot_be_divided:
                    number_of_sensors_in_cluster = len(sensor_cluster_dict[cluster_id_to_break])
                    if number_of_sensors_in_cluster > 1:
                        Kratos.Logger.PrintInfo("", f"\tBreaking the cluster {cluster_id_to_break} having {number_of_sensors_in_cluster} sensor specifications...")
                        new_clusters = self.ClusterListOfSensorViews(2, sensor_cluster_dict[cluster_id_to_break])
                        sensor_cluster_dict[cluster_id_to_break] = new_clusters[1]
                        if 2 in new_clusters.keys():
                            sensor_cluster_dict[max(sensor_cluster_dict.keys()) + 1] = new_clusters[2]
                        else:
                            list_of_cluster_ids_cannot_be_divided.append(cluster_id_to_break)
                    else:
                        list_of_cluster_ids_cannot_be_divided.append(cluster_id_to_break)
            entity_index_cluster_dict = self.ComputeEntityClustering(entity_best_sensors_list, sensor_cluster_dict)

            list_of_clusters_to_break.clear()
            for cluster_id, entity_indices in entity_index_cluster_dict.items():
                # get the domain size of the cluster
                cluster_domain_size = np.sum(np.take(entity_domain_size_np_exp, entity_indices))

                if cluster_domain_size > maximum_allowed_domain_size_per_cluster:
                    list_of_clusters_to_break.append(cluster_id)

            list_of_best_sensor_views.clear()
            entity_cluster_data = np.array([-1] * total_number_of_entities, dtype=np.int32)

            for cluster_id, list_of_sensor_views in sensor_cluster_dict.items():
                for sensor_view in list_of_sensor_views:
                    sensor_view.GetSensor().SetValue(KratosDT.SENSOR_CLUSTER_ID, cluster_id)
                entity_indices = entity_index_cluster_dict[cluster_id]
                entity_cluster_data[entity_indices] = cluster_id

            # now try to find a representative sensor for each cluster
            if self.best_sensor_identification_method == "highest_minimum":
                # now look for the best sensor in each cluster which has the highest minimum for all the
                # entities which it relates to
                for cluster_id, list_of_sensors in sensor_cluster_dict.items():
                    entity_indices = entity_index_cluster_dict[cluster_id]
                    best_sensor = max(list_of_sensors, key=lambda x: np.min(np.take(x.GetContainerExpression().Evaluate(), entity_indices)))
                    list_of_best_sensor_views.append(best_sensor)
            elif self.best_sensor_identification_method == "most_similar":
                for cluster_id, list_of_sensors in sensor_cluster_dict.items():
                    average_array = list_of_sensors[0].GetContainerExpression().Evaluate()
                    for sensor in list_of_sensors[1:]:
                        average_array += sensor.GetContainerExpression().Evaluate()
                    best_sensor = max(list_of_sensors, key=lambda x: np.inner(average_array, x.GetContainerExpression().Evaluate()))
                    list_of_best_sensor_views.append(best_sensor)
            elif self.best_sensor_identification_method == "most_orthogonal" or self.best_sensor_identification_method == "most_distanced":
                # now look for the best sensor in each cluster is orthogonal to each other, starting with the most covered cluster in the sorted cluster ids
                cluster_values = {}
                overall_updating_exp = dummy_cexp.Clone()
                Kratos.Expression.CArrayExpressionIO.Read(overall_updating_exp, np.array([1.0] * total_number_of_entities))
                overall_updating_exp /= KratosOA.ExpressionUtils.NormL2(overall_updating_exp)
                for cluster_id, list_of_sensors in sensor_cluster_dict.items():
                    entity_indices = entity_index_cluster_dict[cluster_id]
                    best_sensor = max(list_of_sensors, key=lambda x: np.min(np.take(x.GetContainerExpression().Evaluate(), entity_indices)))
                    cluster_values[cluster_id] = KratosOA.ExpressionUtils.InnerProduct(best_sensor.GetContainerExpression(), overall_updating_exp)

                sorted_cluster_ids = sorted(sensor_cluster_dict.keys(), key=lambda x: cluster_values[x], reverse=True)

                best_cluster_id = sorted_cluster_ids[0]
                list_of_best_sensor_views.append(max(sensor_cluster_dict[best_cluster_id], key=lambda x: np.min(np.take(x.GetContainerExpression().Evaluate(), entity_index_cluster_dict[best_cluster_id]))))
                if self.best_sensor_identification_method == "most_orthogonal":
                    # first one is taken from the following
                    for cluster_id in sorted_cluster_ids[1:]:
                        list_of_sensors = sensor_cluster_dict[cluster_id]
                        list_of_best_sensor_views.append(max(list_of_sensors, key=lambda y: KratosOA.ExpressionUtils.InnerProduct(min(list_of_best_sensor_views, key=lambda x: KratosOA.ExpressionUtils.InnerProduct(x.GetContainerExpression(), y.GetContainerExpression())).GetContainerExpression(), y.GetContainerExpression())))
                elif self.best_sensor_identification_method == "most_distanced":
                    # first one is taken from the following
                    for cluster_id in sorted_cluster_ids[1:]:
                        list_of_sensors = sensor_cluster_dict[cluster_id]
                        list_of_best_sensor_views.append(max(list_of_sensors, key=lambda y: np.linalg.norm(min(list_of_best_sensor_views, key=lambda x: np.linalg.norm(x.GetSensor().GetLocation() - y.GetSensor().GetLocation())).GetSensor().GetLocation() - y.GetSensor().GetLocation())))
            else:
                raise RuntimeError(f"Unsupported best_sensor_identification_method = {self.best_sensor_identification_method}")

            vtu_output.ClearCellContainerExpressions()
            vtu_output.ClearNodalContainerExpressions()

            cexp = dummy_cexp.Clone()
            Kratos.Expression.CArrayExpressionIO.Read(cexp, entity_cluster_data)
            vtu_output.AddContainerExpression("entity_cluster_id", cexp.Clone())

            if len(list_of_best_sensor_views) > 0:
                heat_map = list_of_best_sensor_views[0].GetContainerExpression().Clone()
                for spec_view in list_of_best_sensor_views[1:]:
                    heat_map += spec_view.GetContainerExpression()
                heat_map /= KratosOA.ExpressionUtils.NormL2(heat_map)
                vtu_output.AddContainerExpression("heat_map", heat_map)

            PrintSensorViewsListToCSV(output_path / f"best_placement_{clustering_iteration:05d}.csv", [s_view.GetSensor() for s_view in list_of_best_sensor_views], ["type", "name", "location", "value", "SENSOR_CLUSTER_ID"])
            vtu_output.PrintOutput(str(output_path / f"heat_map_{clustering_iteration:05d}"))

        Kratos.Logger.PrintInfo(self.__class__.__name__, f"Found {len(list_of_best_sensor_views)} clusters with sensors.")

        if len(list_of_best_sensor_views) > 0:
            # write the json file for best sensor placement.
            list_of_vars: 'list[tuple[SupportedVariableUnionType, typing.Callable[[SupportedValueUnionType], typing.Union[bool, int, float, str, list[float]]]]]' = []
            for var_name in list_of_best_sensor_views[0].GetSensor().GetDataVariableNames():
                var: SupportedVariableUnionType = Kratos.KratosGlobals.GetVariable(var_name)
                list_of_vars.append((var, GetKratosValueToPythonValueConverter(list_of_best_sensor_views[0].GetSensor().GetValue(var))))
            list_of_vars.append((KratosDT.SENSOR_ENTITY_IDS, GetKratosValueToPythonValueConverter(Kratos.Vector())))

            # now write the sensors with data
            json_sensors = {"list_of_sensors": []}
            with open(str(output_path / "best_sensor_data.json"), "w") as file_output:
                for sensor_view in list_of_best_sensor_views:
                    sensor = sensor_view.GetSensor()
                    json_params = json.loads(sensor.GetSensorParameters().WriteJsonString())
                    json_params["variable_data"] = {}

                    # add cluster ids
                    list_of_entity_indices = entity_index_cluster_dict[sensor.GetValue(KratosDT.SENSOR_CLUSTER_ID)]
                    sensor.SetValue(KratosDT.SENSOR_ENTITY_IDS, Kratos.Vector([float(list_of_entity_ids[i]) for i in list_of_entity_indices]))
                    for var, func in list_of_vars:
                        json_params["variable_data"][var.Name()] = func(sensor.GetValue(var))
                    json_sensors["list_of_sensors"].append(json_params)
                file_output.write(json.dumps(json_sensors, indent=4))

        if len(list_of_clusters_to_break) > 0:
            Kratos.Logger.PrintWarning(self.__class__.__name__, f"The max_clustering_iterations and/or unbreakable clusters limit reached without breaking {len(list_of_clusters_to_break)} clusters.\n\t Following clusters does not satisfy the prescribed coverage area:")
            for cluster_id_to_break in list_of_clusters_to_break:
                entity_indices = entity_index_cluster_dict[cluster_id_to_break]
                cluster_domain_size = np.sum(np.take(entity_domain_size_np_exp, entity_indices))
                Kratos.Logger.PrintInfo("", f"\t\t Cluster {cluster_id_to_break} - domain size = {cluster_domain_size * 100.0 / total_domain_size:0.3f} %")

        PrintSensorViewsListToCSV(output_path / f"clusters.csv", [s_view.GetSensor() for s_view in unique_normalized_list], ["type", "name", "location", "value", "SENSOR_CLUSTER_ID"])

    def ClusterListOfSensorViews(self, number_of_cluster: int, list_of_sensor_views: 'list[SensorViewUnionType]') -> 'dict[int, list[SensorViewUnionType]]':
        sensor_cosine_distances = GetCosineDistances(list_of_sensor_views)
        sensor_cosine_linkage = sch.linkage(sensor_cosine_distances, method=self.parameters["clustering_method"].GetString())
        sensor_cosine_clusters = sch.fcluster(sensor_cosine_linkage, number_of_cluster, 'maxclust')

        cluster_dict: 'dict[int, list[SensorViewUnionType]]' = {}
        for i, cluster_id in enumerate(sensor_cosine_clusters):
            if cluster_id not in cluster_dict.keys():
                cluster_dict[cluster_id] = []
            cluster_dict[cluster_id].append(list_of_sensor_views[i])
        return cluster_dict

    def ComputeEntityClustering(self, entity_best_sensors_list: 'list[SensorViewUnionType]', cluster_dict: 'dict[int, list[SensorViewUnionType]]'):
        entity_clusters:'dict[int, list[int]]' = {}
        for cluster_id, sensors_list in cluster_dict.items():
            entity_clusters[cluster_id] = []
            for entity_index, entity_best_spec in enumerate(entity_best_sensors_list):
                if entity_best_spec in sensors_list:
                    entity_clusters[cluster_id].append(entity_index)
        return entity_clusters

    def SmoothenSensitivityFields(self) -> None:
        mp_nodal, nodal_filter = self.__GetFilter(KratosOA.NodalExplicitFilter)
        mp_condition, condition_filter = self.__GetFilter(KratosOA.ConditionExplicitFilter)
        mp_element, element_filter = self.__GetFilter(KratosOA.ElementExplicitFilter)

        output_path = Path(self.parameters["output_folder"].GetString())

        model_part = mp_nodal
        if model_part is None:
            model_part = mp_condition
        if model_part is None:
            model_part = mp_element

        if mp_nodal not in [None, model_part] or mp_condition not in [None, model_part] or mp_element not in [None, model_part]:
            raise RuntimeError(f"Found mismatching model parts.")

        if self.is_vtu_output:
            vtu_output_path = output_path / "sensitivities"
            vtu_output_path.mkdir(parents=True, exist_ok=True)
            vtu_output = Kratos.VtuOutput(model_part)

        for sensor in self.list_of_sensors:
            if self.is_vtu_output:
                vtu_output.ClearCellContainerExpressions()
                vtu_output.ClearNodalContainerExpressions()
            for k, v in sensor.GetNodalExpressionsMap().items():
                filtered_field = nodal_filter.FilterIntegratedField(v.Abs())
                if self.is_vtu_output:
                    vtu_output.AddContainerExpression(f"{k}_raw", v.Clone())
                    vtu_output.AddContainerExpression(f"{k}_filtered", filtered_field.Clone())

                v.SetExpression(filtered_field.GetExpression())

            for k, v in sensor.GetConditionExpressionsMap().items():
                filtered_field = condition_filter.FilterIntegratedField(v.Abs())
                if self.is_vtu_output:
                    vtu_output.AddContainerExpression(f"{k}_raw", v.Clone())
                    vtu_output.AddContainerExpression(f"{k}_filtered", filtered_field.Clone())

                v.SetExpression(filtered_field.GetExpression())

            for k, v in sensor.GetElementExpressionsMap().items():
                filtered_field = element_filter.FilterIntegratedField(v.Abs())
                if self.is_vtu_output:
                    vtu_output.AddContainerExpression(f"{k}_raw", v.Clone())
                    vtu_output.AddContainerExpression(f"{k}_filtered", filtered_field.Clone())

                v.SetExpression(filtered_field.GetExpression())

            if self.is_vtu_output:
                vtu_output.PrintOutput(str(vtu_output_path / f"{sensor.GetName()}"))

    def __GetFilter(self, filter_type: 'typing.Union[typing.Type[KratosOA.NodalExplicitFilter], typing.Type[KratosOA.ConditionExplicitFilter], typing.Type[KratosOA.ElementExplicitFilter]]') -> 'tuple[Kratos.ModelPart, typing.Union[KratosOA.NodalExplicitFilter, KratosOA.ConditionExplicitFilter, KratosOA.ElementExplicitFilter]]':
        filter_settings = self.parameters["filtering"]
        filter_settings.ValidateAndAssignDefaults(self.GetDefaultParameters()["filtering"])

        fixed_model_part_name = filter_settings["fixed_model_part_name"].GetString()
        filter_function_type = filter_settings["filter_function_type"].GetString()
        damping_function_type = filter_settings["damping_function_type"].GetString()
        filter_radius = filter_settings["filter_radius"].GetDouble()
        max_filtering_nodes = filter_settings["max_nodes_in_filter_radius"].GetInt()

        # first create the nodal filter
        vm_filter: filter_type = None
        model_part: Kratos.ModelPart = None
        for sensor in self.list_of_sensors:
            if filter_type == KratosOA.NodalExplicitFilter:
                expressions_list = sensor.GetNodalExpressionsMap().values()
                filter_radius_exp_type = Kratos.Expression.NodalExpression
            elif filter_type == KratosOA.ConditionExplicitFilter:
                expressions_list = sensor.GetConditionExpressionsMap().values()
                filter_radius_exp_type = Kratos.Expression.ConditionExpression
            elif filter_type == KratosOA.ElementExplicitFilter:
                expressions_list = sensor.GetElementExpressionsMap().values()
                filter_radius_exp_type = Kratos.Expression.ElementExpression
            for v in expressions_list:
                if vm_filter is None:
                    model_part = v.GetModelPart()
                    if fixed_model_part_name == "":
                        vm_filter = filter_type(model_part, filter_function_type, max_filtering_nodes)
                    else:
                        vm_filter = filter_type(model_part, self.model[fixed_model_part_name], filter_function_type, damping_function_type, max_filtering_nodes)
                    filter_radius_exp = filter_radius_exp_type(model_part)
                    Kratos.Expression.LiteralExpressionIO.SetData(filter_radius_exp, filter_radius)
                    vm_filter.SetFilterRadius(filter_radius_exp)

                    Kratos.Logger.PrintInfo(self.__class__.__name__, f"Created a {filter_type.__name__} filter for {model_part.FullName()}.")
                    break

            if vm_filter is not None:
                break
        return model_part, vm_filter





