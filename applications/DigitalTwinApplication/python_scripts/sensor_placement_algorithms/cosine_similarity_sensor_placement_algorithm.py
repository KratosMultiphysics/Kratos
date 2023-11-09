import typing
import scipy.cluster.hierarchy as sch
from pathlib import Path
from functools import reduce
import numpy as np

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.DigitalTwinApplication as KratosDT
from KratosMultiphysics.DigitalTwinApplication.sensor_placement_algorithms.sensor_placement_algorithm import SensorPlacementAlgorithm
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_specification_utils import GetNormalizedSpecifications
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_specification_utils import GetCosineDistances
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_specification_utils import PrintSpecificationDataToCSV

class CosineSimilaritySensorPlacementAlgorithm(SensorPlacementAlgorithm):
    @classmethod
    def GetDefaultParameters(cls) -> Kratos.Parameters:
        return Kratos.Parameters("""{
            "type"                      : "cosine_similarity_sensor_placement",
            "output_to_vtu"             : true,
            "output_to_csv"             : true,
            "output_folder"             : "sensor_placement/cluster_average",
            "clustering_method"         : "average",
            "sensor_coverage_percentage": 5.0,
            "max_clustering_iterations" : 100,
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
        self.list_of_specifications:'list[KratosDT.Sensors.SensorSpecification]' = []
        self.parameters.RecursivelyValidateAndAssignDefaults(self.GetDefaultParameters())
        self.is_vtu_output = self.parameters["output_to_vtu"].GetBool()
        self.is_csv_output = self.parameters["output_to_csv"].GetBool()

    def Execute(self, list_of_specifications: 'list[KratosDT.Sensors.SensorSpecification]') -> None:
        self.list_of_specifications = list_of_specifications

        if len(self.list_of_specifications) == 0:
            raise RuntimeError("No specifications are given.")

        self.SmoothenSensitivityFields()
        first_specification = self.list_of_specifications[0]

        for k in first_specification.GetNodalExpressionsMap().keys():
            self.ComputeSensorPlacement(f"nodal/{k}", lambda x: x.GetNodalExpression(k))

        for k in first_specification.GetConditionExpressionsMap().keys():
            self.ComputeSensorPlacement(f"condition/{k}", lambda x: x.GetConditionExpression(k))

        for k in first_specification.GetElementExpressionsMap().keys():
            self.ComputeSensorPlacement(f"element/{k}", lambda x: x.GetElementExpression(k))

    def ComputeSensorPlacement(self, name: str, retrieval_method: 'typing.Callable[[KratosDT.Sensors.SensorSpecification], typing.Union[Kratos.Expression.NodalExpression, Kratos.Expression.ConditionExpression, Kratos.Expression.ElementExpression]]'):
        normalized_specifications = GetNormalizedSpecifications(self.list_of_specifications, retrieval_method)

        if len(normalized_specifications) == 0:
            raise RuntimeError("No specifications with zero vectors found.")

        dummy_cexp = retrieval_method(normalized_specifications[0])
        vtu_output = Kratos.VtuOutput(dummy_cexp.GetModelPart())

        # compute best sensor for each entity
        total_number_of_entities = len(retrieval_method(normalized_specifications[0]).GetContainer())

        entity_best_specs_list:'list[KratosDT.Sensors.SensorSpecification]' = []
        for entity_index in range(total_number_of_entities):
            entity_best_spec = max(normalized_specifications, key=lambda x: retrieval_method(x).Evaluate()[entity_index])
            entity_best_specs_list.append(entity_best_spec)

        unique_normalized_list = list(set(entity_best_specs_list))

        # get the entity domain sizes
        entity_domain_size_exp = retrieval_method(unique_normalized_list[0]).Clone()
        Kratos.Expression.EntityDomainSizeExpressionIO.Read(entity_domain_size_exp)
        total_domain_size = KratosOA.ExpressionUtils.Sum(entity_domain_size_exp)
        entity_domain_size_np_exp = entity_domain_size_exp.Evaluate()
        maximum_allowed_domain_size_per_cluster = total_domain_size * self.parameters["sensor_coverage_percentage"].GetDouble() / 100.0

        spec_cluster_dict = {1: unique_normalized_list}
        entity_index_cluster_dict = {1: range(total_number_of_entities)}

        output_path = Path(self.parameters["output_folder"].GetString()) / f"sensors/{name}"

        list_of_cluster_ids_cannot_be_divided: 'list[int]' = []

        clustering_iteration = 0
        clusters_to_break: 'list[int]' = [1]
        while clustering_iteration < self.parameters["max_clustering_iterations"].GetInt() and set(list_of_cluster_ids_cannot_be_divided) != set(clusters_to_break) and len(clusters_to_break) > 0:
            Kratos.Logger.PrintInfo(self.__class__.__name__,f"Clustering iteration {clustering_iteration}")
            clustering_iteration += 1

            for cluster_id_to_break in clusters_to_break:
                if not cluster_id_to_break in list_of_cluster_ids_cannot_be_divided:
                    number_of_sensors_in_cluster = len(spec_cluster_dict[cluster_id_to_break])
                    if number_of_sensors_in_cluster > 1:
                        Kratos.Logger.PrintInfo("", f"\tBreaking the cluster {cluster_id_to_break} having {number_of_sensors_in_cluster} sensor specifications...")
                        new_clusters = self.ClusterListOfSpecifications(2, spec_cluster_dict[cluster_id_to_break], retrieval_method)
                        spec_cluster_dict[cluster_id_to_break] = new_clusters[1]
                        if 2 in new_clusters.keys():
                            spec_cluster_dict[max(spec_cluster_dict.keys()) + 1] = new_clusters[2]
                        else:
                            list_of_cluster_ids_cannot_be_divided.append(cluster_id_to_break)
                    else:
                        list_of_cluster_ids_cannot_be_divided.append(cluster_id_to_break)
            entity_index_cluster_dict = self.ComputeEntityClustering(entity_best_specs_list, spec_cluster_dict)

            clusters_to_break = []
            for cluster_id, _ in spec_cluster_dict.items():
                # get the domain size of the cluster
                entity_indices = entity_index_cluster_dict[cluster_id]
                cluster_domain_size = np.sum(np.take(entity_domain_size_np_exp, entity_indices))

                if cluster_domain_size > maximum_allowed_domain_size_per_cluster:
                    clusters_to_break.append(cluster_id)

            # now look for the best sensor in each cluster which has the highest minimum for all the
            # entities which it relates to
            list_of_best_sensors: 'list[KratosDT.Sensors.SensorSpecification]' = []
            entity_cluster_data = [-1] * total_number_of_entities
            for cluster_id, list_of_specs in spec_cluster_dict.items():
                for spec in list_of_specs:
                    spec.SetValue(KratosDT.SENSOR_CLUSTER_ID, cluster_id)
                entity_indices = entity_index_cluster_dict[cluster_id]
                for entity_index in entity_indices:
                    entity_cluster_data[entity_index] = cluster_id
                best_spec = max(list_of_specs, key=lambda x: np.min(np.take(retrieval_method(x).Evaluate(), entity_indices)))
                list_of_best_sensors.append(best_spec)

            vtu_output.ClearCellContainerExpressions()
            vtu_output.ClearNodalContainerExpressions()

            cexp = dummy_cexp.Clone()
            Kratos.Expression.CArrayExpressionIO.Read(cexp, np.array(entity_cluster_data, dtype=np.int32))
            vtu_output.AddContainerExpression("entity_cluster_id", cexp.Clone())

            if len(list_of_best_sensors) > 0:
                heat_map = retrieval_method(list_of_best_sensors[0]).Clone()
                for spec in list_of_best_sensors[1:]:
                    heat_map += retrieval_method(spec)
                heat_map /= KratosOA.ExpressionUtils.NormL2(heat_map)
                vtu_output.AddContainerExpression("heat_map", heat_map)

            PrintSpecificationDataToCSV(list_of_best_sensors, output_path / f"best_placement_{clustering_iteration:05d}.csv")
            vtu_output.PrintOutput(str(output_path / f"heat_map_{clustering_iteration:05d}"))

        Kratos.Logger.PrintInfo(self.__class__.__name__, f"Found {len(list_of_best_sensors)} clusters with sensors.")

        if len(clusters_to_break) > 0:
            Kratos.Logger.PrintWarning(self.__class__.__name__, f"The max_clustering_iterations and/or unbreakable clusters limit reached without breaking {len(clusters_to_break)} clusters.\n\t Following clusters does not satisfy the prescribed coverage area:")
            for cluster_id_to_break in clusters_to_break:
                entity_indices = entity_index_cluster_dict[cluster_id_to_break]
                cluster_domain_size = np.sum(np.take(entity_domain_size_np_exp, entity_indices))
                Kratos.Logger.PrintInfo("", f"\t\t Cluster {cluster_id_to_break} - domain size = {cluster_domain_size * 100.0 / total_domain_size:0.3f} %")

        PrintSpecificationDataToCSV(unique_normalized_list, output_path / f"clusters.csv")

    def ClusterListOfSpecifications(self, number_of_cluster: int, list_of_specifications: 'list[KratosDT.Sensors.SensorSpecification]', retrieval_method: 'typing.Callable[[KratosDT.Sensors.SensorSpecification], typing.Union[Kratos.Expression.NodalExpression, Kratos.Expression.ConditionExpression, Kratos.Expression.ElementExpression]]') -> 'dict[int, list[KratosDT.Sensors.SensorSpecification]]':
        spec_cosine_distances = GetCosineDistances(list_of_specifications, retrieval_method)
        spec_cosine_linkage = sch.linkage(spec_cosine_distances, method=self.parameters["clustering_method"].GetString())
        spec_cosine_clusters = sch.fcluster(spec_cosine_linkage, number_of_cluster, 'maxclust')

        cluster_dict: 'dict[int, list[KratosDT.Sensors.SensorSpecification]]' = {}
        for i, cluster_id in enumerate(spec_cosine_clusters):
            if cluster_id not in cluster_dict.keys():
                cluster_dict[cluster_id] = []
            cluster_dict[cluster_id].append(list_of_specifications[i])

        return cluster_dict

    def ComputeEntityClustering(self, entity_best_specs_list: 'list[KratosDT.Sensors.SensorSpecification]', cluster_dict: 'dict[int, list[KratosDT.Sensors.SensorSpecification]]'):
        entity_clusters:'dict[int, list[int]]' = {}
        for cluster_id, specs_list in cluster_dict.items():
            entity_clusters[cluster_id] = []
            for entity_index, entity_best_spec in enumerate(entity_best_specs_list):
                if entity_best_spec in specs_list:
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

        for specification in self.list_of_specifications:
            if self.is_vtu_output:
                vtu_output.ClearCellContainerExpressions()
                vtu_output.ClearNodalContainerExpressions()
            for k, v in specification.GetNodalExpressionsMap().items():
                filtered_field = nodal_filter.FilterIntegratedField(v.Abs())
                if self.is_vtu_output:
                    vtu_output.AddContainerExpression(f"{k}_raw", v.Clone())
                    vtu_output.AddContainerExpression(f"{k}_filtered", filtered_field.Clone())

                v.SetExpression(filtered_field.GetExpression())

            for k, v in specification.GetConditionExpressionsMap().items():
                filtered_field = condition_filter.FilterIntegratedField(v.Abs())
                if self.is_vtu_output:
                    vtu_output.AddContainerExpression(f"{k}_raw", v.Clone())
                    vtu_output.AddContainerExpression(f"{k}_filtered", filtered_field.Clone())

                v.SetExpression(filtered_field.GetExpression())

            for k, v in specification.GetElementExpressionsMap().items():
                filtered_field = element_filter.FilterIntegratedField(v.Abs())
                if self.is_vtu_output:
                    vtu_output.AddContainerExpression(f"{k}_raw", v.Clone())
                    vtu_output.AddContainerExpression(f"{k}_filtered", filtered_field.Clone())

                v.SetExpression(filtered_field.GetExpression())

            if self.is_vtu_output:
                vtu_output.PrintOutput(str(vtu_output_path / f"{specification.GetName()}_{specification.Id:05d}"))

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
        for specification in self.list_of_specifications:
            if filter_type == KratosOA.NodalExplicitFilter:
                expressions_list = specification.GetNodalExpressionsMap().values()
                filter_radius_exp_type = Kratos.Expression.NodalExpression
            elif filter_type == KratosOA.ConditionExplicitFilter:
                expressions_list = specification.GetConditionExpressionsMap().values()
                filter_radius_exp_type = Kratos.Expression.ConditionExpression
            elif filter_type == KratosOA.ElementExplicitFilter:
                expressions_list = specification.GetElementExpressionsMap().values()
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





