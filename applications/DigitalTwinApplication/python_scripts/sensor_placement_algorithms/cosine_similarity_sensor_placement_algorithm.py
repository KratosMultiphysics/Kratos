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
            "type"                                  : "cosine_similarity_sensor_placement",
            "number_of_sensors"                     : 20,
            "output_to_vtu"                         : true,
            "output_to_csv"                         : true,
            "output_folder"                         : "sensor_placement/cluster_average",
            "clustering_method"                     : "average",
            "percentage_of_best_sensors_to_consider": 0.90,
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
        # first cluster them based on the location of the sensor
        number_of_sensors = self.parameters["number_of_sensors"].GetInt()

        normalized_specifications = GetNormalizedSpecifications(self.list_of_specifications, retrieval_method)

        if len(normalized_specifications) == 0:
            raise RuntimeError("No specifications with zero vectors found.")

        # compute best sensor for each entity
        total_number_of_entities = len(retrieval_method(normalized_specifications[0]).GetContainer())

        entity_best_specs_list:'list[KratosDT.Sensors.SensorSpecification]' = []
        for entity_index in range(total_number_of_entities):
            entity_best_spec = max(normalized_specifications, key=lambda x: retrieval_method(x).Evaluate()[entity_index])
            entity_best_specs_list.append(entity_best_spec)

        unique_normalized_list = list(set(entity_best_specs_list))

        # now sensors based on their sensitvities which are used for updating
        spec_cosine_distances = GetCosineDistances(unique_normalized_list, retrieval_method)
        spec_cosine_linkage = sch.linkage(spec_cosine_distances, method=self.parameters["clustering_method"].GetString())
        spec_cosine_clusters = sch.fcluster(spec_cosine_linkage, number_of_sensors, 'maxclust')

        # now group the sensors in each cluster and assign SENSOR_CLUSTER_ID
        # assign the cluster id for the entity as well
        spec_clusters_dict:'dict[int, tuple[list[int], list[KratosDT.Sensors.SensorSpecification]]]' = {}
        entity_cluster_data: 'list[int]' = [-1] * total_number_of_entities
        for i, cluster_id in enumerate(spec_cosine_clusters):
            if cluster_id not in spec_clusters_dict.keys():
                spec_clusters_dict[cluster_id] = ([], [])
            spec = unique_normalized_list[i]
            spec.SetValue(KratosDT.SENSOR_CLUSTER_ID, cluster_id)
            spec_clusters_dict[cluster_id][1].append(spec)
            for entity_index in range(total_number_of_entities):
                if entity_best_specs_list[entity_index] == spec:
                    entity_cluster_data[entity_index] = cluster_id
                    spec_clusters_dict[cluster_id][0].append(entity_index)

        # now look for the best sensor in each cluster which has the highest minimum for all the
        # entities which it relates to
        list_of_best_sensors: 'list[KratosDT.Sensors.SensorSpecification]' = []
        for cluster_id, (entity_indices, list_of_specs) in spec_clusters_dict.items():
            best_spec = max(list_of_specs, key=lambda x: np.min(np.take(retrieval_method(x).Evaluate(), entity_indices)))
            list_of_best_sensors.append(best_spec)

        output_path = Path(self.parameters["output_folder"].GetString()) / f"sensors/{name}"
        PrintSpecificationDataToCSV(unique_normalized_list, output_path / f"clusters_{number_of_sensors:05d}.csv")
        PrintSpecificationDataToCSV(list_of_best_sensors, output_path / f"best_placement_{number_of_sensors:05d}.csv")

        if self.is_vtu_output:
            dummy_cexp = retrieval_method(entity_best_specs_list[0])
            vtu_output = Kratos.VtuOutput(dummy_cexp.GetModelPart())

            cexp = dummy_cexp.Clone()
            Kratos.Expression.CArrayExpressionIO.Read(cexp, np.array(entity_cluster_data, dtype=np.int32))
            vtu_output.AddContainerExpression("entity_cluster_id", cexp.Clone())

            Kratos.Expression.CArrayExpressionIO.Read(cexp, np.array([spec.Id for spec in entity_best_specs_list], dtype=np.int32))
            vtu_output.AddContainerExpression("entity_sensor_id", cexp.Clone())

            if len(list_of_best_sensors) > 0:
                heat_map = retrieval_method(list_of_best_sensors[0])
                for spec in list_of_best_sensors[1:]:
                    heat_map += retrieval_method(spec)
                heat_map /= KratosOA.ExpressionUtils.NormL2(heat_map)
                vtu_output.AddContainerExpression("heat_map", heat_map)

            vtu_output.PrintOutput(str(output_path / f"heat_map_{number_of_sensors:05d}"))


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





