import typing
import scipy.cluster.hierarchy as sch
from pathlib import Path
from functools import reduce

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
        # categorize specs based on name
        dict_of_specs: 'dict[str, list[KratosDT.Sensors.SensorSpecification]]' = {}
        for spec in list_of_specifications:
            if not spec.GetName() in dict_of_specs.keys():
                dict_of_specs[spec.GetName()] = []
            dict_of_specs[spec.GetName()].append(spec)

        # now find the best specs from each category
        self.list_of_specifications: 'list[KratosDT.Sensors.SensorSpecification]' = []
        for _, specs in dict_of_specs.items():
            sorted_specs = sorted(specs, key=lambda x: abs(x.GetSensorValue()), reverse=True)
            self.list_of_specifications.extend(sorted_specs[:int(len(sorted_specs) * self.parameters["percentage_of_best_sensors_to_consider"].GetDouble())])

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
        number_of_sensors = self.parameters["number_of_sensors"].GetInt()
        normalized_specifications = GetNormalizedSpecifications(self.list_of_specifications, retrieval_method)
        distances = GetCosineDistances(normalized_specifications, retrieval_method)
        linkage = sch.linkage(distances, method=self.parameters["clustering_method"].GetString())
        clusters = sch.fcluster(linkage, number_of_sensors, 'maxclust')

        # group clusters
        clustered_specs:'dict[int, list[tuple[int, KratosDT.Sensors.SensorSpecification]]]' = {}
        for i, cluster_id in enumerate(clusters):
            if cluster_id not in clustered_specs.keys():
                clustered_specs[cluster_id] = []
            spec = normalized_specifications[i]
            spec.SetValue(KratosDT.SENSOR_CLUSTER_ID, cluster_id)
            clustered_specs[cluster_id].append((i, spec))

        # order clusters with max average sensor values
        sorted_cluster_ids = sorted(
                                range(1, number_of_sensors + 1),
                                key=lambda x: (reduce(lambda accumulator, spec_data: accumulator + abs(spec_data[1].GetSensorValue()), clustered_specs[x], 0.0)) / len(clustered_specs[x]),
                                reverse=True)

        # now get the first sensor from the cluster which has the highest average sensor value, from that the sensor with the
        # highest sensor value
        best_sensors: 'list[tuple[int, KratosDT.Sensors.SensorSpecification]]' = [max(clustered_specs[sorted_cluster_ids[0]], key=lambda x: abs(x[1].GetSensorValue()))]
        # now find the rest of the sensors which has the highest cosine distance to the already found ones.
        number_of_normalized_specs = len(normalized_specifications)
        for cluster_id in sorted_cluster_ids[1:]:
            # get the cluster spec data
            cluster_sensor_spec_data = clustered_specs[cluster_id]

            max_combined_distance = 0.0
            best_cluster_spec_data:'tuple[int, KratosDT.Sensors.SensorSpecification]' = (-1, None)

            # iterate through every cluster spec to find the best orthogonal
            # spec for already found specs
            for cluster_spec_index, cluster_spec in cluster_sensor_spec_data:
                min_combined_distance = 1e+16
                for best_spec_index, _ in best_sensors:
                    if best_spec_index > cluster_spec_index:
                        i = cluster_spec_index
                        j = best_spec_index
                    elif best_spec_index < cluster_spec_index:
                        i = best_spec_index
                        j = cluster_spec_index
                    else:
                        raise RuntimeError("Trying to add a duplicate sensor")
                    current_combined_distance = distances[number_of_normalized_specs * i + j - ((i + 2) * (i + 1)) // 2]

                    if min_combined_distance > current_combined_distance:
                        min_combined_distance = current_combined_distance

                if max_combined_distance < min_combined_distance:
                    max_combined_distance = min_combined_distance
                    best_cluster_spec_data = (cluster_spec_index, cluster_spec)

            best_sensors.append(best_cluster_spec_data)

        output_path = Path(self.parameters["output_folder"].GetString()) / f"sensors/{name}"
        PrintSpecificationDataToCSV(normalized_specifications, output_path / f"clusters_{number_of_sensors:05d}.csv")
        PrintSpecificationDataToCSV([spec for _, spec in best_sensors], output_path / f"best_placement_{number_of_sensors:05d}.csv")

        if self.is_vtu_output and len(best_sensors) > 0:
            model_part = retrieval_method(best_sensors[0][1]).GetModelPart()
            vtu_output = Kratos.VtuOutput(model_part)
            heat_map = retrieval_method(best_sensors[0][1])
            for _, spec in best_sensors[1:]:
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





