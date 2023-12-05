import typing
from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.DigitalTwinApplication as KratosDT
from KratosMultiphysics.DigitalTwinApplication.sensor_placement_algorithms.sensor_placement_algorithm import SensorPlacementAlgorithm
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import GetNormalizedSensorViews
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import PrintSensorListToCSV
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import PrintSensorListToJson
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import GetFilter, GetSensorCoverageMasks, ComputeMinimumDistanceSquare
from KratosMultiphysics.DigitalTwinApplication.utilities.data_utils import SensorViewUnionType
from KratosMultiphysics.DigitalTwinApplication.utilities.data_utils import SensorViewUnionType
from KratosMultiphysics.DigitalTwinApplication.utilities.expression_utils import ExpressionFilterUnionType, ExpressionUnionType

ClusterUnionType = typing.Union[
                            KratosDT.ClusterUtils.NodalSensorViewCluster,
                            KratosDT.ClusterUtils.ConditionSensorViewCluster,
                            KratosDT.ClusterUtils.ElementSensorViewCluster
                        ]
class MinimumRedundancySensorPlacementAlgorithm(SensorPlacementAlgorithm):
    @classmethod
    def GetDefaultParameters(cls) -> Kratos.Parameters:
        return Kratos.Parameters("""{
            "type"                             : "cosine_similarity_sensor_placement",
            "output_to_vtu"                    : true,
            "output_to_csv"                    : true,
            "output_folder"                    : "sensor_placement",
            "required_maximum_overall_coverage": 0.9,
            "required_minimum_sensor_coverage" : 0.1,
            "number_of_redundant_sensors"      : 5,
            "maximum_number_of_iterations"     : 100,
            "filtering"                        : {}
        }""")

    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters) -> None:
        self.model = model
        self.parameters = parameters
        self.list_of_sensors:'list[KratosDT.Sensors.Sensor]' = []
        self.parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())
        self.is_vtu_output = self.parameters["output_to_vtu"].GetBool()
        self.is_csv_output = self.parameters["output_to_csv"].GetBool()
        self.required_maximum_overall_coverage = self.parameters["required_maximum_overall_coverage"].GetDouble()
        self.required_minimum_sensor_coverage = self.parameters["required_minimum_sensor_coverage"].GetDouble()
        self.maximum_number_of_iterations = self.parameters["maximum_number_of_iterations"].GetInt()
        self.number_of_redundant_sensors = self.parameters["number_of_redundant_sensors"].GetInt()

    def Execute(self, list_of_sensors: 'list[KratosDT.Sensors.Sensor]') -> None:
        self.list_of_sensors = list_of_sensors

        if len(self.list_of_sensors) == 0:
            raise RuntimeError("No specifications are given.")

        self.SmoothenSensitivityFields()

        first_specification = self.list_of_sensors[0]

        for k in first_specification.GetNodalExpressionsMap().keys():
            self.domain_sensor_view_cluster_type = KratosDT.ClusterUtils.NodalSensorViewClusterData
            self.cluster_type = KratosDT.ClusterUtils.NodalSensorViewCluster
            self.ComputeSensorPlacement(f"nodal/{k}", [KratosDT.Sensors.NodalSensorView(sensor, k) for sensor in list_of_sensors])

        for k in first_specification.GetConditionExpressionsMap().keys():
            self.domain_sensor_view_cluster_type = KratosDT.ClusterUtils.ConditionSensorViewClusterData
            self.cluster_type = KratosDT.ClusterUtils.ConditionSensorViewCluster
            self.ComputeSensorPlacement(f"condition/{k}", [KratosDT.Sensors.ConditionSensorView(sensor, k) for sensor in list_of_sensors])

        for k in first_specification.GetElementExpressionsMap().keys():
            self.domain_sensor_view_cluster_type = KratosDT.ClusterUtils.ElementSensorViewClusterData
            self.cluster_type = KratosDT.ClusterUtils.ElementSensorViewCluster
            self.ComputeSensorPlacement(f"element/{k}", [KratosDT.Sensors.ElementSensorView(sensor, k) for sensor in list_of_sensors])

    def ComputeSensorPlacement(self, name: str, list_of_sensor_views: 'list[SensorViewUnionType]'):
        # prepare for output
        if self.is_vtu_output or self.is_csv_output:
            output_path = Path(self.parameters["output_folder"].GetString()) / f"sensors/{name}"
            output_path.mkdir(parents=True, exist_ok=True)

        normalized_sensor_views = GetNormalizedSensorViews(list_of_sensor_views)
        dummy_exp = normalized_sensor_views[0].GetContainerExpression()
        domain_size_exp = dummy_exp.Clone()
        Kratos.Expression.EntityDomainSizeExpressionIO.Read(domain_size_exp)
        coverage_masks = GetSensorCoverageMasks(normalized_sensor_views, domain_size_exp)
        total_domain_size = KratosDT.SensorUtils.Sum(domain_size_exp)

        if self.is_vtu_output:
            vtu_output = Kratos.VtuOutput(normalized_sensor_views[0].GetContainerExpression().GetModelPart())
            for coverage_mask, sensor_view in zip(coverage_masks, normalized_sensor_views):
                vtu_output.AddContainerExpression(sensor_view.GetSensor().GetName(), coverage_mask)
            vtu_output.PrintOutput(str(output_path / "sensor_coverage"))

        if self.is_csv_output:
            PrintSensorListToCSV(output_path / "sensor_coverage.csv", [sensor_view.GetSensor() for sensor_view in normalized_sensor_views], ["type", "name", "location", "value", "SENSOR_COVERAGE"])

        overall_coverage_mask = (dummy_exp * 0.0).Flatten()
        max_coverage = 0.0
        iteration = 1
        list_of_selected_views: 'list[SensorViewUnionType]' = []
        list_of_selected_coverage_masks: 'list[ExpressionUnionType]' = []
        while (max_coverage < self.required_maximum_overall_coverage and iteration <= self.maximum_number_of_iterations):
            # find the sensor mask which increases coverage to the max
            max_coverage = 0
            best_coverage_index = -1
            max_min_distance = 1e+9
            for i, coverage_mask in enumerate(coverage_masks):
                # computing the additional area which is covered by this coverage mask
                modified_overall_coverage_mask = overall_coverage_mask + coverage_mask
                KratosDT.ControlUtils.ClipContainerExpression(modified_overall_coverage_mask, 0, 1)
                coverage = KratosDT.SensorUtils.Sum(domain_size_exp.Scale(modified_overall_coverage_mask)) / total_domain_size

                if abs(coverage - max_coverage) < 1e-9:
                    # both the sensors cover the same additional area
                    # hence we choose the most distanced one to have better redundancy
                    current_min_distance = ComputeMinimumDistanceSquare(normalized_sensor_views[i], list_of_selected_views)
                    if current_min_distance > max_min_distance:
                        max_min_distance = current_min_distance
                        best_coverage_index = i
                        max_coverage = coverage

                elif max_coverage < coverage:
                    # current coverage is better than the max coverage we found so far
                    # hence, getting that max coverage
                    max_coverage = coverage
                    best_coverage_index = i
                    max_min_distance = ComputeMinimumDistanceSquare(normalized_sensor_views[i], list_of_selected_views)

            if best_coverage_index != -1:
                list_of_selected_views.append(normalized_sensor_views[best_coverage_index])
                list_of_selected_coverage_masks.append(coverage_masks[best_coverage_index])
                overall_coverage_mask += coverage_masks[best_coverage_index]
                overall_coverage_mask = overall_coverage_mask.Flatten()

                if self.is_vtu_output:
                    vtu_output.ClearCellContainerExpressions()
                    vtu_output.ClearNodalContainerExpressions()
                    vtu_output.AddContainerExpression("sensor_coverage", coverage_masks[best_coverage_index].Clone())
                    vtu_output.AddContainerExpression("redundancy", overall_coverage_mask.Clone())
                    vtu_output.PrintOutput(str(output_path / f"redundancy_iteration_{iteration:05d}"))

                if self.is_csv_output:
                    PrintSensorListToCSV(output_path / f"redundancy_iteration_{iteration:05d}.csv", [sensor_view.GetSensor() for sensor_view in list_of_selected_views], ["type", "name", "location", "value"])

                Kratos.Logger.PrintInfo(self.__class__.__name__, f"Found sensor {normalized_sensor_views[best_coverage_index].GetSensor().GetName()} to increase coverage to {max_coverage * 100.0:6.3f}%.")

                # now remove them from the remaining searches
                del normalized_sensor_views[best_coverage_index]
                del coverage_masks[best_coverage_index]

            iteration += 1

        Kratos.Logger.PrintInfo(self.__class__.__name__, f"Found {len(list_of_selected_views)} sensors covering {max_coverage * 100.0:6.3f}% domain.")

        # now we have found the base sensors which has a minimum coverage of self.required_maximum_overall_coverage
        # now we can sub divide all the sensor coverage areas to sub areas so it will help localizing
        # the damage areas
        sensor_improvement_gap = total_domain_size / self.number_of_redundant_sensors
        index = -1
        number_of_base_sensor_views = len(list_of_selected_views)
        while index < number_of_base_sensor_views - 1:
            index += 1
            sensor_view = list_of_selected_views[index]
            coverage_mask = list_of_selected_coverage_masks[index]
            sensor_area = KratosDT.SensorUtils.Sum(domain_size_exp.Scale(coverage_mask))
            Kratos.Logger.PrintInfo("", f"--- Finding redundant sensors for {sensor_view.GetSensor().GetName()}...")

            required_max_coverage_area = (sensor_area - sensor_improvement_gap) / total_domain_size
            required_min_coverage_area = required_max_coverage_area - sensor_improvement_gap / total_domain_size
            while (required_min_coverage_area >= 0.0):
                # get the list of possibe sensors covering the chose gap between
                # min and max
                list_of_potential_sensor_views: 'list[int]' = []
                for i, coverage_mask in enumerate(coverage_masks):
                    redundant_coverage_domain = KratosDT.SensorUtils.Sum(domain_size_exp.Scale(coverage_mask)) / total_domain_size
                    if redundant_coverage_domain <= required_max_coverage_area and redundant_coverage_domain >= required_min_coverage_area and redundant_coverage_domain > self.required_minimum_sensor_coverage:
                        list_of_potential_sensor_views.append(i)

                # now find the most distanced sensor from the potential sensors
                max_distance = 0
                most_distance_sensor_view_index = -1
                for i in list_of_potential_sensor_views:
                    min_distance = ComputeMinimumDistanceSquare(normalized_sensor_views[i], list_of_selected_views)
                    if min_distance > max_distance:
                        max_distance = min_distance
                        most_distance_sensor_view_index = i

                if most_distance_sensor_view_index != -1:
                    list_of_selected_views.append(normalized_sensor_views[most_distance_sensor_view_index])
                    list_of_selected_coverage_masks.append(coverage_masks[most_distance_sensor_view_index])
                    overall_coverage_mask += coverage_masks[most_distance_sensor_view_index]
                    overall_coverage_mask = overall_coverage_mask.Flatten()

                    Kratos.Logger.PrintInfo("", f"------ Found redundant sensor {normalized_sensor_views[most_distance_sensor_view_index].GetSensor().GetName()} for {sensor_view.GetSensor().GetName()} with the coverage in [{required_min_coverage_area * 100.0}, {required_max_coverage_area * 100.0}].")

                    if self.is_vtu_output:
                        vtu_output.ClearCellContainerExpressions()
                        vtu_output.ClearNodalContainerExpressions()
                        vtu_output.AddContainerExpression("sensor_coverage", coverage_masks[most_distance_sensor_view_index].Clone())
                        vtu_output.AddContainerExpression("redundancy", overall_coverage_mask.Clone())
                        vtu_output.PrintOutput(str(output_path / f"redundancy_iteration_{iteration:05d}"))

                    if self.is_csv_output:
                        PrintSensorListToCSV(output_path / f"redundancy_iteration_{iteration:05d}.csv", [sensor_view.GetSensor() for sensor_view in list_of_selected_views], ["type", "name", "location", "value"])

                    del coverage_masks[most_distance_sensor_view_index]
                    del normalized_sensor_views[most_distance_sensor_view_index]

                    iteration += 1

                required_max_coverage_area -= sensor_improvement_gap / total_domain_size
                required_min_coverage_area -= sensor_improvement_gap / total_domain_size

        list_of_best_sensors = [sensor_view.GetSensor() for sensor_view in list_of_selected_views]
        PrintSensorListToJson(output_path / "best_sensor_data.json", list_of_best_sensors)
        PrintSensorListToCSV(output_path / "best_sensor_data.csv", list_of_best_sensors, ["type", "name", "location", "value"])

    def SmoothenSensitivityFields(self) -> None:
        mp_nodal, nodal_filter = self.__GetFilter(Kratos.Globals.DataLocation.NodeHistorical)
        mp_condition, condition_filter = self.__GetFilter(Kratos.Globals.DataLocation.Condition)
        mp_element, element_filter = self.__GetFilter(Kratos.Globals.DataLocation.Element)

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

    def __GetFilter(self, filter_field_type: Kratos.Globals.DataLocation) -> 'tuple[Kratos.ModelPart, ExpressionFilterUnionType]':
        field_filter: ExpressionFilterUnionType = None
        model_part: Kratos.ModelPart = None
        for sensor in self.list_of_sensors:
            if filter_field_type in [Kratos.Globals.DataLocation.NodeHistorical, Kratos.Globals.DataLocation.NodeNonHistorical]:
                expressions_dict = sensor.GetNodalExpressionsMap()
            elif filter_field_type == Kratos.Globals.DataLocation.Condition:
                expressions_dict = sensor.GetConditionExpressionsMap()
            elif filter_field_type == Kratos.Globals.DataLocation.Element:
                expressions_dict = sensor.GetElementExpressionsMap()

            for v in expressions_dict.values():
                if field_filter is None:
                    model_part = v.GetModelPart()
                    field_filter = GetFilter(model_part, filter_field_type, self.parameters["filtering"])
                    break

            if field_filter is not None:
                break

        return model_part, field_filter





