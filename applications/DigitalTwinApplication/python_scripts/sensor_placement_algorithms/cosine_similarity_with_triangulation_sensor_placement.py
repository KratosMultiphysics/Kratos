import typing
from pathlib import Path
from scipy.spatial.distance import squareform


import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.DigitalTwinApplication as KratosDT
from KratosMultiphysics.DigitalTwinApplication.sensor_placement_algorithms.sensor_placement_algorithm import SensorPlacementAlgorithm
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import GetNormalizedSensorViews
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import GetCosineDistances
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import GetEuclideanDistances
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import PrintSensorListToCSV
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import PrintSensorListToJson
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import GetBestCoverageSensorView
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import ComputeHeatMap
from KratosMultiphysics.DigitalTwinApplication.utilities.data_utils import SensorViewUnionType

class CosineSimilarityWithTriangulationSensorPlacement(SensorPlacementAlgorithm):
    @classmethod
    def GetDefaultParameters(cls) -> Kratos.Parameters:
        return Kratos.Parameters("""{
            "type"                             : "least_cosine_euclidean_similarity_sensor_placement_algorithm",
            "output_to_vtu"                    : true,
            "output_to_csv"                    : true,
            "output_folder"                    : "sensor_placement/",
            "minimum_distance"                 : 0.0,
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
        self.minimum_distance = self.parameters["minimum_distance"].GetDouble()

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

        # get the best coverage sensor as the first sensor
        most_covered_sensor = GetBestCoverageSensorView(normalized_sensor_views)
        Kratos.Logger.PrintInfo(self.__class__.__name__, f"Found most covered sensor {most_covered_sensor.GetSensor().GetName()}")
        list_of_best_sensor_view_indices = [normalized_sensor_views.index(most_covered_sensor)]
        remaining_possible_sensor_indices = list(range(0, len(normalized_sensor_views)))
        del remaining_possible_sensor_indices[list_of_best_sensor_view_indices[0]]

        # calculate distances
        cosine_distances = GetCosineDistances(normalized_sensor_views)
        euclidean_distances = GetEuclideanDistances(normalized_sensor_views)
        cosine_distances = squareform(cosine_distances)
        euclidean_distances = squareform(euclidean_distances)

        output_path = Path(self.parameters["output_folder"].GetString()) / f"sensors/{name}"

        vtu_output = Kratos.VtuOutput(most_covered_sensor.GetContainerExpression().GetModelPart())
        heat_map = most_covered_sensor.GetContainerExpression().Clone()
        vtu_output.AddContainerExpression("heat_map", heat_map)

        found_sensor = True
        iteration = 0
        while found_sensor:
            iteration += 1
            # sort remaing sensors such that the min cosine distance between already found sensors is maximized
            remaining_possible_sensor_indices = sorted(remaining_possible_sensor_indices, key=lambda y: min(list_of_best_sensor_view_indices, key=lambda x: cosine_distances[x, y]), reverse=True)

            # now iterate and find the next best outside the minimum distance
            found_sensor = False
            for i, possible_sensor_index in enumerate(remaining_possible_sensor_indices):
                valid_possible_sensor = True
                for best_sensor_index in list_of_best_sensor_view_indices:
                    if euclidean_distances[best_sensor_index, possible_sensor_index] <= self.minimum_distance:
                        valid_possible_sensor = False

                if valid_possible_sensor:
                    del remaining_possible_sensor_indices[i]
                    list_of_best_sensor_view_indices.append(possible_sensor_index)
                    found_sensor = True
                    Kratos.Logger.PrintInfo(self.__class__.__name__, f"Found sensor {normalized_sensor_views[possible_sensor_index].GetSensor().GetName()}")
                    break

            list_of_best_sensors = [normalized_sensor_views[i].GetSensor() for i in list_of_best_sensor_view_indices]
            PrintSensorListToCSV(output_path / f"sensor_{iteration}.csv", list_of_best_sensors, ["type", "name", "location", "value"])
            heat_map.SetExpression(ComputeHeatMap([normalized_sensor_views[i] for i in list_of_best_sensor_view_indices]).GetExpression())
            vtu_output.PrintOutput(str(output_path / f"sensor_{iteration}"))

        list_of_best_sensors = [normalized_sensor_views[i].GetSensor() for i in list_of_best_sensor_view_indices]
        Kratos.Logger.PrintInfo(self.__class__.__name__, f"Found {len(list_of_best_sensors)} sensors.")
        PrintSensorListToJson(output_path / "best_sensor_data.json", list_of_best_sensors)
        PrintSensorListToCSV(output_path / "best_sensor_data.csv", list_of_best_sensors, ["type", "name", "location", "value"])

        heat_map.SetExpression(ComputeHeatMap([normalized_sensor_views[i] for i in list_of_best_sensor_view_indices]).GetExpression())
        vtu_output.PrintOutput(str(output_path / "heat_map"))

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





