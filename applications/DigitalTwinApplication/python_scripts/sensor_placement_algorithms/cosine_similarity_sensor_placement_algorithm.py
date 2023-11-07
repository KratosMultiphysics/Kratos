import typing
import scipy.cluster.hierarchy as sch
import numpy as np
from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.DigitalTwinApplication as KratosDT
from KratosMultiphysics.DigitalTwinApplication.sensor_placement_algorithms.sensor_placement_algorithm import SensorPlacementAlgorithm

class CosineSimilaritySensorPlacementAlgorithm(SensorPlacementAlgorithm):
    @classmethod
    def GetDefaultParameters(cls) -> Kratos.Parameters:
        return Kratos.Parameters("""{
            "type"             : "cosine_similarity_sensor_placement",
            "number_of_sensors": 20,
            "output_to_vtu"    : true,
            "output_to_csv"    : true,
            "output_folder"    : "sensor_placement/cluster_average",
            "clustering_method": "average",
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

        self.SmoothenSensitivityFields()

        retrieval_mathod: 'typing.Callable[[KratosDT.Sensors.SensorSpecification], typing.Union[Kratos.Expression.NodalExpression, Kratos.Expression.ConditionExpression, Kratos.Expression.ElementExpression]]' = lambda x: x.GetElementExpression("YOUNG_MODULUS_SENSITIVITY")

        data: 'list[np.ndarray]' = []
        indices: 'list[int]' = []
        for i, spec in enumerate(list_of_specifications):
            cexp = retrieval_mathod(spec)
            l2_norm = KratosOA.ExpressionUtils.NormL2(cexp)
            if l2_norm > 0.0:
                indices.append(i)
                data.append(cexp.Evaluate() / l2_norm)

        number_of_sensors = self.parameters["number_of_sensors"].GetInt()
        distances = sch.distance.pdist(data, metric="cosine")
        linkage = sch.linkage(distances, method=self.parameters["clustering_method"].GetString())
        clusters = sch.fcluster(linkage, number_of_sensors, 'maxclust')

        # group clusters
        clustered_specs:'dict[int, list[KratosDT.Sensors.SensorSpecification]]' = {}
        for i, cluster_id in enumerate(clusters):
            if cluster_id not in clustered_specs.keys():
                clustered_specs[cluster_id] = []
            spec = list_of_specifications[indices[i]]
            spec.SetValue(KratosDT.SENSOR_CLUSTER_ID, cluster_id)
            clustered_specs[cluster_id].append(spec)

        best_sensors: 'list[KratosDT.Sensors.SensorSpecification]' = []
        for _, list_of_cluster_specs in clustered_specs.items():
            best_sensors.append(max(list_of_cluster_specs, key=lambda x: abs(x.GetSensorValue())))

        output_path = Path(self.parameters["output_folder"].GetString())
        self.PrintSpecificationDataToCSV([list_of_specifications[i] for i in indices], output_path / f"sensors/clusters_{number_of_sensors:05d}.csv")
        self.PrintSpecificationDataToCSV(best_sensors, output_path / f"sensors/best_placement_{number_of_sensors:05d}.csv")

        if self.is_vtu_output and len(best_sensors) > 0:
            model_part = retrieval_mathod(best_sensors[0]).GetModelPart()
            vtu_output = Kratos.VtuOutput(model_part)
            heat_map = retrieval_mathod(best_sensors[0])
            for spec in best_sensors[1:]:
                heat_map += retrieval_mathod(spec)
            heat_map /= KratosOA.ExpressionUtils.NormL2(heat_map)
            vtu_output.AddContainerExpression("heat_map", heat_map)
            vtu_output.PrintOutput(str(output_path / f"sensors/heat_map_{number_of_sensors:05d}"))

    def PrintSpecificationDataToCSV(self, list_of_specifications:'list[KratosDT.Sensors.SensorSpecification]', output_file_name: Path) -> None:
        number_of_clusters = len(list_of_specifications)

        # do nothing if number of clusters is zero
        if number_of_clusters == 0:
            return

        output_file_name.parent.mkdir(exist_ok=True, parents=True)

        with open(str(output_file_name), "w") as csv_output:
            # write the headers of spec
            csv_output.write("#; name; value; location_x; location_y; location_z")
            # now write the headers of data containers
            list_of_vars = []
            for var_name in list_of_specifications[0].GetDataVariableNames():
                csv_output.write(f"; {var_name}")
                list_of_vars.append(Kratos.KratosGlobals.GetVariable(var_name))
            csv_output.write("\n")

            # now write the data
            sensor_id = 0
            for spec in list_of_specifications:
                sensor_id += 1
                loc = spec.GetLocation()
                csv_output.write(f"{sensor_id}; {spec.GetName()}; {spec.GetSensorValue()}; {loc[0]}; {loc[1]}; {loc[2]}")
                for var in list_of_vars:
                    csv_output.write(f"; {spec.GetValue(var)}")
                csv_output.write("\n")

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





