import math
import numpy
from itertools import combinations
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm import Algorithm
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import time_decorator
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import OptimizationAlgorithmTimeLogger
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import GetSensors

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
    return AlgorithmBruteForceOptimalSensorPlacement(model, parameters, optimization_problem)

class AlgorithmBruteForceOptimalSensorPlacement(Algorithm):
    """
        A brute force algorithm to find optimal sensors from a given set of sensors.
        This algorithm tries all the available combinations to find the best sensor combination.
    """

    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "module"            : "KratosMultiphysics.SystemIdentificationApplication.algorithms",
            "type"              : "algorithm_brute_force_optimal_sensor_placement",
            "settings"          : {
                "sensor_group_name"                      : "",
                "sensor_mask_name"                       : "",
                "number_of_expected_sensors"             : 4,
                "cluster_size_ratio_comparison_tolerance": 1e-9
            }
        }""")

    def __init__(self, model:Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        self.model = model
        self.parameters = parameters
        self._optimization_problem = optimization_problem

        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())

        settings = parameters["settings"]
        settings.ValidateAndAssignDefaults(self.GetDefaultParameters()["settings"])

        self.sensor_group_name = settings["sensor_group_name"].GetString()
        self.sensor_mask_name = settings["sensor_mask_name"].GetString()
        self.number_of_expected_sensors = settings["number_of_expected_sensors"].GetInt()
        self.cluster_size_ratio_comparison_tolerance = settings["cluster_size_ratio_comparison_tolerance"].GetDouble()

        ComponentDataView("algorithm", self._optimization_problem).SetDataBuffer(self.GetMinimumBufferSize())

    def GetMinimumBufferSize(self) -> int:
        return 2

    def Check(self):
        pass

    def Initialize(self):
        self.algorithm_data = ComponentDataView("algorithm", self._optimization_problem)

    def Finalize(self):
        pass

    @time_decorator()
    def Output(self) -> KratosOA.CollectiveExpression:
        for process in self._optimization_problem.GetListOfProcesses("output_processes"):
            if process.IsOutputStep():
                process.PrintOutput()

    @time_decorator()
    def Solve(self):
        buffered_data = self.algorithm_data.GetBufferedData()
        sensor_group_data = ComponentDataView(self.sensor_group_name, self._optimization_problem)
        list_of_sensors = GetSensors(sensor_group_data)

        domain_size_expression = list_of_sensors[0].GetContainerExpression(self.sensor_mask_name).Clone()
        Kratos.Expression.DomainSizeExpressionIO.Read(domain_size_expression)

        total_domain_size = float(numpy.sum(domain_size_expression.Evaluate()))
        buffered_data["best_maximum_cluster_size"] = total_domain_size
        buffered_data["current_maximum_cluster_size"] = total_domain_size
        buffered_data["best_maximum_cluster_size_ratio"] = 1.0
        buffered_data["current_maximum_cluster_size_ratio"] = 1.0
        buffered_data["best_sensor_ids"] = "[]"
        buffered_data["current_sensor_ids"] = "[]"

        list_of_sensor_combinations = combinations(list_of_sensors, self.number_of_expected_sensors)
        number_of_combinations = math.comb(len(list_of_sensors), self.number_of_expected_sensors)
        for i, sensor_combination in enumerate(list_of_sensor_combinations):
            with OptimizationAlgorithmTimeLogger("AlgorithmBruteForceOptimalSensorPlacement",self._optimization_problem.GetStep()):
                self._InitializeIteration()

                Kratos.VariableUtils().SetNonHistoricalVariableToZero(KratosSI.SENSOR_STATUS, self.model[self.sensor_group_name].Nodes)

                buffered_data.SetValue("current_sensor_ids", "[" + ":".join([str(sensor.GetNode().Id) for sensor in sensor_combination]) + "]", overwrite=True)

                Kratos.Logger.PrintInfo(self.__class__.__name__, f"Trying sensor combination {i + 1} / {number_of_combinations}: " + buffered_data[f"current_sensor_ids"])

                list_of_masks = [sensor.GetContainerExpression(self.sensor_mask_name) for sensor in sensor_combination]

                current_max_cluster_size = self.__CalculateMaxClusterSize(domain_size_expression, list_of_masks)
                current_max_cluster_size_ratio = current_max_cluster_size / total_domain_size
                buffered_data.SetValue("current_maximum_cluster_size", current_max_cluster_size, overwrite=True)
                buffered_data.SetValue("current_maximum_cluster_size_ratio", current_max_cluster_size_ratio, overwrite=True)

                if (i > 0):
                    buffered_data.SetValue("best_sensor_ids", buffered_data.GetValue("best_sensor_ids", 1), overwrite=True)
                    buffered_data.SetValue("best_maximum_cluster_size", buffered_data.GetValue("best_maximum_cluster_size", 1), overwrite=True)
                    buffered_data.SetValue("best_maximum_cluster_size_ratio", buffered_data.GetValue("best_maximum_cluster_size_ratio", 1), overwrite=True)
                if (i == 0) or (buffered_data.GetValue("best_maximum_cluster_size_ratio", 1) - current_max_cluster_size_ratio > self.cluster_size_ratio_comparison_tolerance):
                    buffered_data.SetValue("best_sensor_ids", buffered_data[f"current_sensor_ids"], overwrite=True)
                    buffered_data.SetValue("best_maximum_cluster_size", current_max_cluster_size, overwrite=True)
                    buffered_data.SetValue("best_maximum_cluster_size_ratio", current_max_cluster_size_ratio, overwrite=True)
                elif abs(buffered_data.GetValue("best_maximum_cluster_size_ratio", 1) - current_max_cluster_size_ratio) < self.cluster_size_ratio_comparison_tolerance:
                    buffered_data.SetValue("best_sensor_ids", "[" + buffered_data["best_sensor_ids"][1:-1] + ":" + buffered_data[f"current_sensor_ids"][1:-1] + "]", overwrite=True)
                    buffered_data.SetValue("best_maximum_cluster_size", buffered_data.GetValue("best_maximum_cluster_size", 0), overwrite=True)
                    buffered_data.SetValue("best_maximum_cluster_size_ratio", buffered_data.GetValue("best_maximum_cluster_size_ratio", 0), overwrite=True)

                for sensor_id_str in buffered_data["best_sensor_ids"][1:-1].split(":"):
                    self.model[self.sensor_group_name].GetNode(int(sensor_id_str)).SetValue(KratosSI.SENSOR_STATUS, 1.0)

                self._FinalizeIteration()

                self.Output()

                self._optimization_problem.AdvanceStep()

        return True

    @staticmethod
    def __GenerateCombination(current_item: 'list[int]', N: int, position: int, max_positions: int):
        if position < max_positions:
            for i in range(position, N):
                current_item[position] = i
                AlgorithmBruteForceOptimalSensorPlacement.__GenerateCombination(current_item, N, position + 1, max_positions)

    def __CalculateMaxClusterSize(self, domain_size_expression: ContainerExpressionTypes, list_of_masks: 'list[ContainerExpressionTypes]') -> float:
        # first calculate the coverage exp
        coverage_exp = list_of_masks[0].Clone()
        for mask in list_of_masks:
            coverage_exp = KratosSI.MaskUtils.Union(coverage_exp, mask)
        self.algorithm_data.GetUnBufferedData().SetValue("coverage", coverage_exp.Clone(), overwrite=True)

        # now calculate the clustering
        cluster_info_list: 'list[tuple[list[int], ContainerExpressionTypes]]' = KratosSI.MaskUtils.ClusterMasks(list_of_masks)
        sorted_cluster_info = sorted(cluster_info_list, key=lambda x: numpy.sum((x[1] * domain_size_expression).Evaluate()))
        overall_clustering = coverage_exp * 0.0
        for i, (_, cluster_exp) in enumerate(sorted_cluster_info):
            overall_clustering += cluster_exp * i
        self.algorithm_data.GetUnBufferedData().SetValue("clustering", overall_clustering.Clone(), overwrite=True)

        # now return the max cluster size
        return float(numpy.sum((sorted_cluster_info[-1][1] * domain_size_expression).Evaluate()))
