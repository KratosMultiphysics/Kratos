from typing import Any
from math import sqrt

import KratosMultiphysics as Kratos

import KratosMultiphysics.StatisticsApplication as KratosStats
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.testing.utilities import ReadModelPart
from KratosMultiphysics.StatisticsApplication.test_utilities import CheckValues

class SpatialMethodTests(KratosUnittest.TestCase):
    class ValueNorm:
        def Evaluate(self, v: Any) -> Any:
            return v

    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 2

        # add all variable types solution step
        cls.model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.INITIAL_STRAIN)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.LOCAL_TANGENT_MATRIX)

        with KratosUnittest.WorkFolderScope(".", __file__):
            ReadModelPart("spatial_statistics_process/spatial_statistics_process", cls.model_part)

        def populate(container, setter_method, offset):
            for item in container:
                item_id = item.Id
                setter_method(item, Kratos.PRESSURE, -item_id - 1 - offset)
                setter_method(item, Kratos.VELOCITY, Kratos.Array3([item_id + 1 + offset, -item_id - 2 - offset, item_id + 3 + offset]))
                setter_method(item, Kratos.INITIAL_STRAIN, Kratos.Vector([item_id + 1 + offset, -item_id - 2 - offset, item_id + 3 + offset, -item_id - 4 - offset, item_id + 5 + offset]))
                setter_method(item, Kratos.LOCAL_TANGENT_MATRIX, Kratos.Matrix([[item_id + 1 + offset, -item_id - 2 - offset], [item_id + 3 + offset, -item_id - 4 - offset]]))

        populate(cls.model_part.Nodes, lambda x, y, z: x.SetSolutionStepValue(y, z), 1)
        populate(cls.model_part.Nodes, lambda x, y, z: x.SetValue(y, z), 2)
        populate(cls.model_part.Conditions, lambda x, y, z: x.SetValue(y, z), 3)
        populate(cls.model_part.Elements, lambda x, y, z: x.SetValue(y, z), 4)

        cls.data_location = Kratos.Globals.DataLocation

        cls.data_locations = [Kratos.Globals.DataLocation.NodeHistorical, Kratos.Globals.DataLocation.NodeNonHistorical, Kratos.Globals.DataLocation.Condition, Kratos.Globals.DataLocation.Element]
        cls.norms_dict = {
            Kratos.PRESSURE: [SpatialMethodTests.ValueNorm(), KratosStats.Norms.L2(), KratosStats.Norms.Infinity()],
            Kratos.VELOCITY: [SpatialMethodTests.ValueNorm(), KratosStats.Norms.L2(), KratosStats.Norms.Infinity(), KratosStats.Norms.P(2.5)],
            Kratos.INITIAL_STRAIN: [SpatialMethodTests.ValueNorm(), KratosStats.Norms.L2(), KratosStats.Norms.Infinity(), KratosStats.Norms.P(2.5)],
            Kratos.LOCAL_TANGENT_MATRIX: [SpatialMethodTests.ValueNorm(), KratosStats.Norms.L2(), KratosStats.Norms.Infinity(), KratosStats.Norms.P(2.5), KratosStats.Norms.Trace(), KratosStats.Norms.LPQ(1.3, 2.3)]
        }

    def __GetDataRetrievalInfo(self, data_location: Kratos.Globals.DataLocation) -> Any:
        if data_location == Kratos.Globals.DataLocation.NodeHistorical:
            return self.model_part.GetCommunicator().LocalMesh().Nodes, lambda x, y: x.GetSolutionStepValue(y)
        elif data_location == Kratos.Globals.DataLocation.NodeNonHistorical:
            return self.model_part.GetCommunicator().LocalMesh().Nodes,lambda x, y: x.GetValue(y)
        elif data_location == Kratos.Globals.DataLocation.Condition:
            return self.model_part.GetCommunicator().LocalMesh().Conditions,lambda x, y: x.GetValue(y)
        elif data_location == Kratos.Globals.DataLocation.Element:
            return self.model_part.GetCommunicator().LocalMesh().Elements,lambda x, y: x.GetValue(y)
        else:
            raise RuntimeError("Unsupported data location provided.")

    @staticmethod
    def __ConvertValueToList(value: Any) -> 'list[float]':
        if isinstance(value, float) or isinstance(value, int):
            return [value]
        elif isinstance(value, Kratos.Array3):
            return [value[0], value[1], value[2]]
        elif isinstance(value, Kratos.Vector):
            return [value[0], value[1], value[2], value[3], value[4]]
        elif isinstance(value, Kratos.Matrix):
            return [value[0, 0], value[0, 1], value[1, 0], value[1, 1]]
        else:
            raise RuntimeError("Unsupported value type.")

    @staticmethod
    def __ConvertListToValue(values_list: 'list[float]') -> Any:
        if  len(values_list) == 1:
            return values_list[0]
        elif len(values_list) == 3:
            return Kratos.Array3([values_list[0], values_list[1], values_list[2]])
        elif len(values_list) == 5:
            return Kratos.Vector([values_list[0], values_list[1], values_list[2], values_list[3], values_list[4]])
        elif len(values_list) == 4:
            return Kratos.Matrix([[values_list[0], values_list[1]], [values_list[2], values_list[3]]])
        else:
            raise RuntimeError("Unsupported values list length.")

    @staticmethod
    def __GetInitializedValue(variable: Any, norm: Any, initialization_value):
        if isinstance(variable, Kratos.DoubleVariable):
            initialized_variable_value = initialization_value
        elif isinstance(variable, Kratos.Array1DVariable3):
            initialized_variable_value = Kratos.Array3([initialization_value, initialization_value, initialization_value])
        elif isinstance(variable, Kratos.VectorVariable):
            initialized_variable_value = Kratos.Vector([initialization_value, initialization_value, initialization_value, initialization_value, initialization_value])
        elif isinstance(variable, Kratos.MatrixVariable):
            initialized_variable_value = Kratos.Matrix([[initialization_value, initialization_value], [initialization_value, initialization_value]])
        else:
            raise RuntimeError("Unsupported variable type.")

        if isinstance(norm, SpatialMethodTests.ValueNorm):
            return initialized_variable_value
        else:
            return initialization_value

    def __RunTest(self, stat_method: Any, analytical_method: Any) -> None:
        def run_granular_test(stat_method: Any, variable: Any, data_location: Kratos.Globals.DataLocation, ref_value: Any, norm: Any = None):
            if isinstance(norm, SpatialMethodTests.ValueNorm):
                stat_results = stat_method(self.model_part, variable, data_location)
            else:
                stat_results = stat_method(self.model_part, variable, data_location, norm)

            CheckValues(self, stat_results, ref_value, 12)

        for data_location in self.data_locations:
            container, data_retrieval_method = self.__GetDataRetrievalInfo(data_location)
            for variable, norms_list in self.norms_dict.items():
                for norm in norms_list:
                    ref_value = analytical_method(container, variable, norm, data_retrieval_method)
                    run_granular_test(stat_method, variable, data_location, ref_value, norm)

    def testSumMethod(self):
        def analytical_method(container, variable, norm, value_retriever):
            v = SpatialMethodTests.__GetInitializedValue(variable, norm, 0.0)
            for item in container:
                v += norm.Evaluate(value_retriever(item, variable))
            v_list = SpatialMethodTests.__ConvertValueToList(v)
            global_v_list = []
            for v in v_list:
                global_v_list.append(self.model_part.GetCommunicator().GetDataCommunicator().SumAll(v))
            return SpatialMethodTests.__ConvertListToValue(global_v_list)

        self.__RunTest(KratosStats.SpatialMethods.Sum, analytical_method)

    def testMeanMethod(self):
        def analytical_method(container, variable, norm, value_retriever):
            v = SpatialMethodTests.__GetInitializedValue(variable, norm, 0.0)
            for item in container:
                v += norm.Evaluate(value_retriever(item, variable))
            v_list = SpatialMethodTests.__ConvertValueToList(v)
            global_v_list = []
            for v in v_list:
                global_v_list.append(self.model_part.GetCommunicator().GetDataCommunicator().SumAll(v))
            global_sum = SpatialMethodTests.__ConvertListToValue(global_v_list)

            global_n = self.model_part.GetCommunicator().GetDataCommunicator().SumAll(len(container))
            return global_sum / global_n

        self.__RunTest(KratosStats.SpatialMethods.Mean, analytical_method)

    def testRootMeanSquareMethod(self):
        def analytical_method(container, variable, norm, value_retriever):
            global_n = self.model_part.GetCommunicator().GetDataCommunicator().SumAll(len(container))

            v_list = self.__ConvertValueToList(SpatialMethodTests.__GetInitializedValue(variable, norm, 0.0))
            for item in container:
                temp = self.__ConvertValueToList(norm.Evaluate(value_retriever(item, variable)))
                for i, t_i in enumerate(temp):
                    v_list[i] += t_i ** 2
            global_v_list = []
            for v in v_list:
                global_v_list.append(sqrt(self.model_part.GetCommunicator().GetDataCommunicator().SumAll(v) / global_n))
            return SpatialMethodTests.__ConvertListToValue(global_v_list)

        self.__RunTest(KratosStats.SpatialMethods.RootMeanSquare, analytical_method)

    def testVarianceMethod(self):
        def analytical_method(container, variable, norm, value_retriever):
            global_n = self.model_part.GetCommunicator().GetDataCommunicator().SumAll(len(container))

            v_mean_list = self.__ConvertValueToList(SpatialMethodTests.__GetInitializedValue(variable, norm, 0.0))
            v_variance_list = self.__ConvertValueToList(SpatialMethodTests.__GetInitializedValue(variable, norm, 0.0))
            for item in container:
                temp = self.__ConvertValueToList(norm.Evaluate(value_retriever(item, variable)))
                for i, t_i in enumerate(temp):
                    v_mean_list[i] += t_i
                    v_variance_list[i] += t_i ** 2
            global_v_mean_list = []
            global_v_variance_list = []
            for v_mean, v_variance in zip(v_mean_list, v_variance_list):
                current_mean = self.model_part.GetCommunicator().GetDataCommunicator().SumAll(v_mean) / global_n
                global_v_mean_list.append(current_mean)
                global_v_variance_list.append(self.model_part.GetCommunicator().GetDataCommunicator().SumAll(v_variance) / global_n - current_mean ** 2)
            return SpatialMethodTests.__ConvertListToValue(global_v_mean_list), SpatialMethodTests.__ConvertListToValue(global_v_variance_list)

        self.__RunTest(KratosStats.SpatialMethods.Variance, analytical_method)

    def testMinMethod(self):
        data_communicator: Kratos.DataCommunicator = self.model_part.GetCommunicator().GetDataCommunicator()
        def analytical_method(container, variable, norm, value_retriever):
            v_list = SpatialMethodTests.__ConvertValueToList(SpatialMethodTests.__GetInitializedValue(variable, norm, 1e+16))
            v_id_list = [1e+5] * len(v_list)
            for item in container:
                current_v = SpatialMethodTests.__ConvertValueToList(norm.Evaluate(value_retriever(item, variable)))
                for i, c_v in enumerate(current_v):
                    if v_list[i] > c_v:
                        v_list[i] = c_v
                        v_id_list[i] = item.Id
                    elif v_list[i] == c_v:
                        v_id_list[i] = min(v_id_list[i], item.Id)

            all_v_list = data_communicator.AllGathervDoubles(v_list)
            all_v_id_list = data_communicator.AllGathervInts(v_id_list)

            for global_v_list, global_v_id_list in zip(all_v_list, all_v_id_list):
                for i, global_v in enumerate(global_v_list):
                    if v_list[i] > global_v:
                        v_list[i] = global_v
                        v_id_list[i] = global_v_id_list[i]
                    elif v_list[i] == c_v:
                        v_id_list[i] = min(v_id_list[i], global_v_id_list[i])

            if len(v_id_list) == 1:
                return SpatialMethodTests.__ConvertListToValue(v_list), v_id_list[0]
            else:
                return SpatialMethodTests.__ConvertListToValue(v_list), v_id_list

        self.__RunTest(KratosStats.SpatialMethods.Min, analytical_method)

    def testMaxMethod(self):
        data_communicator: Kratos.DataCommunicator = self.model_part.GetCommunicator().GetDataCommunicator()
        def analytical_method(container, variable, norm, value_retriever):
            v_list = SpatialMethodTests.__ConvertValueToList(SpatialMethodTests.__GetInitializedValue(variable, norm, -1e+16))
            v_id_list = [1e+5] * len(v_list)
            for item in container:
                current_v = SpatialMethodTests.__ConvertValueToList(norm.Evaluate(value_retriever(item, variable)))
                for i, c_v in enumerate(current_v):
                    if v_list[i] < c_v:
                        v_list[i] = c_v
                        v_id_list[i] = item.Id
                    elif v_list[i] == c_v:
                        v_id_list[i] = min(v_id_list[i], item.Id)

            all_v_list = data_communicator.AllGathervDoubles(v_list)
            all_v_id_list = data_communicator.AllGathervInts(v_id_list)

            for global_v_list, global_v_id_list in zip(all_v_list, all_v_id_list):
                for i, global_v in enumerate(global_v_list):
                    if v_list[i] < global_v:
                        v_list[i] = global_v
                        v_id_list[i] = global_v_id_list[i]
                    elif v_list[i] == c_v:
                        v_id_list[i] = min(v_id_list[i], global_v_id_list[i])

            if len(v_id_list) == 1:
                return SpatialMethodTests.__ConvertListToValue(v_list), v_id_list[0]
            else:
                return SpatialMethodTests.__ConvertListToValue(v_list), v_id_list

        self.__RunTest(KratosStats.SpatialMethods.Max, analytical_method)

if __name__ == '__main__':
    KratosUnittest.main()