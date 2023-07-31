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
                setter_method(item, Kratos.LOCAL_TANGENT_MATRIX, Kratos.Matrix([[item_id + 1 + offset, -item_id - 2 - offset], [item_id + 3 + offset, item_id + 4 + offset]]))

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

    def __RunTest(self, stat_method: Any, ref_values: Any, *args) -> None:
        def run_granular_test(stat_method: Any, variable: Any, data_location: Kratos.Globals.DataLocation, ref_value: Any, norm: Any = None):
            if isinstance(norm, SpatialMethodTests.ValueNorm):
                stat_results = stat_method(self.model_part, variable, data_location, *args)
            else:
                stat_results = stat_method(self.model_part, variable, data_location, *args, norm)

            if stat_method == KratosStats.SpatialMethods.Distribution:
                stat_results: KratosStats.SpatialMethods.Array3DistributionInfo
                modified_stat_result = (
                        stat_results.GetMin(),
                        stat_results.GetMax(),
                        stat_results.GetGroupUpperValues(),
                        stat_results.GetGroupNumberOfValues(),
                        stat_results.GetGroupValueDistributionPercentage(),
                        stat_results.GetGroupMeans(),
                        stat_results.GetGroupVariances()
                    )
                CheckValues(self, modified_stat_result, ref_value, 12)
            else:
                CheckValues(self, stat_results, ref_value, 12)

        for data_location, info_1 in ref_values.items():
            for variable, info_2 in info_1.items():
                for norm, ref_value in info_2.items():
                    run_granular_test(stat_method, variable, data_location, ref_value, norm)

    def testSumMethod(self):
        ref_values = {
            self.data_location.NodeHistorical: {
                Kratos.PRESSURE: {
                    SpatialMethodTests.ValueNorm(): -63.0,
                    KratosStats.Norms.Infinity()  : 63
                },
                Kratos.VELOCITY: {
                    SpatialMethodTests.ValueNorm(): Kratos.Array3([63, -72, 81]),
                    KratosStats.Norms.Infinity()  : 81
                },
                Kratos.INITIAL_STRAIN: {
                    SpatialMethodTests.ValueNorm(): Kratos.Vector([63, -72, 81, -90, 99]),
                    KratosStats.Norms.Infinity()  : 99
                },
                Kratos.LOCAL_TANGENT_MATRIX: {
                    SpatialMethodTests.ValueNorm(): Kratos.Matrix([[63, -72], [81, 90]]),
                    KratosStats.Norms.Infinity()  : 171
                }
            }
        }

        self.__RunTest(KratosStats.SpatialMethods.Sum, ref_values)

    def testMeanMethod(self):
        ref_values = {
            self.data_location.NodeHistorical: {
                Kratos.PRESSURE: {
                    SpatialMethodTests.ValueNorm(): -63.0/9,
                    KratosStats.Norms.Infinity()  : 63/9
                },
                Kratos.VELOCITY: {
                    SpatialMethodTests.ValueNorm(): Kratos.Array3([63/9, -72/9, 81/9]),
                    KratosStats.Norms.Infinity()  : 81/9
                },
                Kratos.INITIAL_STRAIN: {
                    SpatialMethodTests.ValueNorm(): Kratos.Vector([63/9, -72/9, 81/9, -90/9, 99/9]),
                    KratosStats.Norms.Infinity()  : 99/9
                },
                Kratos.LOCAL_TANGENT_MATRIX: {
                    SpatialMethodTests.ValueNorm(): Kratos.Matrix([[63/9, -72/9], [81/9, 90/9]]),
                    KratosStats.Norms.Infinity()  : 171/9
                }
            }
        }

        self.__RunTest(KratosStats.SpatialMethods.Mean, ref_values)

    def testRootMeanSquareMethod(self):
        ref_values = {
            self.data_location.NodeHistorical: {
                Kratos.PRESSURE: {
                    SpatialMethodTests.ValueNorm(): sqrt((11*12*23/6-5)/9),
                    KratosStats.Norms.Infinity()  : sqrt((11*12*23/6-5)/9)
                },
                Kratos.VELOCITY: {
                    SpatialMethodTests.ValueNorm(): Kratos.Array3([sqrt((11*12*23/6-5)/9), sqrt((12*13*25/6-14)/9), sqrt((13*14*27/6-30)/9)]),
                    KratosStats.Norms.Infinity()  : sqrt((13*14*27/6-30)/9)
                },
                Kratos.INITIAL_STRAIN: {
                    SpatialMethodTests.ValueNorm(): Kratos.Vector([sqrt((11*12*23/6-5)/9), sqrt((12*13*25/6-14)/9), sqrt((13*14*27/6-30)/9), sqrt((14*15*29/6-55)/9), sqrt((15*16*31/6-91)/9)]),
                    KratosStats.Norms.Infinity()  : sqrt((15*16*31/6-91)/9)
                },
                Kratos.LOCAL_TANGENT_MATRIX: {
                    SpatialMethodTests.ValueNorm(): Kratos.Matrix([[sqrt((11*12*23/6-5)/9), sqrt((12*13*25/6-14)/9)], [sqrt((13*14*27/6-30)/9), sqrt((14*15*29/6-55)/9)]]),
                    KratosStats.Norms.Infinity()  : sqrt((14*29*27/3-5*11*9/3)/9)
                }
            }
        }

        self.__RunTest(KratosStats.SpatialMethods.RootMeanSquare, ref_values)

    def testVarianceMethod(self):
        data_locations = [Kratos.Globals.DataLocation.NodeHistorical, Kratos.Globals.DataLocation.NodeNonHistorical, Kratos.Globals.DataLocation.Condition, Kratos.Globals.DataLocation.Element]
        norms_dict = {
            Kratos.PRESSURE: [SpatialMethodTests.ValueNorm(), KratosStats.Norms.L2(), KratosStats.Norms.Infinity()],
            Kratos.VELOCITY: [SpatialMethodTests.ValueNorm(), KratosStats.Norms.L2(), KratosStats.Norms.Infinity(), KratosStats.Norms.P(2.5)],
            Kratos.INITIAL_STRAIN: [SpatialMethodTests.ValueNorm(), KratosStats.Norms.L2(), KratosStats.Norms.Infinity(), KratosStats.Norms.P(2.5)],
            Kratos.LOCAL_TANGENT_MATRIX: [SpatialMethodTests.ValueNorm(), KratosStats.Norms.L2(), KratosStats.Norms.Infinity(), KratosStats.Norms.P(2.5), KratosStats.Norms.Trace(), KratosStats.Norms.LPQ(1.3, 2.5)]
        }

        ref_values = {}
        for data_location in data_locations:
            ref_values[data_location] = {}
            for variable, norms_list in norms_dict.items():
                ref_values[data_location][variable] = {}
                for norm in norms_list:
                    if isinstance(norm, SpatialMethodTests.ValueNorm):
                        mean = KratosStats.SpatialMethods.Mean(self.model_part, variable, data_location)
                        rms = KratosStats.SpatialMethods.RootMeanSquare(self.model_part, variable, data_location)
                    else:
                        mean = KratosStats.SpatialMethods.Mean(self.model_part, variable, data_location, norm)
                        rms = KratosStats.SpatialMethods.RootMeanSquare(self.model_part, variable, data_location, norm)
                    ref_values[data_location][variable][norm] = (mean, KratosStats.MethodUtilities.RaiseToPower(rms, 2) - KratosStats.MethodUtilities.RaiseToPower(mean, 2))

        self.__RunTest(KratosStats.SpatialMethods.Variance, ref_values)

    def testMinMethod(self):
        ref_values = {
            self.data_location.NodeHistorical: {
                Kratos.PRESSURE: {
                    SpatialMethodTests.ValueNorm(): (-11, 9),
                    KratosStats.Norms.Infinity()  : (3, 1)
                },
                Kratos.VELOCITY: {
                    SpatialMethodTests.ValueNorm(): (Kratos.Array3([3, -12, 5]), [1, 9, 1]),
                    KratosStats.Norms.Infinity()  : (5, 1)
                },
                Kratos.INITIAL_STRAIN: {
                    SpatialMethodTests.ValueNorm(): (Kratos.Vector([3, -12, 5, -14, 7]), [1, 9, 1, 9, 1]),
                    KratosStats.Norms.Infinity()  : (7, 1)
                },
                Kratos.LOCAL_TANGENT_MATRIX: {
                    SpatialMethodTests.ValueNorm(): (Kratos.Matrix([[3, -12], [5, 6]]), [1, 9, 1, 1]),
                    KratosStats.Norms.Infinity()  : (11, 1)
                }
            }
        }

        self.__RunTest(KratosStats.SpatialMethods.Min, ref_values)

    def testMaxMethod(self):
        ref_values = {
            self.data_location.NodeHistorical: {
                Kratos.PRESSURE: {
                    SpatialMethodTests.ValueNorm(): (-3, 1),
                    KratosStats.Norms.Infinity()  : (11, 9)
                },
                Kratos.VELOCITY: {
                    SpatialMethodTests.ValueNorm(): (Kratos.Array3([11, -4, 13]), [9, 1, 9]),
                    KratosStats.Norms.Infinity()  : (13, 9)
                },
                Kratos.INITIAL_STRAIN: {
                    SpatialMethodTests.ValueNorm(): (Kratos.Vector([11, -4, 13, -6, 15]), [9, 1, 9, 1, 9]),
                    KratosStats.Norms.Infinity()  : (15, 9)
                },
                Kratos.LOCAL_TANGENT_MATRIX: {
                    SpatialMethodTests.ValueNorm(): (Kratos.Matrix([[11, -4], [13, 14]]), [9, 1, 9, 9]),
                    KratosStats.Norms.Infinity()  : (27, 9)
                }
            }
        }

        self.__RunTest(KratosStats.SpatialMethods.Max, ref_values)

    def testMedianMethod(self):
        ref_values = {
            self.data_location.NodeHistorical: {
                Kratos.PRESSURE: {
                    SpatialMethodTests.ValueNorm(): (-7, 5),
                    KratosStats.Norms.Infinity()  : (7, 5)
                },
                Kratos.VELOCITY: {
                    SpatialMethodTests.ValueNorm(): (Kratos.Array3([7, -8, 9]), [5, 5, 5]),
                    KratosStats.Norms.Infinity()  : (9, 5)
                },
                Kratos.INITIAL_STRAIN: {
                    SpatialMethodTests.ValueNorm(): (Kratos.Vector([7, -8, 9, -10, 11]), [5, 5, 5, 5, 5]),
                    KratosStats.Norms.Infinity()  : (11, 5)
                },
                Kratos.LOCAL_TANGENT_MATRIX: {
                    SpatialMethodTests.ValueNorm(): (Kratos.Matrix([[7, -8], [9, 10]]), [5, 5, 5, 5]),
                    KratosStats.Norms.Infinity()  : (19, 5)
                }
            }
        }

        self.__RunTest(KratosStats.SpatialMethods.Median, ref_values)

    def testDistributionMethod(self):
        ref_values = {
            self.data_location.NodeHistorical: {
                Kratos.PRESSURE: {
                    SpatialMethodTests.ValueNorm(): (
                        -11,
                        -3,
                        [(-11 + 8*i/4) for i in range(5)] + [-3.0],
                        [0,     2,    2,    2,    3, 0],
                        [0,   2/9,  2/9,  2/9,  3/9, 0],
                        [0, -10.5, -8.5, -6.5, -4.0, 0],
                        [0,  0.25, 0.25, 0.25,  2/3, 0]
                    ),
                    KratosStats.Norms.Infinity(): (
                        3,
                        11,
                        [(3 + 8*i/4) for i in range(5)] + [11],
                        [0,    2,    2,    2,    3, 0],
                        [0,  2/9,  2/9,  2/9,  3/9, 0],
                        [0,  3.5,  5.5,  7.5, 10.0, 0],
                        [0, 0.25, 0.25, 0.25,  2/3, 0]
                    ),
                },
                Kratos.VELOCITY: {
                    SpatialMethodTests.ValueNorm(): (
                        Kratos.Array3([ 3, -12,  5]),
                        Kratos.Array3([11,  -4, 13]),
                        [Kratos.Array3([(3 + 8*i/4), (-12 + 8*i/4), (5 + 8*i/4)]) for i in range(5)] + [Kratos.Array3([11, -4, 13])],
                        [               [0,0,0],                        [2,2,2],                       [2,2,2],                       [2,2,2],                    [3,3,3],               [0,0,0]],
                        [Kratos.Array3([0,0,0]),       Kratos.Array3([2,2,2])/9,      Kratos.Array3([2,2,2])/9,      Kratos.Array3([2,2,2])/9,   Kratos.Array3([3,3,3])/9, Kratos.Array3([0,0,0])],
                        [Kratos.Array3([0,0,0]), Kratos.Array3([3.5,-11.5,5.5]), Kratos.Array3([5.5,-9.5,7.5]), Kratos.Array3([7.5,-7.5,9.5]),  Kratos.Array3([10,-5,12]), Kratos.Array3([0,0,0])],
                        [Kratos.Array3([0,0,0]),    Kratos.Array3([1,1,1])*0.25,   Kratos.Array3([1,1,1])*0.25,   Kratos.Array3([1,1,1])*0.25, Kratos.Array3([1,1,1])*2/3, Kratos.Array3([0,0,0])],
                    ),
                    KratosStats.Norms.Infinity(): (
                        5,
                        13,
                        [(5 + 8*i/4) for i in range(5)] + [13],
                        [0,    2,    2,    2,    3, 0],
                        [0,  2/9,  2/9,  2/9,  3/9, 0],
                        [0,  5.5,  7.5,  9.5, 12.0, 0],
                        [0, 0.25, 0.25, 0.25,  2/3, 0]
                    ),
                },
                Kratos.INITIAL_STRAIN: {
                    SpatialMethodTests.ValueNorm(): (
                        Kratos.Vector([ 3, -12,  5, -14,  7]),
                        Kratos.Vector([11,  -4, 13,  -6, 15]),
                        [Kratos.Vector([(3 + 8*i/4), (-12 + 8*i/4), (5 + 8*i/4), (-14 + 8*i/4), (7 + 8*i/4)]) for i in range(5)] + [Kratos.Vector([11, -4, 13, -6, 15])],
                        [               [0,0,0,0,0],                              [2,2,2,2,2],                       [2,2,2,2,2],                       [2,2,2,2,2],                    [3,3,3,3,3],               [0,0,0,0,0]],
                        [Kratos.Vector([0,0,0,0,0]),             Kratos.Vector([2,2,2,2,2])/9,      Kratos.Vector([2,2,2,2,2])/9,      Kratos.Vector([2,2,2,2,2])/9,   Kratos.Vector([3,3,3,3,3])/9, Kratos.Vector([0,0,0,0,0])],
                        [Kratos.Vector([0,0,0,0,0]), Kratos.Vector([3.5,-11.5,5.5,-13.5,7.5]), Kratos.Vector([5.5,-9.5,7.5,-11.5,9.5]), Kratos.Vector([7.5,-7.5,9.5,-9.5,11.5]),  Kratos.Vector([10,-5,12,-7,14]), Kratos.Vector([0,0,0,0,0])],
                        [Kratos.Vector([0,0,0,0,0]),    Kratos.Vector([1,1,1,1,1])*0.25,   Kratos.Vector([1,1,1,1,1])*0.25,   Kratos.Vector([1,1,1,1,1])*0.25, Kratos.Vector([1,1,1,1,1])*2/3, Kratos.Vector([0,0,0,0,0])],
                    ),
                    KratosStats.Norms.Infinity(): (
                        7,
                        15,
                        [(7 + 8*i/4) for i in range(5)] + [15],
                        [0,    2,    2,    2,    3, 0],
                        [0,  2/9,  2/9,  2/9,  3/9, 0],
                        [0,  7.5,  9.5,  11.5, 14.0, 0],
                        [0, 0.25, 0.25, 0.25,  2/3, 0]
                    ),
                },
                Kratos.LOCAL_TANGENT_MATRIX: {
                    SpatialMethodTests.ValueNorm(): (
                        Kratos.Matrix([[ 3, -12],  [5,  6]]),
                        Kratos.Matrix([[11,  -4], [13, 14]]),
                        [Kratos.Matrix([[(3 + 8*i/4), (-12 + 8*i/4)], [(5 + 8*i/4), (6 + 8*i/4)]]) for i in range(5)] + [Kratos.Matrix([[11, -4], [13, 14]])],
                        [                   [0,0,0,0],                              [2,2,2,2],                             [2,2,2,2],                              [2,2,2,2],                        [3,3,3,3],                    [0,0,0,0]],
                        [Kratos.Matrix([[0,0],[0,0]]),         Kratos.Matrix([[2,2],[2,2]])/9,        Kratos.Matrix([[2,2],[2,2]])/9,         Kratos.Matrix([[2,2],[2,2]])/9,   Kratos.Matrix([[3,3],[3,3]])/9, Kratos.Matrix([[0,0],[0,0]])],
                        [Kratos.Matrix([[0,0],[0,0]]), Kratos.Matrix([[3.5,-11.5],[5.5,6.5]]), Kratos.Matrix([[5.5,-9.5],[7.5,8.5]]), Kratos.Matrix([[7.5,-7.5],[9.5,10.5]]), Kratos.Matrix([[10,-5],[12,13]]), Kratos.Matrix([[0,0],[0,0]])],
                        [Kratos.Matrix([[0,0],[0,0]]),      Kratos.Matrix([[1,1],[1,1]])*0.25,     Kratos.Matrix([[1,1],[1,1]])*0.25,      Kratos.Matrix([[1,1],[1,1]])*0.25, Kratos.Matrix([[1,1],[1,1]])*2/3, Kratos.Matrix([[0,0],[0,0]])],
                    ),
                    KratosStats.Norms.Infinity(): (
                        11,
                        27,
                        [(11 + 16*i/4) for i in range(5)] + [27],
                        [0,    2,    2,    2,    3, 0],
                        [0,  2/9,  2/9,  2/9,  3/9, 0],
                        [0,   12,   16,   20, 75/3, 0],
                        [0,  1.0,  1.0,  1.0,  8/3, 0]
                    ),
                }
            }
        }

        params = Kratos.Parameters("""{
            "number_of_value_groups" : 4,
            "min_value"              : "min",
            "max_value"              : "max"
        }""")
        self.__RunTest(KratosStats.SpatialMethods.Distribution, ref_values, params)

if __name__ == '__main__':
    KratosUnittest.main()