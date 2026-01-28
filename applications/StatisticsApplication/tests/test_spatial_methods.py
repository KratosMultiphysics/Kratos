from typing import Any
from math import sqrt

import KratosMultiphysics as Kratos

from KratosMultiphysics import IsDistributedRun
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

    def __GetContainerExpression(self, variable, data_location):
        if data_location == self.data_location.NodeHistorical:
            container_expression = Kratos.Expression.NodalExpression(self.model_part)
            Kratos.Expression.VariableExpressionIO.Read(container_expression, variable, True)
        elif data_location == self.data_location.NodeNonHistorical:
            container_expression = Kratos.Expression.NodalExpression(self.model_part)
            Kratos.Expression.VariableExpressionIO.Read(container_expression, variable, False)
        elif data_location == self.data_location.Condition:
            container_expression = Kratos.Expression.ConditionExpression(self.model_part)
            Kratos.Expression.VariableExpressionIO.Read(container_expression, variable)
        elif data_location == self.data_location.Element:
            container_expression = Kratos.Expression.ElementExpression(self.model_part)
            Kratos.Expression.VariableExpressionIO.Read(container_expression, variable)
        else:
            raise RuntimeError("Unsupported data container type.")

        return container_expression

    def __RunVariableTest(self, stat_method: Any, ref_values: Any, *args) -> None:
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

    def __RunContainerExpressionTest(self, stat_method, *args):
        for data_location in self.data_locations:
            for variable, norms_list in self.norms_dict.items():
                container_expression = self.__GetContainerExpression(variable, data_location)
                for norm in norms_list:
                    if isinstance(norm, SpatialMethodTests.ValueNorm):
                        exp_value = stat_method(container_expression, *args)
                        var_value = stat_method(self.model_part, variable, data_location, *args)
                    else:
                        exp_value = stat_method(container_expression, *args, norm)
                        var_value = stat_method(self.model_part, variable, data_location, *args, norm)

                    if stat_method == KratosStats.SpatialMethods.Distribution:
                        modified_exp_result = (
                                exp_value.GetMin(),
                                exp_value.GetMax(),
                                exp_value.GetGroupUpperValues(),
                                exp_value.GetGroupNumberOfValues(),
                                exp_value.GetGroupValueDistributionPercentage(),
                                exp_value.GetGroupMeans(),
                                exp_value.GetGroupVariances()
                            )
                        modified_var_result = (
                                var_value.GetMin(),
                                var_value.GetMax(),
                                var_value.GetGroupUpperValues(),
                                var_value.GetGroupNumberOfValues(),
                                var_value.GetGroupValueDistributionPercentage(),
                                var_value.GetGroupMeans(),
                                var_value.GetGroupVariances()
                            )
                        CheckValues(self, modified_exp_result, modified_var_result, 12)
                    elif stat_method in [KratosStats.SpatialMethods.Min, KratosStats.SpatialMethods.Max, KratosStats.SpatialMethods.Median]:
                        # since expressions start indices from zero and model part entities indices start from one.

                        CheckValues(self, exp_value[0], var_value[0], 12)
                        if not IsDistributedRun():
                            # ids can only be checked in the shared memory case because,
                            # in the case of expressions, the corresponding index israndom.
                            if isinstance(exp_value[1], int):
                                CheckValues(self, exp_value[1], var_value[1] - 1, 12)
                            else:
                                for i, v in enumerate(exp_value[1]):
                                    CheckValues(self, v, var_value[1][i]-1, 12)
                    else:
                        CheckValues(self, exp_value, var_value, 12)

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

        self.__RunVariableTest(KratosStats.SpatialMethods.Sum, ref_values)

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

        self.__RunVariableTest(KratosStats.SpatialMethods.Mean, ref_values)

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

        self.__RunVariableTest(KratosStats.SpatialMethods.RootMeanSquare, ref_values)

    def testVarianceMethod(self):
        ref_values = {}
        for data_location in self.data_locations:
            ref_values[data_location] = {}
            for variable, norms_list in self.norms_dict.items():
                ref_values[data_location][variable] = {}
                for norm in norms_list:
                    if isinstance(norm, SpatialMethodTests.ValueNorm):
                        mean = KratosStats.SpatialMethods.Mean(self.model_part, variable, data_location)
                        rms = KratosStats.SpatialMethods.RootMeanSquare(self.model_part, variable, data_location)
                    else:
                        mean = KratosStats.SpatialMethods.Mean(self.model_part, variable, data_location, norm)
                        rms = KratosStats.SpatialMethods.RootMeanSquare(self.model_part, variable, data_location, norm)
                    ref_values[data_location][variable][norm] = (mean, KratosStats.MethodUtilities.RaiseToPower(rms, 2) - KratosStats.MethodUtilities.RaiseToPower(mean, 2))

        self.__RunVariableTest(KratosStats.SpatialMethods.Variance, ref_values)

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

        self.__RunVariableTest(KratosStats.SpatialMethods.Min, ref_values)

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

        self.__RunVariableTest(KratosStats.SpatialMethods.Max, ref_values)

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

        self.__RunVariableTest(KratosStats.SpatialMethods.Median, ref_values)

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
        self.__RunVariableTest(KratosStats.SpatialMethods.Distribution, ref_values, params)

    def testSumContainerExpression(self):
        self.__RunContainerExpressionTest(KratosStats.SpatialMethods.Sum)

    def testMeanContainerExpression(self):
        self.__RunContainerExpressionTest(KratosStats.SpatialMethods.Mean)

    def testRootMeanSquareContainerExpression(self):
        self.__RunContainerExpressionTest(KratosStats.SpatialMethods.RootMeanSquare)

    def testVarianceContainerExpression(self):
        self.__RunContainerExpressionTest(KratosStats.SpatialMethods.Variance)

    def testMinExpression(self):
        self.__RunContainerExpressionTest(KratosStats.SpatialMethods.Min)

    def testMaxExpression(self):
        self.__RunContainerExpressionTest(KratosStats.SpatialMethods.Max)

    def testMedianExpression(self):
        self.__RunContainerExpressionTest(KratosStats.SpatialMethods.Median)

    def testDistributionExpression(self):
        params = Kratos.Parameters("""{
            "number_of_value_groups" : 4,
            "min_value"              : "min",
            "max_value"              : "max"
        }""")
        self.__RunContainerExpressionTest(KratosStats.SpatialMethods.Distribution, params)

if __name__ == '__main__':
    KratosUnittest.main()