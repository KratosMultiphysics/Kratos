import KratosMultiphysics as Kratos

import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.StatisticsApplication.spatial_utilities import GetItemContainer
from KratosMultiphysics.StatisticsApplication.method_utilities import GetNormTypeContainer
from KratosMultiphysics.StatisticsApplication.method_utilities import GetMethod

from random import uniform


class SpatialMethodTests(KratosUnittest.TestCase):
    def setUp(self):
        self.model = Kratos.Model()
        self.model_part = self.model.CreateModelPart("test_model_part")
        self.containers_to_test = [
            "nodal_historical", "nodal_non_historical",
            "element_non_historical", "condition_non_historical"
        ]

        self.test_cases = {}
        self.test_cases[Kratos.DENSITY] = ["none", "magnitude", "value"]
        self.test_cases[Kratos.VELOCITY] = [
            "none", "magnitude", "component_x", "component_y", "component_z"
        ]
        self.test_cases[Kratos.LOAD_MESHES] = [
            "magnitude", "index_0", "index_1", "index_2", "index_4"
        ]
        self.test_cases[Kratos.GREEN_LAGRANGE_STRAIN_TENSOR] = [
            "frobenius", "index_(0,0)", "index_(0,1)", "index_(4,1)",
            "index_(1,1)"
        ]

        self.norm_only_methods = ["min", "max", "median", "distribution"]

        SpatialMethodTests.__AddNodalSolutionStepVariables(self.model_part)
        SpatialMethodTests.__CreateModelPart(self.model_part)
        SpatialMethodTests.__InitializeVariables(self.model_part)

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def testSumMethod(self):
        def analytical_method(container, container_type, norm_type, variable):
            analytical_value = SpatialMethodTests.__GetInitialValue(
                variable, norm_type)
            for item in container:
                analytical_value += SpatialMethodTests.__GetNormValue(
                    variable,
                    SpatialMethodTests.__GetValue(item, container_type,
                                                  variable), norm_type)
            return analytical_value

        self.__TestMethod("sum", analytical_method)

    def testMeanMethod(self):
        def analytical_method(container, container_type, norm_type, variable):
            analytical_value = SpatialMethodTests.__GetInitialValue(
                variable, norm_type)
            for item in container:
                analytical_value += SpatialMethodTests.__GetNormValue(
                    variable,
                    SpatialMethodTests.__GetValue(item, container_type,
                                                  variable), norm_type)

            return analytical_value / len(container)

        self.__TestMethod("mean", analytical_method)

    def testVarianceMethod(self):
        def analytical_method(container, container_type, norm_type, variable):
            mean_value = SpatialMethodTests.__GetInitialValue(
                variable, norm_type)
            variance_value = SpatialMethodTests.__GetInitialValue(
                variable, norm_type)
            for item in container:
                current_value = SpatialMethodTests.__GetNormValue(
                    variable,
                    SpatialMethodTests.__GetValue(item, container_type,
                                                  variable), norm_type)
                mean_value += current_value
                variance_value += KratosStats.MethodUtilities.RaiseToPower(
                    current_value, 2)

            n = len(container)
            mean_value /= n
            variance_value = variance_value / n - KratosStats.MethodUtilities.RaiseToPower(
                mean_value, 2)

            return mean_value, variance_value

        self.__TestMethod("variance", analytical_method)

    def testMinMethod(self):
        def analytical_method(container, container_type, norm_type, variable):
            analytical_value = 1e+12
            for item in container:
                current_value = SpatialMethodTests.__GetNormValue(
                    variable,
                    SpatialMethodTests.__GetValue(item, container_type,
                                                  variable), norm_type)
                if (current_value < analytical_value):
                    analytical_value = current_value
                    analytical_id = item.Id

            return analytical_value, analytical_id

        self.__TestMethod("min", analytical_method)

    def testMaxMethod(self):
        def analytical_method(container, container_type, norm_type, variable):
            analytical_value = -1e+12
            for item in container:
                current_value = SpatialMethodTests.__GetNormValue(
                    variable,
                    SpatialMethodTests.__GetValue(item, container_type,
                                                  variable), norm_type)
                if (current_value > analytical_value):
                    analytical_value = current_value
                    analytical_id = item.Id

            return analytical_value, analytical_id

        self.__TestMethod("max", analytical_method)

    def testMedianMethod(self):
        def analytical_method(container, container_type, norm_type, variable):
            item_values = []
            for item in container:
                current_value = SpatialMethodTests.__GetNormValue(
                    variable,
                    SpatialMethodTests.__GetValue(item, container_type,
                                                  variable), norm_type)

                item_values.append(current_value)

            item_values = sorted(item_values)
            n = len(item_values)
            if (n % 2 != 0):
                return item_values[n // 2]
            else:
                return (item_values[(n - 1) // 2] + item_values[n // 2]) * 0.5

        self.__TestMethod("median", analytical_method)

    def testDistributionMethod(self):
        default_parameters = Kratos.Parameters("""
        {
            "number_of_value_groups" : 10,
            "min_value"              : "min",
            "max_value"              : "max"
        }""")

        def analytical_method(container, container_type, norm_type, variable):
            item_values = []
            for item in container:
                current_value = SpatialMethodTests.__GetNormValue(
                    variable,
                    SpatialMethodTests.__GetValue(item, container_type,
                                                  variable), norm_type)
                item_values.append(current_value)

            min_value = min(item_values)
            max_value = max(item_values)
            group_limits = [
                min_value + (max_value - min_value) * i / 10 for i in range(11)
            ]
            group_limits[-1] += 1.0
            group_limits.append(1e+100)

            data_distribution = [0 for i in range(len(group_limits))]
            for value in item_values:
                for i in range(len(group_limits)):
                    if (value < group_limits[i]):
                        data_distribution[i] += 1
                        break
            percentage_data_distribution = []
            for i in range(len(group_limits)):
                percentage_data_distribution.append(data_distribution[i] /
                                                    len(item_values))

            group_limits[-2] -= 1.0
            group_limits[-1] = max_value
            return min_value, max_value, group_limits, data_distribution, percentage_data_distribution

        self.__TestMethod("distribution", analytical_method,
                          default_parameters)

    def __TestMethod(self,
                     test_method_name,
                     analytical_method,
                     method_params=Kratos.Parameters("""{}""")):
        for container_type in self.containers_to_test:
            container = SpatialMethodTests.__GetContainer(
                self.model_part, container_type)
            item_method_container = GetItemContainer(container_type)
            for variable, norm_types in self.test_cases.items():
                for norm_type in norm_types:
                    item_method_norm_container = GetNormTypeContainer(
                        item_method_container, norm_type)

                    if (norm_type == "none"
                            and test_method_name in self.norm_only_methods):
                        continue

                    test_method = GetMethod(item_method_norm_container,
                                            test_method_name)
                    if (norm_type == "none"):
                        method_value = test_method(self.model_part, variable)
                    else:
                        method_value = test_method(norm_type, self.model_part,
                                                   variable, method_params)

                    analytical_value = analytical_method(
                        container, container_type, norm_type, variable)
                    self.__CheckValues(analytical_value, method_value, 10)

    def __CheckValues(self, value_a, value_b, tolerance):
        if (isinstance(value_a, tuple)):
            for i in range(len(value_a)):
                self.__CheckValues(value_a[i], value_b[i], tolerance)
        else:
            if (isinstance(value_a, Kratos.Matrix)):
                self.assertMatrixAlmostEqual(value_a, value_b, tolerance)
            elif (isinstance(value_a, Kratos.Vector)):
                self.assertVectorAlmostEqual(value_a, value_b, tolerance)
            elif (isinstance(value_a, list)):
                for i in range(len(value_a)):
                    self.assertAlmostEqual(value_a[i], value_b[i], tolerance)
            else:
                self.assertAlmostEqual(value_a, value_b, tolerance)

    @staticmethod
    def __AddNodalSolutionStepVariables(model_part):
        model_part.AddNodalSolutionStepVariable(Kratos.DENSITY)
        model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        model_part.AddNodalSolutionStepVariable(Kratos.LOAD_MESHES)
        model_part.AddNodalSolutionStepVariable(
            Kratos.GREEN_LAGRANGE_STRAIN_TENSOR)

    @staticmethod
    def __CreateModelPart(model_part):
        model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        model_part.CreateNewNode(3, 2.0, 0.0, 0.0)
        model_part.CreateNewNode(4, 2.0, 1.0, 0.0)
        model_part.CreateNewNode(5, 1.5, 1.0, 0.0)
        model_part.CreateNewNode(6, 0.5, 1.0, 0.0)
        model_part.CreateNewNode(7, 0.0, 1.0, 0.0)

        prop = model_part.GetProperties()[0]

        model_part.CreateNewElement("Element2D3N", 1, [1, 6, 7], prop)
        model_part.CreateNewElement("Element2D3N", 2, [1, 2, 6], prop)
        model_part.CreateNewElement("Element2D3N", 3, [6, 2, 5], prop)
        model_part.CreateNewElement("Element2D3N", 4, [2, 3, 5], prop)
        model_part.CreateNewElement("Element2D3N", 5, [5, 3, 4], prop)

        model_part.CreateNewCondition("Condition2D2N", 1, [1, 2], prop)
        model_part.CreateNewCondition("Condition2D2N", 2, [2, 3], prop)
        model_part.CreateNewCondition("Condition2D2N", 3, [3, 4], prop)
        model_part.CreateNewCondition("Condition2D2N", 4, [4, 5], prop)
        model_part.CreateNewCondition("Condition2D2N", 5, [5, 6], prop)
        model_part.CreateNewCondition("Condition2D2N", 6, [6, 7], prop)
        model_part.CreateNewCondition("Condition2D2N", 7, [7, 1], prop)

    @staticmethod
    def __InitializeVariables(model_part):

        for node in model_part.Nodes:
            s, a, v, m = SpatialMethodTests.__GetNewValue()
            node.SetValue(Kratos.DENSITY, s)
            node.SetValue(Kratos.VELOCITY, a)
            node.SetValue(Kratos.LOAD_MESHES, v)
            node.SetValue(Kratos.GREEN_LAGRANGE_STRAIN_TENSOR, m)

            s, a, v, m = SpatialMethodTests.__GetNewValue()
            node.SetSolutionStepValue(Kratos.DENSITY, 0, s)
            node.SetSolutionStepValue(Kratos.VELOCITY, 0, a)
            node.SetSolutionStepValue(Kratos.LOAD_MESHES, 0, v)
            node.SetSolutionStepValue(Kratos.GREEN_LAGRANGE_STRAIN_TENSOR, 0,
                                      m)

        for condition in model_part.Conditions:
            s, a, v, m = SpatialMethodTests.__GetNewValue()
            condition.SetValue(Kratos.DENSITY, s)
            condition.SetValue(Kratos.VELOCITY, a)
            condition.SetValue(Kratos.LOAD_MESHES, v)
            condition.SetValue(Kratos.GREEN_LAGRANGE_STRAIN_TENSOR, m)

        for element in model_part.Elements:
            s, a, v, m = SpatialMethodTests.__GetNewValue()
            element.SetValue(Kratos.DENSITY, s)
            element.SetValue(Kratos.VELOCITY, a)
            element.SetValue(Kratos.LOAD_MESHES, v)
            element.SetValue(Kratos.GREEN_LAGRANGE_STRAIN_TENSOR, m)

    @staticmethod
    def __GetValue(item, container_type, variable):
        if (container_type.endswith("non_historical")):
            return item.GetValue(variable)
        else:
            return item.GetSolutionStepValue(variable)

    @staticmethod
    def __GetNormValue(variable, value, norm_type):
        if (norm_type == "none"):
            return value

        norm_method = KratosStats.MethodUtilities.GetNormMethod(
            variable, norm_type)
        return norm_method(value)

    @staticmethod
    def __GetContainer(model_part, container_type):
        if (container_type.startswith("nodal")):
            return model_part.Nodes
        elif (container_type.startswith("element")):
            return model_part.Elements
        elif (container_type.startswith("condition")):
            return model_part.Conditions

    @staticmethod
    def __GetNewValue():
        m = Kratos.Matrix(5, 5)
        v = Kratos.Vector(5)
        a = Kratos.Vector(3)

        for i in range(5):
            v[i] = uniform(-10.0, 10.0)
            a[i % 3] = uniform(-10.0, 10.0)
            for j in range(5):
                m[i, j] = uniform(-10.0, 10.0)

        return uniform(-10.0, 10.0), a, v, m

    @staticmethod
    def __GetInitialValue(variable, norm_type):
        if (norm_type == "none"):
            if (variable == Kratos.DENSITY):
                return 0.0
            elif (variable == Kratos.VELOCITY):
                return Kratos.Vector(3, 0.0)
            elif (variable == Kratos.LOAD_MESHES):
                return Kratos.Vector(5, 0.0)
            elif (variable == Kratos.GREEN_LAGRANGE_STRAIN_TENSOR):
                return Kratos.Matrix(5, 5, 0.0)
        else:
            return 0.0


if __name__ == '__main__':
    KratosUnittest.main()