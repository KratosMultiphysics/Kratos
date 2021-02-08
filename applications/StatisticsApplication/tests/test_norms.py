import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.StatisticsApplication import MethodUtilities
from KratosMultiphysics.StatisticsApplication.test_utilities import GetRandomValue
from KratosMultiphysics.StatisticsApplication.test_utilities import GetRandomVector
from KratosMultiphysics.StatisticsApplication.test_utilities import GetRandomMatrix

from math import fabs


class NormTests(KratosUnittest.TestCase):
    def testScalarNorms(self):
        variable = Kratos.PRESSURE
        value = GetRandomValue()

        norm_type = "value"
        norm_method = MethodUtilities.GetNormMethod(variable, norm_type)
        method_value = norm_method(value)
        analytical_value = value
        self.assertAlmostEqual(method_value, analytical_value, 8)

        norm_type = "magnitude"
        norm_method = MethodUtilities.GetNormMethod(variable, norm_type)
        method_value = norm_method(value)
        analytical_value = fabs(value)
        self.assertAlmostEqual(method_value, analytical_value, 8)

    def testVector3DNorms(self):
        variable = Kratos.VELOCITY
        value = GetRandomVector(3)

        norm_type = "magnitude"
        norm_method = MethodUtilities.GetNormMethod(variable, norm_type)
        method_value = norm_method(value)
        analytical_value = pow(
            pow(value[0], 2) + pow(value[1], 2) + pow(value[2], 2), 0.5)
        self.assertAlmostEqual(method_value, analytical_value, 8)

        norm_type = "euclidean"
        norm_method = MethodUtilities.GetNormMethod(variable, norm_type)
        method_value = norm_method(value)
        analytical_value = pow(
            pow(value[0], 2) + pow(value[1], 2) + pow(value[2], 2), 0.5)
        self.assertAlmostEqual(method_value, analytical_value, 8)

        norm_type = "infinity"
        norm_method = MethodUtilities.GetNormMethod(variable, norm_type)
        method_value = norm_method(value)
        analytical_value = max(fabs(value[0]),
                               max(fabs(value[1]), fabs(value[2])))
        self.assertAlmostEqual(method_value, analytical_value, 8)

        norm_type = "pnorm_2.5"
        norm_method = MethodUtilities.GetNormMethod(variable, norm_type)
        method_value = norm_method(value)
        analytical_value = pow(
            pow(fabs(value[0]), 2.5) + pow(fabs(value[1]), 2.5) +
            pow(fabs(value[2]), 2.5), 1 / 2.5)
        self.assertAlmostEqual(method_value, analytical_value, 8)

        norm_type = "component_x"
        norm_method = MethodUtilities.GetNormMethod(variable, norm_type)
        method_value = norm_method(value)
        analytical_value = value[0]
        self.assertAlmostEqual(method_value, analytical_value, 8)

        norm_type = "component_y"
        norm_method = MethodUtilities.GetNormMethod(variable, norm_type)
        method_value = norm_method(value)
        analytical_value = value[1]
        self.assertAlmostEqual(method_value, analytical_value, 8)

        norm_type = "component_z"
        norm_method = MethodUtilities.GetNormMethod(variable, norm_type)
        method_value = norm_method(value)
        analytical_value = value[2]
        self.assertAlmostEqual(method_value, analytical_value, 8)

    def testVectorNorms(self):
        variable = Kratos.LOAD_MESHES
        value = GetRandomVector(4)

        norm_type = "magnitude"
        norm_method = MethodUtilities.GetNormMethod(variable, norm_type)
        method_value = norm_method(value)
        analytical_value = pow(
            pow(value[0], 2) + pow(value[1], 2) + pow(value[2], 2) +
            pow(value[3], 2), 0.5)
        self.assertAlmostEqual(method_value, analytical_value, 8)

        norm_type = "infinity"
        norm_method = MethodUtilities.GetNormMethod(variable, norm_type)
        method_value = norm_method(value)
        analytical_value = max(
            fabs(value[0]),
            max(fabs(value[1]), max(fabs(value[2]), fabs(value[3]))))
        self.assertAlmostEqual(method_value, analytical_value, 8)

        norm_type = "euclidean"
        norm_method = MethodUtilities.GetNormMethod(variable, norm_type)
        method_value = norm_method(value)
        analytical_value = pow(
            pow(value[0], 2) + pow(value[1], 2) + pow(value[2], 2) +
            pow(value[3], 2), 0.5)
        self.assertAlmostEqual(method_value, analytical_value, 8)

        norm_type = "pnorm_2.5"
        norm_method = MethodUtilities.GetNormMethod(variable, norm_type)
        method_value = norm_method(value)
        analytical_value = pow(
            pow(fabs(value[0]), 2.5) + pow(fabs(value[1]), 2.5) +
            pow(fabs(value[2]), 2.5) + pow(fabs(value[3]), 2.5), 1 / 2.5)
        self.assertAlmostEqual(method_value, analytical_value, 8)

        norm_type = "index_3"
        norm_method = MethodUtilities.GetNormMethod(variable, norm_type)
        method_value = norm_method(value)
        analytical_value = value[3]
        self.assertAlmostEqual(method_value, analytical_value, 8)

    def testMatrixNorms(self):
        variable = Kratos.GREEN_LAGRANGE_STRAIN_TENSOR
        value = GetRandomMatrix(2, 2)

        norm_type = "magnitude"
        norm_method = MethodUtilities.GetNormMethod(variable, norm_type)
        method_value = norm_method(value)
        analytical_value = pow(
            pow(value[0, 0], 2) + pow(value[0, 1], 2) + pow(value[1, 0], 2) +
            pow(value[1, 1], 2), 0.5)
        self.assertAlmostEqual(method_value, analytical_value, 8)

        norm_type = "frobenius"
        norm_method = MethodUtilities.GetNormMethod(variable, norm_type)
        method_value = norm_method(value)
        analytical_value = pow(
            pow(value[0, 0], 2) + pow(value[0, 1], 2) + pow(value[1, 0], 2) +
            pow(value[1, 1], 2), 0.5)
        self.assertAlmostEqual(method_value, analytical_value, 8)

        norm_type = "infinity"
        norm_method = MethodUtilities.GetNormMethod(variable, norm_type)
        method_value = norm_method(value)
        analytical_value = max(
            fabs(value[0, 0]) + fabs(value[0, 1]),
            fabs(value[1, 0]) + fabs(value[1, 1]))
        self.assertAlmostEqual(method_value, analytical_value, 8)

        norm_type = "trace"
        norm_method = MethodUtilities.GetNormMethod(variable, norm_type)
        method_value = norm_method(value)
        analytical_value = value[0, 0] + value[1, 1]
        self.assertAlmostEqual(method_value, analytical_value, 8)

        norm_type = "pnorm_2.5"
        norm_method = MethodUtilities.GetNormMethod(variable, norm_type)
        method_value = norm_method(value)
        analytical_value = pow(
            pow(fabs(value[0, 0]), 2.5) + pow(fabs(value[0, 1]), 2.5) +
            pow(fabs(value[1, 0]), 2.5) + pow(fabs(value[1, 1]), 2.5), 1 / 2.5)
        self.assertAlmostEqual(method_value, analytical_value, 8)

        norm_type = "index_(0,1)"
        norm_method = MethodUtilities.GetNormMethod(variable, norm_type)
        method_value = norm_method(value)
        analytical_value = value[0, 1]
        self.assertAlmostEqual(method_value, analytical_value, 8)

        norm_type = "lpqnorm_(2.5,2.1)"
        norm_method = MethodUtilities.GetNormMethod(variable, norm_type)
        method_value = norm_method(value)
        analytical_value = pow(
            pow(
                pow(fabs(value[0, 0]), 2.5) + pow(fabs(value[1, 0]), 2.5),
                2.1 / 2.5) + pow(
                    pow(fabs(value[0, 1]), 2.5) + pow(fabs(value[1, 1]), 2.5),
                    2.1 / 2.5), 1 / 2.1)
        self.assertAlmostEqual(method_value, analytical_value, 8)


if __name__ == '__main__':
    KratosUnittest.main()