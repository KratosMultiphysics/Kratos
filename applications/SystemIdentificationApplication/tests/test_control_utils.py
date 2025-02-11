import numpy as np

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
from KratosMultiphysics.SystemIdentificationApplication.utilities.expression_utils import ExpressionUnionType

class TestControlUtils(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        """
          (0,0)
            1------2------3------4
            |      |      |      |
            |   1  |  2   |  3   |
            |      |      |      |
            5------6------7------8
            |      |      |      |
            |   4  |  5   |  6   |
            |      |      |      |
            9-----10-----11-----12
                                (6,4)
        """

        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.model_part.AddNodalSolutionStepVariable(Kratos.DISPLACEMENT)

        cls.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(2, 2.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(3, 4.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(4, 6.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(5, 0.0, 2.0, 0.0)
        cls.model_part.CreateNewNode(6, 2.0, 2.0, 0.0)
        cls.model_part.CreateNewNode(7, 4.0, 2.0, 0.0)
        cls.model_part.CreateNewNode(8, 6.0, 2.0, 0.0)
        cls.model_part.CreateNewNode(9, 0.0, 4.0, 0.0)
        cls.model_part.CreateNewNode(10, 2.0, 4.0, 0.0)
        cls.model_part.CreateNewNode(11, 4.0, 4.0, 0.0)
        cls.model_part.CreateNewNode(12, 6.0, 4.0, 0.0)

        prop = cls.model_part.CreateNewProperties(1)

        cls.model_part.CreateNewElement("Element2D4N", 1, [1, 2, 6, 5], prop).SetValue(Kratos.PRESSURE, 1.0)
        cls.model_part.CreateNewElement("Element2D4N", 2, [2, 3, 7, 6], prop).SetValue(Kratos.PRESSURE, 2.0)
        cls.model_part.CreateNewElement("Element2D4N", 3, [3, 4, 8, 7], prop).SetValue(Kratos.PRESSURE, 3.0)
        cls.model_part.CreateNewElement("Element2D4N", 4, [5, 6, 10, 9], prop).SetValue(Kratos.PRESSURE, 4.0)
        cls.model_part.CreateNewElement("Element2D4N", 5, [6, 7, 11, 10], prop).SetValue(Kratos.PRESSURE, 5.0)
        cls.model_part.CreateNewElement("Element2D4N", 6, [7, 8, 12, 11], prop).SetValue(Kratos.PRESSURE, 6.0)

        cls.list_of_points = KratosSI.ControlUtils.GetIntegrationPoints(cls.model_part.Elements)
        cls.list_of_integration_point_areas = KratosSI.ControlUtils.GetIntegrationPointAreas(cls.model_part.Elements)

        for node in cls.model_part.Nodes:
            node.SetValue(Kratos.DISTANCE, node.Id)

    def test_GetIntegrationPoints(self):
        ref_local_gauss_points = np.array([
            [0.42265, 0.42265, 0],
            [1.57735, 0.42265, 0],
            [1.57735, 1.57735, 0],
            [0.42265, 1.57735, 0]
        ])

        for i, element in enumerate(self.model_part.Elements):
            global_point = ref_local_gauss_points + element.GetGeometry().Center() - self.model_part.GetElement(1).GetGeometry().Center()
            for j in range(4):
                point = self.list_of_points[i * 4 + j]
                self.assertVectorAlmostEqual(global_point[j, :], [point.X, point.Y, point.Z], 4)

    def test_GetIntegrationPointAreas(self):
        self.assertEqual(len(self.list_of_points), len(self.list_of_integration_point_areas))
        for i, element in enumerate(self.model_part.Elements):
            self.assertAlmostEqual(self.list_of_integration_point_areas[i], element.GetGeometry().DomainSize() / 4)

    def test_EvaluateElementValuesAtPoints(self):
        values = KratosSI.ControlUtils.EvaluateElementValuesAtPoints(Kratos.PRESSURE, self.model_part, self.list_of_points)
        for i, element in enumerate(self.model_part.Elements):
            for j in range(4):
                self.assertEqual(element.GetValue(Kratos.PRESSURE), values[i * 4 + j])

    def test_EvaluateNodalNonHistoricalValuesAtPoints(self):
        reference_shape_func_values = np.array([
            (0.622008,0.166667,0.0446582,0.166667),
            (0.166667,0.622008,0.166667,0.0446582),
            (0.0446582,0.166667,0.622008,0.166667),
            (0.166667,0.0446582,0.166667,0.622008)
        ])

        values = KratosSI.ControlUtils.EvaluateNodalNonHistoricalValuesAtPoints(Kratos.DISTANCE, self.model_part, self.list_of_points)
        for i, element in enumerate(self.model_part.Elements):
            for j in range(4):
                p = 0.0
                for k, node in enumerate(element.GetGeometry()):
                    p += node.GetValue(Kratos.DISTANCE) * reference_shape_func_values[j, k]
                self.assertAlmostEqual(p, values[i * 4 + j], 5)

if __name__ == '__main__':
    UnitTest.main()