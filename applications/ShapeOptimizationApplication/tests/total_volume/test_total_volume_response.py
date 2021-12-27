import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as UnitTest
from KratosMultiphysics.ShapeOptimizationApplication.response_functions.total_volume import TotalVolume


class TestTotalVolumeResponseFunction(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(3, 1.0, 1.0, 0.0)

        prop = cls.model_part.GetProperties()[0]
        cls.model_part.CreateNewElement("Element2D3N", 1, [2, 3, 1], prop)

        cls.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 2

        cls.response_function = TotalVolume(1, Kratos.Parameters("""{
            "response_type"            : "total_volume",
            "model_part_name"          : "test",
            "output_variable_name"     : "SHAPE_SENSITIVITY"
        }"""),  cls.model)

        cls.response_function.Initialize()
        cls.response_function.CalculateValue()
        cls.ref_value = cls.response_function.GetValue()
        cls.response_function.CalculateGradient()


        cls.element = cls.model_part.GetElement(1)

    def testShapeSensitivity2D(self):
        analytical_gradient = self.response_function.GetNodalGradient(Kratos.SHAPE_SENSITIVITY)

        geometry = self.element.GetGeometry()
        delta = 1e-9
        for c, node in enumerate(geometry):
            node.X += delta
            self.response_function.CalculateValue()
            fd_shape_sensitivity = (self.response_function.GetValue() - self.ref_value) / delta
            self.assertAlmostEqual(fd_shape_sensitivity, analytical_gradient[node.Id][0], 6)
            node.X -= delta

            node.Y += delta
            self.response_function.CalculateValue()
            fd_shape_sensitivity = (self.response_function.GetValue() - self.ref_value) / delta
            self.assertAlmostEqual(fd_shape_sensitivity, analytical_gradient[node.Id][1], 6)
            node.Y -= delta

if __name__ == '__main__':
    UnitTest.main()