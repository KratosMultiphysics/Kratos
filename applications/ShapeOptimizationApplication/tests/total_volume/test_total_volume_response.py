import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as UnitTest
from KratosMultiphysics.ShapeOptimizationApplication.response_functions.total_volume import TotalVolume


class TestTotalVolumeResponseFunction(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls, model, model_part):
        cls.response_function = TotalVolume(1, Kratos.Parameters("""{
            "response_type"            : "total_volume",
            "model_part_name"          : "test"
        }"""),  model)

        cls.response_function.Initialize()
        cls.response_function.CalculateValue()
        cls.ref_value = cls.response_function.GetValue()
        cls.response_function.CalculateGradient()
        cls.element = model_part.GetElement(1)

    def CalculateFiniteDifferenceSensitivity(self, node, direction_index, delta):
        if direction_index == 0:
            def nodal_perturbation_method(node, delta):
                node.X += delta
        elif direction_index == 1:
            def nodal_perturbation_method(node, delta):
                node.Y += delta
        elif direction_index == 2:
            def nodal_perturbation_method(node, delta):
                node.Z += delta
        else:
            raise Exception("Unsupported direction index")

        nodal_perturbation_method(node, delta)
        self.response_function.CalculateValue()
        fd_shape_sensitivity = (self.response_function.GetValue() - self.ref_value) / delta
        nodal_perturbation_method(node, -delta)
        return fd_shape_sensitivity

class TestTotalVolumeResponseFunction2D(TestTotalVolumeResponseFunction):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.model_part.AddNodalSolutionStepVariable(Kratos.SHAPE_SENSITIVITY)
        cls.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
        prop = cls.model_part.GetProperties()[0]
        cls.model_part.CreateNewElement("Element2D3N", 1, [2, 3, 1], prop)
        cls.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 2

        super().setUpClass(cls.model, cls.model_part)

    def testShapeSensitivity2D(self):
        analytical_gradient = self.response_function.GetNodalGradient(Kratos.SHAPE_SENSITIVITY)

        geometry = self.element.GetGeometry()
        delta = 1e-9
        for node in geometry:
            for k in range(3):
                fd_shape_sensitivity = self.CalculateFiniteDifferenceSensitivity(node, k, delta)
                self.assertAlmostEqual(fd_shape_sensitivity, analytical_gradient[node.Id][k], 6)

            self.assertAlmostEqual(fd_shape_sensitivity, 0.0, 6)

class TestTotalVolumeResponseFunction3D(TestTotalVolumeResponseFunction):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.model_part.AddNodalSolutionStepVariable(Kratos.SHAPE_SENSITIVITY)
        cls.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
        cls.model_part.CreateNewNode(4, 0.0, 0.0, 1.0)
        prop = cls.model_part.GetProperties()[0]
        cls.model_part.CreateNewElement("Element3D4N", 1, [2, 3, 1, 4], prop)
        cls.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3

        super().setUpClass(cls.model, cls.model_part)

    def testShapeSensitivity3D(self):
        analytical_gradient = self.response_function.GetNodalGradient(Kratos.SHAPE_SENSITIVITY)

        geometry = self.element.GetGeometry()
        delta = 1e-9
        for node in geometry:
            for k in range(3):
                fd_shape_sensitivity = self.CalculateFiniteDifferenceSensitivity(node, k, delta)
                self.assertAlmostEqual(fd_shape_sensitivity, analytical_gradient[node.Id][k], 6)

if __name__ == '__main__':
    UnitTest.main()