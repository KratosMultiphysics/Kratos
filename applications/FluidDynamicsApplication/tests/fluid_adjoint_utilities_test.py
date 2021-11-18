import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

class FluidAdjointUtilities():
    @staticmethod
    def PerturbNode(node, direction, delta):
        if direction == 0:
            node.X += delta
        elif direction == 1:
            node.Y += delta
        elif direction == 2:
            node.Z += delta
        else:
            raise Exception("Unsupported direction is requested at nodal perturbation.")

class FluidAdjointUtilities2D(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")

        cls.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(2, 0.5, 0.5, 0.0)
        cls.model_part.CreateNewNode(3, 1.5, 2.0, 0.0)

        prop = cls.model_part.GetProperties()[0]
        cls.element = cls.model_part.CreateNewElement("Element2D3N", 1, [2, 3, 1], prop)

    def testCalculateTriangleAreaDerivative(self):
        geometry = self.element.GetGeometry()

        ref_area = geometry.DomainSize()

        delta = 1e-6
        for c, node in enumerate(geometry):
            for k in range(2):
                # calculate analytical derivatives
                analytical_derivative = KratosCFD.FluidAdjointUtilities2D.CalculateTriangleAreaDerivative(geometry, c, k)

                FluidAdjointUtilities.PerturbNode(node, k, delta)
                fd_derivative = (geometry.DomainSize() - ref_area) / delta
                FluidAdjointUtilities.PerturbNode(node, k, -delta)

                self.assertAlmostEqual(analytical_derivative, fd_derivative, 8)

class FluidAdjointUtilities3D(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")

        cls.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(2, 0.5, 0.5, 1.0)
        cls.model_part.CreateNewNode(3, 1.5, 2.0, 4.0)

        prop = cls.model_part.GetProperties()[0]
        cls.condition = cls.model_part.CreateNewCondition("SurfaceCondition3D3N", 1, [2, 3, 1], prop)

    def testCalculateTriangleAreaDerivative(self):
        geometry = self.condition.GetGeometry()

        ref_area = geometry.DomainSize()

        delta = 1e-8
        for c, node in enumerate(geometry):
            for k in range(3):
                # calculate analytical derivatives
                analytical_derivative = KratosCFD.FluidAdjointUtilities3D.CalculateTriangleAreaDerivative(geometry, c, k)

                FluidAdjointUtilities.PerturbNode(node, k, delta)
                fd_derivative = (geometry.DomainSize() - ref_area) / delta
                FluidAdjointUtilities.PerturbNode(node, k, -delta)

                self.assertAlmostEqual(analytical_derivative, fd_derivative, 5)

if __name__ == '__main__':
    UnitTest.main()
