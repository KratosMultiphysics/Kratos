import numpy
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics as Kratos

class TestGeometryContainerDomainSize(KratosUnittest.TestCase):
    """
    Tests the DomainSize() method exposed on the Python interface of
    ModelPart::GeometryContainerType, i.e. the array of geometries stored
    in a ModelPart.

    The method computes the domain size of every geometry in the container
    (length, area or volume depending on the geometry dimension) and returns
    a numpy array with one entry per geometry, following the container
    iteration order.
    """

    def test_domain_size_mixed_geometries(self):
        current_model = Kratos.Model()
        model_part = current_model.CreateModelPart("Main")

        node_a = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        node_b = model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        node_c = model_part.CreateNewNode(3, 0.0, 1.0, 0.0)
        node_d = model_part.CreateNewNode(4, 0.0, 0.0, 1.0)
        node_e = model_part.CreateNewNode(5, 1.0, 1.0, 0.0)
        node_f = model_part.CreateNewNode(6, 1.0, 0.0, 1.0)
        node_g = model_part.CreateNewNode(7, 1.0, 1.0, 1.0)
        node_h = model_part.CreateNewNode(8, 0.0, 1.0, 1.0)

        # Line: length 1, triangle: area 0.5, tetra: volume 1/6, hexa: volume 1
        line = Kratos.Line2D2(node_a, node_b)
        triangle = Kratos.Triangle2D3(node_a, node_b, node_c)
        tetra = Kratos.Tetrahedra3D4(node_a, node_b, node_c, node_d)
        hexa = Kratos.Hexahedra3D8(node_a, node_b, node_e, node_c, node_d, node_f, node_g, node_h)

        for geometry in (line, triangle, tetra, hexa):
            model_part.AddGeometry(geometry)

        expected_size = {line.Id: 1.0, triangle.Id: 0.5, tetra.Id: 1.0 / 6.0, hexa.Id: 1.0}

        domain_sizes = model_part.Geometries.DomainSize()

        self.assertIsInstance(domain_sizes, numpy.ndarray)
        self.assertEqual(domain_sizes.shape, (4,))
        for index, geometry in enumerate(model_part.Geometries):
            self.assertAlmostEqual(domain_sizes[index], expected_size[geometry.Id], places=13)
            self.assertAlmostEqual(domain_sizes[index], geometry.DomainSize(), places=13)

    def test_domain_size_empty_container(self):
        current_model = Kratos.Model()
        model_part = current_model.CreateModelPart("Main")

        domain_sizes = model_part.Geometries.DomainSize()

        self.assertIsInstance(domain_sizes, numpy.ndarray)
        self.assertEqual(domain_sizes.shape, (0,))

if __name__ == '__main__':
    KratosUnittest.main()
