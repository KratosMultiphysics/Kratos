from __future__ import print_function, absolute_import, division

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics import *

class TestGeometry(KratosUnittest.TestCase):
    def setUp(self):
        pass

    def test_tetrahedra_3D4N(self):
        current_model = Model()
        model_part= current_model.CreateModelPart("Main")
        tester = GeometryTesterUtility()
        self.assertTrue( tester.TestTetrahedra3D4N(model_part) )

    def test_tetrahedra_2D3N(self):
        current_model = Model()
        model_part= current_model.CreateModelPart("Main")
        tester = GeometryTesterUtility()
        self.assertTrue( tester.TestTriangle2D3N(model_part) )

    def test_tetrahedra_2D6N(self):
        current_model = Model()
        model_part= current_model.CreateModelPart("Main")
        tester = GeometryTesterUtility()
        self.assertTrue( tester.TestTriangle2D6N(model_part) )

    def test_tetrahedra_3D10N(self):
        current_model = Model()
        model_part= current_model.CreateModelPart("Main")
        tester = GeometryTesterUtility()
        self.assertTrue( tester.TestTetrahedra3D10N(model_part) )

    def test_tetrahedra_3D8N(self):
        current_model = Model()
        model_part= current_model.CreateModelPart("Main")
        tester = GeometryTesterUtility()
        self.assertTrue( tester.TestHexahedra3D8N(model_part) )

    def test_tetrahedra_3D27N(self):
        current_model = Model()
        model_part= current_model.CreateModelPart("Main")
        tester = GeometryTesterUtility()
        self.assertTrue( tester.TestHexahedra3D27N(model_part) )
        
    def test_tetrahedra_3D20N(self):
        current_model = Model()
        model_part= current_model.CreateModelPart("Main")
        tester = GeometryTesterUtility()
        self.assertTrue( tester.TestHexahedra3D20N(model_part) )

    def test_quadrilateral_interface_2D4N(self):
        current_model = Model()
        model_part= current_model.CreateModelPart("Main")
        tester = GeometryTesterUtility()
        self.assertTrue( tester.TestQuadrilateralInterface2D4N(model_part) )

    def test_prism_interface_2D4N(self):
        current_model = Model()
        model_part= current_model.CreateModelPart("Main")
        tester = GeometryTesterUtility()
        self.assertTrue( tester.TestPrismInterface3D6N(model_part) )

    def test_hexahedra_interface_2D4N(self):
        current_model = Model()
        model_part= current_model.CreateModelPart("Main")
        tester = GeometryTesterUtility()
        self.assertTrue( tester.TestHexahedraInterface3D8N(model_part) )

    def create_geometry(model_part):
        current_model = Model()
        model_part= current_model.CreateModelPart("Main")

        model_part.CreateNewNode(1, 0.0, 0.0, 0)
        model_part.CreateNewNode(2, 2.5, 0.0, 0)
        model_part.CreateNewNode(3, 5.0, 0.0, 0)

        model_part.CreateNewNode(4, 0.0, 0.5, 0)
        model_part.CreateNewNode(5, 2.5, 0.5, 0)
        model_part.CreateNewNode(6, 5.0, 0.5, 0)

        model_part.CreateNewNode(7, 0.0, 1.0, 0)
        model_part.CreateNewNode(8, 2.5, 1.0, 0)
        model_part.CreateNewNode(9, 5.0, 1.0, 0)

        nodes = model_part.NodesArray

        knots_u = Vector(4)
        knots_u[0] = 0.0
        knots_u[1] = 0.0
        knots_u[2] = 5.0
        knots_u[3] = 5.0

        knots_v = Vector(4)
        knots_v[0] = 0.0
        knots_v[1] = 0.0
        knots_v[2] = 1.0
        knots_v[3] = 1.0

        surface = NurbsSurfaceGeometry3D(
            nodes,
            2,
            2,
            knots_u,
            knots_v)

if __name__ == '__main__':
    KratosUnittest.main()
