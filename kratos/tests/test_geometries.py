from __future__ import print_function, absolute_import, division

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics import *

class TestGeometry(KratosUnittest.TestCase):
    def setUp(self):
        pass

    def test_triangle_2D3N(self):
        current_model = Model()
        model_part= current_model.CreateModelPart("Main")
        tester = GeometryTesterUtility()
        self.assertTrue( tester.TestTriangle2D3N(model_part) )

    def test_triangle_2D6N(self):
        current_model = Model()
        model_part= current_model.CreateModelPart("Main")
        tester = GeometryTesterUtility()
        self.assertTrue( tester.TestTriangle2D6N(model_part) )

    def test_quadrilateral_2D4N(self):
        current_model = Model()
        model_part= current_model.CreateModelPart("Main")
        tester = GeometryTesterUtility()
        self.assertTrue( tester.TestQuadrilateral2D4N(model_part) )

    def test_quadrilateral_2D9N(self):
        current_model = Model()
        model_part= current_model.CreateModelPart("Main")
        tester = GeometryTesterUtility()
        self.assertTrue( tester.TestQuadrilateral2D9N(model_part) )

    def test_tetrahedra_3D4N(self):
        current_model = Model()
        model_part= current_model.CreateModelPart("Main")
        tester = GeometryTesterUtility()
        self.assertTrue( tester.TestTetrahedra3D4N(model_part) )

    def test_tetrahedra_3D10N(self):
        current_model = Model()
        model_part= current_model.CreateModelPart("Main")
        tester = GeometryTesterUtility()
        self.assertTrue( tester.TestTetrahedra3D10N(model_part) )

    def test_prism_3D6N(self):
        current_model = Model()
        model_part= current_model.CreateModelPart("Main")
        tester = GeometryTesterUtility()
        self.assertTrue( tester.TestPrism3D6N(model_part) )

    def test_hexahedra_3D8N(self):
        current_model = Model()
        model_part= current_model.CreateModelPart("Main")
        tester = GeometryTesterUtility()
        self.assertTrue( tester.TestHexahedra3D8N(model_part) )

    def test_hexahedra_3D20N(self):
        current_model = Model()
        model_part= current_model.CreateModelPart("Main")
        tester = GeometryTesterUtility()
        self.assertTrue( tester.TestHexahedra3D20N(model_part) )

    def test_hexahedra_3D27N(self):
        current_model = Model()
        model_part= current_model.CreateModelPart("Main")
        tester = GeometryTesterUtility()
        self.assertTrue( tester.TestHexahedra3D27N(model_part) )

if __name__ == '__main__':
    KratosUnittest.main()
