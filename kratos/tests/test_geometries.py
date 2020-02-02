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

    def test_nurbs_surface_3d(model_part):
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
            nodes, 2, 2, knots_u, knots_v)

        # check nodes
        self.assertAlmostEqual(surface[0].X, 0.0)
        self.assertAlmostEqual(surface[1].X, 2.5)

        # check knot spans
        self.assertAlmostEqual(surface.KnotsU()[3], 5.0)
        self.assertEqual(surface.NumberOfKnotsU(), 4)
        self.assertEqual(surface.NumberOfControlPointsU(), 3)

        # check rational
        self.assertFalse(surface.IsRational())

    def test_nurbs_curve_3d(model_part):
        current_model = Model()
        model_part= current_model.CreateModelPart("Main")

        model_part.CreateNewNode(1, 0, 0, 0)
        model_part.CreateNewNode(2, 3.3333333333333335, 1.6666666666666667, 0)
        model_part.CreateNewNode(3, 6.6666666666666661, 3.333333333333333, 0)
        model_part.CreateNewNode(4, 10, 5, 0)

        nodes = model_part.NodesArray

        knots = Vector(4)
        knots[0] = 0.0
        knots[1] = 0.0
        knots[2] = 0.0
        knots[3] = 11.180339887498949
        knots[4] = 11.180339887498949
        knots[5] = 11.180339887498949

        curve = NurbsCurveGeometry3D(
            nodes, 3, knots)

        # check nodes
        self.assertAlmostEqual(surface[0].X, 0.0)
        self.assertAlmostEqual(surface[1].X, 2.5)

        # check knot spans
        self.assertAlmostEqual(surface.Knots()[4], 11.180339887498949)
        self.assertEqual(surface.NumberOfKnots(), 6)
        self.assertEqual(surface.NumberOfControlPoints(), 4)

        # check rational
        self.assertFalse(curve.IsRational())

        # check general information
        self.assertEqual(curve.Dimension(), 1);
        self.assertEqual(curve.WorkingSpaceDimension(), 3);
        self.assertEqual(curve.LocalSpaceDimension(), 1);

if __name__ == '__main__':
    KratosUnittest.main()
