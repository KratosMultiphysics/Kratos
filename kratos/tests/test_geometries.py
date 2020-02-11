from __future__ import print_function, absolute_import, division

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics as KM

class TestGeometry(KratosUnittest.TestCase):
    def setUp(self):
        pass

    def test_tetrahedra_3D4N(self):
        current_model = KM.Model()
        model_part= current_model.CreateModelPart("Main")
        tester = KM.GeometryTesterUtility()
        self.assertTrue( tester.TestTetrahedra3D4N(model_part) )

    def test_tetrahedra_2D3N(self):
        current_model = KM.Model()
        model_part= current_model.CreateModelPart("Main")
        tester = KM.GeometryTesterUtility()
        self.assertTrue( tester.TestTriangle2D3N(model_part) )

    def test_tetrahedra_2D6N(self):
        current_model = KM.Model()
        model_part= current_model.CreateModelPart("Main")
        tester = KM.GeometryTesterUtility()
        self.assertTrue( tester.TestTriangle2D6N(model_part) )

    def test_tetrahedra_3D10N(self):
        current_model = KM.Model()
        model_part= current_model.CreateModelPart("Main")
        tester = KM.GeometryTesterUtility()
        self.assertTrue( tester.TestTetrahedra3D10N(model_part) )

    def test_tetrahedra_3D8N(self):
        current_model = KM.Model()
        model_part= current_model.CreateModelPart("Main")
        tester = KM.GeometryTesterUtility()
        self.assertTrue( tester.TestHexahedra3D8N(model_part) )

    def test_tetrahedra_3D27N(self):
        current_model = KM.Model()
        model_part= current_model.CreateModelPart("Main")
        tester = KM.GeometryTesterUtility()
        self.assertTrue( tester.TestHexahedra3D27N(model_part) )
        
    def test_tetrahedra_3D20N(self):
        current_model = KM.Model()
        model_part= current_model.CreateModelPart("Main")
        tester = KM.GeometryTesterUtility()
        self.assertTrue( tester.TestHexahedra3D20N(model_part) )

    def test_quadrilateral_interface_2D4N(self):
        current_model = KM.Model()
        model_part= current_model.CreateModelPart("Main")
        tester = KM.GeometryTesterUtility()
        self.assertTrue( tester.TestQuadrilateralInterface2D4N(model_part) )

    def test_prism_interface_2D4N(self):
        current_model = KM.Model()
        model_part= current_model.CreateModelPart("Main")
        tester = KM.GeometryTesterUtility()
        self.assertTrue( tester.TestPrismInterface3D6N(model_part) )

    def test_hexahedra_interface_2D4N(self):
        current_model = KM.Model()
        model_part= current_model.CreateModelPart("Main")
        tester = KM.GeometryTesterUtility()
        self.assertTrue( tester.TestHexahedraInterface3D8N(model_part) )

    def test_nurbs_surface_3d(self):
        current_model = KM.Model()
        model_part= current_model.CreateModelPart("Main")

        node1 = model_part.CreateNewNode(1, 0.0, 0.0, 0)
        node2 = model_part.CreateNewNode(2, 2.5, 0.0, 0)
        node3 = model_part.CreateNewNode(3, 5.0, 0.0, 0)

        node4 = model_part.CreateNewNode(4, 0.0, 0.5, 0)
        node5 = model_part.CreateNewNode(5, 2.5, 0.5, 0)
        node6 = model_part.CreateNewNode(6, 5.0, 0.5, 0)

        node7 = model_part.CreateNewNode(7, 0.0, 1.0, 0)
        node8 = model_part.CreateNewNode(8, 2.5, 1.0, 0)
        node9 = model_part.CreateNewNode(9, 5.0, 1.0, 0)

        nodes = KM.NodesVector()
        nodes.append(node1)
        nodes.append(node2)
        nodes.append(node3)
        nodes.append(node4)
        nodes.append(node5)
        nodes.append(node6)
        nodes.append(node7)
        nodes.append(node8)
        nodes.append(node9)

        knots_u = KM.Vector(4)
        knots_u[0] = 0.0
        knots_u[1] = 0.0
        knots_u[2] = 5.0
        knots_u[3] = 5.0

        knots_v = KM.Vector(4)
        knots_v[0] = 0.0
        knots_v[1] = 0.0
        knots_v[2] = 1.0
        knots_v[3] = 1.0

        surface = KM.NurbsSurfaceGeometry3D(
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

    def test_nurbs_curve_3d(self):
        current_model = KM.Model()
        model_part= current_model.CreateModelPart("Main")

        node1 = model_part.CreateNewNode(1, 0, 0, 0)
        node2 = model_part.CreateNewNode(2, 3.3333333333333335, 1.6666666666666667, 0)
        node3 = model_part.CreateNewNode(3, 6.6666666666666661, 3.333333333333333, 0)
        node4 = model_part.CreateNewNode(4, 10, 5, 0)

        nodes = KM.NodesVector()
        nodes.append(node1)
        nodes.append(node2)
        nodes.append(node3)
        nodes.append(node4)

        knots = KM.Vector(6)
        knots[0] = 0.0
        knots[1] = 0.0
        knots[2] = 0.0
        knots[3] = 11.180339887498949
        knots[4] = 11.180339887498949
        knots[5] = 11.180339887498949

        curve = KM.NurbsCurveGeometry3D(
            nodes, 3, knots)

        # check nodes
        self.assertAlmostEqual(curve[0].X, 0.0)
        self.assertAlmostEqual(curve[1].X, 3.3333333333333335)

        # checking right pointers/ sources of nodes
        curve[1].X = 2
        self.assertAlmostEqual(curve[1].X, node2.X)
        self.assertAlmostEqual(node2.X, 2)

        # check knot spans
        self.assertAlmostEqual(curve.Knots()[4], 11.180339887498949)
        self.assertEqual(curve.NumberOfKnots(), 6)
        self.assertEqual(curve.NumberOfControlPoints(), 4)

        # check rational
        self.assertFalse(curve.IsRational())

        # check general information
        self.assertEqual(curve.Dimension(), 1);
        self.assertEqual(curve.WorkingSpaceDimension(), 3);
        self.assertEqual(curve.LocalSpaceDimension(), 1);

    def test_nodes_vector_iterators(self):
        current_model = KM.Model()
        model_part= current_model.CreateModelPart("Main")

        node1 = model_part.CreateNewNode(1, 0, 0, 0)
        node2 = model_part.CreateNewNode(2, 3.3333333333333335, 1.6666666666666667, 0)
        node3 = model_part.CreateNewNode(3, 6.6666666666666661, 3.333333333333333, 0)
        node4 = model_part.CreateNewNode(4, 10, 5, 0)

        nodes = KM.NodesVector()
        nodes.append(node1)
        nodes.append(node2)
        nodes.append(node3)
        nodes.append(node4)

        for node in nodes:
            node.X = node.X + 1

        # check new coordinates
        self.assertAlmostEqual(nodes[0].X, node1.X)
        self.assertAlmostEqual(nodes[0].X, 1.0)

        # check length of nodes vector
        self.assertEqual(len(nodes), 4)

    def test_nodes_vector_geometry_iterators(self):
        current_model = KM.Model()
        model_part= current_model.CreateModelPart("Main")

        node1 = model_part.CreateNewNode(1, 0, 0, 0)
        node2 = model_part.CreateNewNode(2, 3.3333333333333335, 1.6666666666666667, 0)
        node3 = model_part.CreateNewNode(3, 6.6666666666666661, 3.333333333333333, 0)
        node4 = model_part.CreateNewNode(4, 10, 5, 0)

        nodes = KM.NodesVector()
        nodes.append(node1)
        nodes.append(node2)
        nodes.append(node3)
        nodes.append(node4)

        knots = KM.Vector(6)
        knots[0] = 0.0
        knots[1] = 0.0
        knots[2] = 0.0
        knots[3] = 11.180339887498949
        knots[4] = 11.180339887498949
        knots[5] = 11.180339887498949

        curve = KM.NurbsCurveGeometry3D(
            nodes, 3, knots)

        for node in curve:
            node.X = node.X + 1

        # check new coordinates
        self.assertAlmostEqual(nodes[0].X, node1.X)
        self.assertAlmostEqual(nodes[0].X, 1.0)

        # check length of geometry node vector
        self.assertEqual(len(curve), 4)

if __name__ == '__main__':
    KratosUnittest.main()
