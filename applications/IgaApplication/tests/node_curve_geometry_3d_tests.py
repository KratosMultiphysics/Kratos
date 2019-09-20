from KratosMultiphysics import *
from KratosMultiphysics.IgaApplication import *

import KratosMultiphysics.KratosUnittest as KratosUnittest


def GetFilePath(fileName):
    return os.path.dirname(__file__) + "/" + fileName


class NodeCurveGeometry3DTests(KratosUnittest.TestCase):

    def setUp(self):
        # generate nodes = control points
        node_0 = Node(1, 0.0, 0.0, 0.0)
        node_1 = Node(2, 1.0, 1.0, 0.0)
        node_2 = Node(3, 2.0, 0.0, 0.0)

        # set weights of control points
        node_0.SetValue(NURBS_CONTROL_POINT_WEIGHT, 1.0)
        node_1.SetValue(NURBS_CONTROL_POINT_WEIGHT, 1.0)
        node_2.SetValue(NURBS_CONTROL_POINT_WEIGHT, 1.0)

        # create a nurbs curve
        curve = NodeCurveGeometry3D(Degree=2, NumberOfNodes=3)

        # set knots
        curve.SetKnot(Index=0, Value=0.0)
        curve.SetKnot(Index=1, Value=0.0)
        curve.SetKnot(Index=2, Value=1.0)
        curve.SetKnot(Index=3, Value=1.0)

        # set nodes
        curve.SetNode(Index=0, Value=node_0)
        curve.SetNode(Index=1, Value=node_1)
        curve.SetNode(Index=2, Value=node_2)

        # store data
        self.curve = curve

    def tearDown(self):
        pass

    def testPointAt(self):
        curve = self.curve

        point = curve.PointAt(T=0.5)

        self.assertAlmostEqual(point[0], 1.0)
        self.assertAlmostEqual(point[1], 0.5)
        self.assertAlmostEqual(point[2], 0.0)

    def testMoveNode(self):
        curve = self.curve

        node_1 = curve.Node(Index=1)
        node_1.Y = 2

        point = curve.PointAt(T=0.5)

        self.assertAlmostEqual(point[0], 1.0)
        self.assertAlmostEqual(point[1], 1.0)
        self.assertAlmostEqual(point[2], 0.0)

    def testSetWeight(self):
        curve = self.curve

        node_1 = curve.Node(Index=1)
        node_1.SetValue(NURBS_CONTROL_POINT_WEIGHT, 2.0)

        point = curve.PointAt(T=0.5)

        self.assertAlmostEqual(point[0], 1.0)
        self.assertAlmostEqual(point[1], 2.0 / 3.0)
        self.assertAlmostEqual(point[2], 0.0)

    def testValueAt(self):
        curve = self.curve

        node_1 = curve.Node(Index=1)
        node_1.SetValue(DISPLACEMENT_Y, 1.0)

        n = 10

        for i in range(n + 1):
            t = 1.0 / n * i

            displacement_at_t = curve.ValueAt(Variable=DISPLACEMENT, T=t)

            self.assertAlmostEqual(displacement_at_t[0], 0.0)
            self.assertAlmostEqual(displacement_at_t[1], 2.0 * t - 2.0 * t**2)
            self.assertAlmostEqual(displacement_at_t[2], 0.0)

        for i in range(n + 1):
            t = 1.0 / n * i

            displacement_y_at_t = curve.ValueAt(Variable=DISPLACEMENT_Y, T=t)

            self.assertAlmostEqual(displacement_y_at_t, 2.0 * t - 2.0 * t**2)

    def testNodeIndexOutOfRange(self):
        if Kernel.BuildType() == 'Release':
            self.skipTest(reason='Index checked only in Debug mode')

        with self.assertRaises(Exception):
            self.curve.Node(Index=-1)

        with self.assertRaises(Exception):
            self.curve.Node(Index=3)
