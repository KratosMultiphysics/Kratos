from KratosMultiphysics import *
from KratosMultiphysics.IgaApplication import *

import KratosMultiphysics.KratosUnittest as KratosUnittest


def GetFilePath(fileName):
    return os.path.dirname(__file__) + "/" + fileName


class NodeSurfaceGeometry3DTests(KratosUnittest.TestCase):

    def setUp(self):
        # generate nodes = control points
        node_00 = Node(1, 0.0, 0.0, 0.0)
        node_01 = Node(2, 0.0, 1.0, 0.0)
        node_02 = Node(3, 0.0, 2.0, 0.0)
        node_10 = Node(1, 1.0, 0.0, 0.0)
        node_11 = Node(2, 1.0, 1.0, 1.0)
        node_12 = Node(3, 1.0, 2.0, 0.0)
        node_20 = Node(1, 2.0, 0.0, 0.0)
        node_21 = Node(2, 2.0, 1.0, 0.0)
        node_22 = Node(3, 2.0, 2.0, 0.0)

        # set weights of control points
        node_00.SetValue(NURBS_CONTROL_POINT_WEIGHT, 1.0)
        node_01.SetValue(NURBS_CONTROL_POINT_WEIGHT, 1.0)
        node_02.SetValue(NURBS_CONTROL_POINT_WEIGHT, 1.0)
        node_10.SetValue(NURBS_CONTROL_POINT_WEIGHT, 1.0)
        node_11.SetValue(NURBS_CONTROL_POINT_WEIGHT, 1.0)
        node_12.SetValue(NURBS_CONTROL_POINT_WEIGHT, 1.0)
        node_20.SetValue(NURBS_CONTROL_POINT_WEIGHT, 1.0)
        node_21.SetValue(NURBS_CONTROL_POINT_WEIGHT, 1.0)
        node_22.SetValue(NURBS_CONTROL_POINT_WEIGHT, 1.0)

        # create a nurbs surface
        surface = NodeSurfaceGeometry3D(DegreeU=2, DegreeV=2, NumberOfNodesU=3,
            NumberOfNodesV=3)

        # set knots in u-direction
        surface.SetKnotU(Index=0, Value=0.0)
        surface.SetKnotU(Index=1, Value=0.0)
        surface.SetKnotU(Index=2, Value=1.0)
        surface.SetKnotU(Index=3, Value=1.0)

        # set knots in v-direction
        surface.SetKnotV(Index=0, Value=0.0)
        surface.SetKnotV(Index=1, Value=0.0)
        surface.SetKnotV(Index=2, Value=1.0)
        surface.SetKnotV(Index=3, Value=1.0)

        # set nodes
        surface.SetNode(IndexU=0, IndexV=0, Value=node_00)
        surface.SetNode(IndexU=0, IndexV=1, Value=node_01)
        surface.SetNode(IndexU=0, IndexV=2, Value=node_02)
        surface.SetNode(IndexU=1, IndexV=0, Value=node_10)
        surface.SetNode(IndexU=1, IndexV=1, Value=node_11)
        surface.SetNode(IndexU=1, IndexV=2, Value=node_12)
        surface.SetNode(IndexU=2, IndexV=0, Value=node_20)
        surface.SetNode(IndexU=2, IndexV=1, Value=node_21)
        surface.SetNode(IndexU=2, IndexV=2, Value=node_22)

        # store data
        self.surface = surface

    def tearDown(self):
        pass

    def testPointAt(self):
        surface = self.surface

        point = surface.PointAt(U=0.5, V=0.5)

        self.assertAlmostEqual(point[0], 1.00)
        self.assertAlmostEqual(point[1], 1.00)
        self.assertAlmostEqual(point[2], 0.25)

    def testMoveNode(self):
        surface = self.surface

        node_11 = surface.Node(IndexU=1, IndexV=1)
        node_11.Z = 2

        point = surface.PointAt(U=0.5, V=0.5)

        self.assertAlmostEqual(point[0], 1.00)
        self.assertAlmostEqual(point[1], 1.00)
        self.assertAlmostEqual(point[2], 0.50)

    def testSetWeight(self):
        surface = self.surface

        node_11 = surface.Node(IndexU=1, IndexV=1)
        node_11.SetValue(NURBS_CONTROL_POINT_WEIGHT, 2.0)

        point = surface.PointAt(U=0.5, V=0.5)

        self.assertAlmostEqual(point[0], 1.0)
        self.assertAlmostEqual(point[1], 1.0)
        self.assertAlmostEqual(point[2], 0.4)

    def testValueAt(self):
        surface = self.surface

        node_11 = surface.Node(IndexU=1, IndexV=1)
        node_11.SetValue(DISPLACEMENT_Z, 1.0)

        n = 10

        for i in range(n + 1):
            t = 1.0 / n * i

            displacement_at_t = surface.ValueAt(Variable=DISPLACEMENT, U=t,
                V=0.5)

            self.assertAlmostEqual(displacement_at_t[0], 0.0)
            self.assertAlmostEqual(displacement_at_t[1], 0.0)
            self.assertAlmostEqual(displacement_at_t[2], t - t**2)

        for i in range(n + 1):
            t = 1.0 / n * i

            displacement_z_at_t = surface.ValueAt(Variable=DISPLACEMENT_Z, U=t,
                V=0.5)

            self.assertAlmostEqual(displacement_z_at_t, t - t**2)

    def testNodeIndexOutOfRange(self):
        if Kernel.BuildType() == 'Release':
            self.skipTest(reason='Index checked only in Debug mode')

        with self.assertRaises(Exception):
            self.surface.Node(IndexU=-1, IndexV=0)

        with self.assertRaises(Exception):
            self.surface.Node(IndexU=3, IndexV=0)

        with self.assertRaises(Exception):
            self.surface.Node(IndexU=0, IndexV=-1)

        with self.assertRaises(Exception):
            self.surface.Node(IndexU=0, IndexV=3)
