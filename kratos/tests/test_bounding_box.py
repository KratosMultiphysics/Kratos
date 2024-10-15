import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestBoundingBox(KratosUnittest.TestCase):
    def CreateModelPart(self):
        current_model = KM.Model()
        model_part = current_model.CreateModelPart("TestModelPart")
        model_part.AddNodalSolutionStepVariable(KM.VORTICITY)
        model_part.CreateNewNode(1,  0.0, 10.0,  0.0)
        model_part.CreateNewNode(2, -1.0,  0.0,  0.0)
        model_part.CreateNewNode(3,  1.0, -5.0, 80.0)

        model_part.AddProperties(KM.Properties(1))
        model_part.CreateNewElement("Element2D3N", 1, [1,2,3], model_part.GetProperties()[1])
        model_part.CreateNewCondition("SurfaceCondition3D3N", 1, [1,2,3], model_part.GetProperties()[1])
        return model_part

    def test_BoundingBoxWithNodesArray(self):
        model_part = self.CreateModelPart()
        bb = KM.BoundingBox(model_part.Nodes)
        self._CheckBoundingBoxPoints(bb)

    def test_BoundingBoxWithPoints(self):
        min_point = KM.Point(-1.0, -5.0,  0.0)
        max_point = KM.Point( 1.0, 10.0, 80.0)
        bb = KM.BoundingBox(min_point, max_point)
        self._CheckBoundingBoxPoints(bb)

    def test_BoundingBoxSet(self):
        model_part = self.CreateModelPart()
        min_point = KM.Point(-100.0, -100.0, -100.0)
        max_point = KM.Point( 100.0,  100.0,  100.0)
        bb = KM.BoundingBox(min_point, max_point)
        bb.Set(model_part.Nodes)
        self._CheckBoundingBoxPoints(bb)

    def test_BoundingBoxExtend(self):
        model_part = self.CreateModelPart()
        min_point = KM.Point(0.0, 0.0, 0.0)
        max_point = KM.Point(0.0, 0.0, 0.0)
        bb = KM.BoundingBox(min_point, max_point)
        bb.Extend(model_part.Nodes)
        self._CheckBoundingBoxPoints(bb)

    def test_BoundingBoxExtendWithDouble(self):
        min_point = KM.Point(0.0, -4.0,  1.0)
        max_point = KM.Point(0.0,  9.0, 79.0)
        bb = KM.BoundingBox(min_point, max_point)
        bb.Extend(1.0)
        self._CheckBoundingBoxPoints(bb)

    def test_BoundingBoxIsInside(self):
        min_point = KM.Point(-1.0, -5.0,  0.0)
        max_point = KM.Point( 1.0, 10.0, 80.0)
        bb = KM.BoundingBox(min_point, max_point)
        point_inside = KM.Point(0.0, 0.0, 0.0)
        self.assertTrue(bb.IsInside(point_inside))
        point_outside = KM.Point(0.0, 0.0, 100.0)
        self.assertFalse(bb.IsInside(point_outside))

    def _CheckBoundingBoxPoints(self, bounding_box):
        min_point, max_point = bounding_box.GetPoints()

        self.assertAlmostEqual(min_point.X, -1.0)
        self.assertAlmostEqual(min_point.Y, -5.0)
        self.assertAlmostEqual(min_point.Z,  0.0)

        self.assertAlmostEqual(max_point.X,  1.0)
        self.assertAlmostEqual(max_point.Y, 10.0)
        self.assertAlmostEqual(max_point.Z, 80.0)

if __name__ == "__main__":
    KratosUnittest.main()
