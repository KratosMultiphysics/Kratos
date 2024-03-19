import KratosMultiphysics
from KratosMultiphysics import Parameters
import KratosMultiphysics.MeshMovingApplication as MeshMoving
import KratosMultiphysics.KratosUnittest as UnitTest

import math

class ParametricAffineTransformTest(UnitTest.TestCase):

    @staticmethod
    def GeneratePoints():
        Array3 = KratosMultiphysics.Array3
        return [Array3([2.0, 0.0, 0.0]),
                Array3([0.0, 2.0, 0.0]),
                Array3([0.0, 0.0, 2.0])]

    def CheckDefaultTransformedPoints(self, points: list):
        self.assertVectorAlmostEqual(points[0], [0.0, 2.0, 6.0], places=5)
        self.assertVectorAlmostEqual(points[1], [0.0, 4.0, 4.0], places=5)
        self.assertVectorAlmostEqual(points[2], [-2.0, 2.0, 4.0], places=5)

    def test_ConstantTransform(self):
        points = self.GeneratePoints()
        transform = MeshMoving.ParametricAffineTransform(
            self.GetConstantAxis(),
            self.GetConstantAngle(),
            self.GetConstantReferencePoint(),
            self.GetConstantTranslationVector())

        self.CheckDefaultTransformedPoints([transform.Apply(point, 0.0, 0.0, 0.0, 0.0) for point in points])

    def test_ConstantTransformEuler(self):
        points = self.GeneratePoints()
        transform = MeshMoving.ParametricAffineTransform(
            self.GetConstantEulerAngles(),
            self.GetConstantReferencePoint(),
            self.GetConstantTranslationVector())

        self.CheckDefaultTransformedPoints([transform.Apply(point, 0.0, 0.0, 0.0, 0.0) for point in points])

    def test_ParametricAngle(self):
        points = self.GeneratePoints()
        transform = MeshMoving.ParametricAffineTransform(
            self.GetConstantAxis(),
            Parameters(f""" "-(t+1.0) * {math.pi} / 2.0" """),
            self.GetConstantReferencePoint(),
            self.GetConstantTranslationVector())

        self.CheckDefaultTransformedPoints([transform.Apply(point, 0.0, 0.0, 0.0, 0.0) for point in points])

        self.assertVectorAlmostEqual(transform.Apply(points[0], 1.0, 0.0, 0.0, 0.0), [-3.0, 2.0, 3.0], places=5)
        self.assertVectorAlmostEqual(transform.Apply(points[1], 1.0, 0.0, 0.0, 0.0), [-1.0, 4.0, 3.0], places=5)
        self.assertVectorAlmostEqual(transform.Apply(points[2], 1.0, 0.0, 0.0, 0.0), [-1.0, 2.0, 1.0], places=5)

    def test_ParametricAxis(self):
        points = self.GeneratePoints()
        transform = MeshMoving.ParametricAffineTransform(
            Parameters(""" ["t", "1.0-t", 0.0] """),
            self.GetConstantAngle(),
            self.GetConstantReferencePoint(),
            self.GetConstantTranslationVector())

        self.CheckDefaultTransformedPoints([transform.Apply(point, 0.0, 0.0, 0.0, 0.0) for point in points])

        self.assertVectorAlmostEqual(transform.Apply(points[0], 1.0, 0.0, 0.0, 0.0), [3.0, 2.0, 3.0], places=5)
        self.assertVectorAlmostEqual(transform.Apply(points[1], 1.0, 0.0, 0.0, 0.0), [1.0, 2.0, 1.0], places=5)
        self.assertVectorAlmostEqual(transform.Apply(points[2], 1.0, 0.0, 0.0, 0.0), [1.0, 4.0, 3.0], places=5)

    def test_ParametricEulerAngles(self):
        points = self.GeneratePoints()
        transform = MeshMoving.ParametricAffineTransform(
            Parameters(f""" [{-math.pi/2.0}, {-math.pi/2.0}, "(1.0-t)*{math.pi/2.0}"] """),
            self.GetConstantReferencePoint(),
            self.GetConstantTranslationVector())

        self.CheckDefaultTransformedPoints([transform.Apply(point, 0.0, 0.0, 0.0, 0.0) for point in points])

        self.assertVectorAlmostEqual(transform.Apply(points[0], 1.0, 0.0, 0.0, 0.0), [0.0, -1.0, 3.0], places=5)
        self.assertVectorAlmostEqual(transform.Apply(points[1], 1.0, 0.0, 0.0, 0.0), [0.0, 1.0, 5.0], places=5)
        self.assertVectorAlmostEqual(transform.Apply(points[2], 1.0, 0.0, 0.0, 0.0), [-2.0, 1.0, 3.0], places=5)

    def test_ParametricReferencePoint(self):
        points = self.GeneratePoints()
        transform = MeshMoving.ParametricAffineTransform(
            self.GetConstantAxis(),
            self.GetConstantAngle(),
            Parameters(""" ["t-1.0", 0.0, 0.0] """),
            self.GetConstantTranslationVector())

        self.CheckDefaultTransformedPoints([transform.Apply(point, 0.0, 0.0, 0.0, 0.0) for point in points])

        self.assertVectorAlmostEqual(transform.Apply(points[0], 1.0, 0.0, 0.0, 0.0), [1.0, 2.0, 5.0], places=5)
        self.assertVectorAlmostEqual(transform.Apply(points[1], 1.0, 0.0, 0.0, 0.0), [1.0, 4.0, 3.0], places=5)
        self.assertVectorAlmostEqual(transform.Apply(points[2], 1.0, 0.0, 0.0, 0.0), [-1.0, 2.0, 3.0], places=5)

    def test_ParametricTranslationVector(self):
        points = self.GeneratePoints()
        transform = MeshMoving.ParametricAffineTransform(
            self.GetConstantAxis(),
            self.GetConstantAngle(),
            self.GetConstantReferencePoint(),
            Parameters(""" ["1.0-t", "2.0*(1.0-t)", "3.0*(1.0-t)"] """))

        self.CheckDefaultTransformedPoints([transform.Apply(point, 0.0, 0.0, 0.0, 0.0) for point in points])

        self.assertVectorAlmostEqual(transform.Apply(points[0], 1.0, 0.0, 0.0, 0.0), [-1.0, 0.0, 3.0], places=5)
        self.assertVectorAlmostEqual(transform.Apply(points[1], 1.0, 0.0, 0.0, 0.0), [-1.0, 2.0, 1.0], places=5)
        self.assertVectorAlmostEqual(transform.Apply(points[2], 1.0, 0.0, 0.0, 0.0), [-3.0, 0.0, 1.0], places=5)

    @staticmethod
    def GetConstantAxis():
        return Parameters("[0.0, 1.0, 0.0]")

    @staticmethod
    def GetConstantAngle():
        return Parameters(f"{-math.pi/2.0}")

    @staticmethod
    def GetConstantReferencePoint():
        return Parameters("[-1.0, 0.0, 0.0]")

    @staticmethod
    def GetConstantTranslationVector():
        return Parameters("[1.0, 2.0, 3.0]")

    @staticmethod
    def GetConstantEulerAngles():
        return Parameters(f"[{-math.pi/2},{-math.pi/2},{math.pi/2}]")

if __name__ == "__main__":
    UnitTest.main()