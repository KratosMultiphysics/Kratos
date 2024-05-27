import KratosMultiphysics
import KratosMultiphysics.MeshMovingApplication as MeshMoving
import KratosMultiphysics.KratosUnittest as UnitTest

import math

class AffineTransformTest(UnitTest.TestCase):

    @staticmethod
    def GeneratePoints():
        Array3 = KratosMultiphysics.Array3
        return [Array3([2.0, 0.0, 0.0]),
                Array3([0.0, 2.0, 0.0]),
                Array3([0.0, 0.0, 2.0])]

    def CheckPoints(self, points: list, **kwargs):
        self.assertVectorAlmostEqual(points[0], [0.0, 2.0, 6.0], **kwargs)
        self.assertVectorAlmostEqual(points[1], [0.0, 4.0, 4.0], **kwargs)
        self.assertVectorAlmostEqual(points[2], [-2.0, 2.0, 4.0], **kwargs)

    def test_AxisAngle(self):
        axis = [0.0, 1.0, 0.0]
        angle = -math.pi / 2.0
        reference_point = [-1.0, 0.0, 0.0]
        translation_vector = [1.0, 2.0, 3.0]

        transform = MeshMoving.AffineTransform(
            axis,
            angle,
            reference_point,
            translation_vector)

        self.CheckPoints([transform.Apply(point) for point in self.GeneratePoints()])

    def test_EulerAngles(self):
        euler_angles = [-math.pi/2.0, -math.pi/2.0, math.pi/2.0]
        reference_point = [-1.0, 0.0, 0.0]
        translation_vector = [1.0, 2.0, 3.0]

        transform = MeshMoving.AffineTransform(
            euler_angles,
            reference_point,
            translation_vector)

        self.CheckPoints([transform.Apply(point) for point in self.GeneratePoints()])

    def test_Quaternion(self):
        quaternion = KratosMultiphysics.Quaternion()
        quaternion.X = 0.0
        quaternion.Y = -math.sqrt(2.0) / 2.0
        quaternion.Z = 0.0
        quaternion.W = math.sqrt(2.0) / 2.0
        reference_point = [-1.0, 0.0, 0.0]
        translation_vector = [1.0, 2.0, 3.0]

        transform = MeshMoving.AffineTransform(
            quaternion,
            reference_point,
            translation_vector)

        self.CheckPoints([transform.Apply(point) for point in self.GeneratePoints()])

if __name__ == "__main__":
    UnitTest.main()