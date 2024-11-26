import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
import KratosMultiphysics.KratosUnittest as UnitTest

class TestDistanceMatrix(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")

        for i in range(101):
            cls.model_part.CreateNewNode(i + 1, i + 1, i + 2, i + 3)

        cls.distance_matrix = KratosSI.DistanceMatrix()

        exp = Kratos.Expression.NodalExpression(cls.model_part)
        Kratos.Expression.NodalPositionExpressionIO.Read(exp, Kratos.Configuration.Initial)
        cls.distance_matrix.Update(exp)

    def test_GetDistance1(self) -> None:
        for i, node_i in enumerate(self.model_part.Nodes):
            for j, node_j in enumerate(self.model_part.Nodes):
                distance = ((node_i.X - node_j.X) ** 2 + (node_i.Y - node_j.Y) ** 2 + (node_i.Z - node_j.Z) ** 2) ** (0.5)
                self.assertAlmostEqual(self.distance_matrix.GetDistance(i, j), distance)

if __name__ == '__main__':
    UnitTest.main()