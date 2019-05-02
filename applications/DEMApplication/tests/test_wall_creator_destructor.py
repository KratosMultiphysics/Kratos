import KratosMultiphysics as Kratos
from KratosMultiphysics.DEMApplication import *
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestWallCreatorDestructor(KratosUnittest.TestCase):

    def setUp(self):
        self.current_model = Kratos.Model()
        self.walls_model_part = self.current_model.CreateModelPart("Walls")
        self.walls_model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)

        properties = Kratos.Properties(0)
        properties[Kratos.YOUNG_MODULUS] = 3.331

        self.walls_model_part.AddProperties(properties)


    def test_CreateWallTriangle(self):

        self.CreateNodes()

        condition_name = "RigidFace3D3N"

        self.walls_model_part.CreateNewCondition(condition_name, 1, [1, 2, 3], self.walls_model_part.GetProperties()[0])

        counter = 0
        for condition in self.walls_model_part.Conditions:
            counter += 1
            self.assertEqual(condition.Properties[Kratos.YOUNG_MODULUS], 3.331)

        self.assertEqual(counter, 1)

    def CreateNodes(self):

        self.walls_model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        self.walls_model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        self.walls_model_part.CreateNewNode(3, 0.0, 1.0, 0.0)

        counter = 0
        for node in self.walls_model_part.Nodes:
            counter += 1

        self.assertEqual(counter, 3)


    def test_CreateWallQuadrilateral(self):
        pass


if __name__ == '__main__':
    KratosUnittest.main()
