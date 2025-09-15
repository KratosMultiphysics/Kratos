# Import Kratos core and apps
import KratosMultiphysics as KM

# Additional imports
from KratosMultiphysics.KratosUnittest import TestCase
from KratosMultiphysics.ShapeOptimizationApplication.custom_ios.wrl_io import WrlIO


class WrlIOTest(TestCase):

    def test_simple_model_part(self):
        with KM.KratosUnittest.WorkFolderScope(".", __file__):
            model = KM.Model()
            model_part = model.CreateModelPart("wrl_model_part")
            model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
            model_part_io = WrlIO("joined")
            model_part_io.ReadModelPart(model_part)

        self.assertEqual(model_part.NumberOfNodes(), 20)
        self.assertEqual(model_part.NumberOfConditions(), 18)

        coords = [
            [-10, 0, 0],
            [-10, 0, 5],
            [-10, 0, 10],
            [-5, 0, 5],
            [-5, 0, 10],
            [0, 0, 0],
            [0, 0, 5],
            [0, 0, 10],
            [-6.666666507720947, 0, 0],
            [-3.333333253860474, 0, 0],
            [-10, 0, 0],
            [-10, 5, 0],
            [-10, 10, 0],
            [-5, 5, 0],
            [-5, 10, 0],
            [0, 0, 0],
            [0, 5, 0],
            [0, 10, 0],
            [-6.666666507720947, 0, 0],
            [-3.333333253860474, 0, 0]
        ]

        i = 0
        for node in model_part.Nodes:
            self.assertAlmostEqual(node.X, coords[i][0])
            self.assertAlmostEqual(node.Y, coords[i][1])
            self.assertAlmostEqual(node.Z, coords[i][2])
            i += 1

        topology = [
               [ 4,    2,    1],
			   [ 4,    1,    3],
			   [ 7,    4,    3],
			   [ 7,    3,    6],
			   [ 1,    0,    8],
			   [ 3,    1,    8],
			   [ 9,    5,    6],
			   [ 9,    6,    3],
			   [ 8,    9,    3],
			   [14,   13,   11],
			   [14,   11,   12],
			   [17,   16,   13],
			   [17,   13,   14],
			   [11,   18,   10],
			   [13,   18,   11],
			   [19,   13,   16],
			   [19,   16,   15],
			   [18,   13,   19]
            ]

        i = 0
        for condition in model_part.Conditions:
            j = 0
            for node in condition.GetNodes():
                self.assertEqual(node.Id, topology[i][j])
                j += 1
            i+=1

    def test_split_model_part(self):
        with KM.KratosUnittest.WorkFolderScope(".", __file__):
            model = KM.Model()
            model_part = model.CreateModelPart("wrl_model_part")
            model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
            model_part_io = WrlIO("split")
            model_part_io.ReadModelPart(model_part)

        self.assertEqual(model_part.NumberOfNodes(), 18)
        self.assertEqual(model_part.NumberOfConditions(), 16)

        sub_model_part_1 = model.GetModelPart("wrl_model_part.triangle mesh")

        self.assertEqual(sub_model_part_1.NumberOfNodes(), 9)
        self.assertEqual(sub_model_part_1.NumberOfConditions(), 8)

        coords = [
            [-10, 0, 0],
            [-10, 0, 5],
            [-10, 0, 10],
            [-5, 0, 0],
            [-5, 0, 5],
            [-5, 0, 10],
            [0, 0, 0],
            [0, 0, 5],
            [0, 0, 10]
        ]

        i = 0
        for node in sub_model_part_1.Nodes:
            self.assertAlmostEqual(node.X, coords[i][0])
            self.assertAlmostEqual(node.Y, coords[i][1])
            self.assertAlmostEqual(node.Z, coords[i][2])
            i += 1

        topology = [
            [4,    1,    0],
            [4,    0,    3],
            [5,    2,    1],
            [5,    1,    4],
            [7,    4,    3],
            [7,    3,    6],
            [8,    5,    4],
            [8,    4,    7],
        ]

        i = 0
        for condition in sub_model_part_1.Conditions:
            j = 0
            for node in condition.GetNodes():
                self.assertEqual(node.Id, topology[i][j])
                j += 1
            i+=1

        sub_model_part_2 = model.GetModelPart("wrl_model_part.triangle mesh_2")

        self.assertEqual(sub_model_part_2.NumberOfNodes(), 9)
        self.assertEqual(sub_model_part_2.NumberOfConditions(), 8)

if __name__ == '__main__':
    KM.KratosUnittest.main()
