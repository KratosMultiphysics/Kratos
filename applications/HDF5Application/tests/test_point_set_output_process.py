import KratosMultiphysics
from KratosMultiphysics import KratosUnittest as UnitTest
from KratosMultiphysics.HDF5Application.point_set_output_process import Factory as PointSetOuputProcessFactory


class TestPointSetOutputProcess(UnitTest.TestCase):

    @staticmethod
    def MakeModelPart(model: KratosMultiphysics.Model):
        model_part = model.CreateModelPart("main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)

        node = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, [0.0, 0.0, 0.0])
        node.SetSolutionStepValue(KratosMultiphysics.REACTION, [0.0, 0.0, 0.0])

        node = model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, [1.0, 0.0, 0.0])
        node.SetSolutionStepValue(KratosMultiphysics.REACTION, [0.0, 0.0, 1.0])

        node = model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
        node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, [2.0, 0.0, 0.0])
        node.SetSolutionStepValue(KratosMultiphysics.REACTION, [0.0, 0.0, 2.0])

        node = model_part.CreateNewNode(4, 0.0, 1.0, 0.0)
        node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, [3.0, 0.0, 0.0])
        node.SetSolutionStepValue(KratosMultiphysics.REACTION, [0.0, 0.0, 3.0])

        node = model_part.CreateNewNode(5, 2.0, 0.0, 0.0)
        node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, [4.0, 0.0, 0.0])
        node.SetSolutionStepValue(KratosMultiphysics.REACTION, [0.0, 0.0, 4.0])

        node = model_part.CreateNewNode(6, 2.0, 1.0, 0.0)
        node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, [5.0, 0.0, 0.0])
        node.SetSolutionStepValue(KratosMultiphysics.REACTION, [0.0, 0.0, 5.0])

        model_part.CreateNewElement(
            "Element2D4N",
            0,
            [1, 2, 3, 4],
            model_part.GetProperties()[1])

        model_part.CreateNewElement(
            "Element2D4N",
            1,
            [2, 5, 6, 3],
            model_part.GetProperties()[1])

        return model_part


    def test_PointSetOutputProcess(self):
        model = KratosMultiphysics.Model()
        model_part = self.MakeModelPart(model)

        parameters = KratosMultiphysics.Parameters("""{
            "model_part_name" : "main",
            "positions" : [
                [0.0, 0.0, 0.0],
                [0.5, 0.5, 0.0],
                [1.0, 1.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.5, 0.5, 0.0],
                [2.0, 1.0, 0.0]
            ],
            "output_variables" : ["DISPLACEMENT", "REACTION"],
            "file_path" : "test.h5"
        }""")

        point_set_output_process = PointSetOuputProcessFactory(parameters, model)
        point_set_output_process.ExecuteInitialize()
        point_set_output_process.ExecuteFinalizeSolutionStep()





if __name__ == "__main__":
    UnitTest.main()