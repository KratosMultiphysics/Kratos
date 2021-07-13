# Kratos Imports
import KratosMultiphysics
from  KratosMultiphysics.MeshMovingApplication.impose_mesh_motion_process import ImposeMeshMotionProcess
import KratosMultiphysics.KratosUnittest as KratosUnittest


def GenerateModel():
    model = KratosMultiphysics.Model()
    model_part = model.CreateModelPart("Main")

    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_DISPLACEMENT)

    # Generate nodes
    domain = KratosMultiphysics.Quadrilateral2D4(
        KratosMultiphysics.Node(1, 0.0, 0.0, 0.0),
        KratosMultiphysics.Node(2, 0.0, 2.0, 0.0),
        KratosMultiphysics.Node(3, 2.0, 2.0, 0.0),
        KratosMultiphysics.Node(4, 2.0, 0.0, 0.0))

    parameters = KratosMultiphysics.Parameters("""{
        "element_name" : "Element2D3N",
        "number_of_divisions" : 2,
        "create_skin_sub_model_part" : false
    }""")

    KratosMultiphysics.StructuredMeshGeneratorProcess(
        domain,
        model_part,
        parameters
    ).Execute()

    # Reset current node positions
    # (mesh motion should use the initial coordinates)
    for node in model_part.Nodes:
        node.X = 0.0
        node.Y = 0.0
        node.Z = 0.0

    return model, model_part


class TestImposeMeshMotionProcess(KratosUnittest.TestCase):

    def CheckNodes(self, model_part: KratosMultiphysics.ModelPart):
        self.assertEqual(len(model_part.Nodes), 9)

        self.assertVectorAlmostEqual(model_part.GetNode(1).GetSolutionStepValue(KratosMultiphysics.MESH_DISPLACEMENT), [0.0, 2.0, 4.0])
        self.assertVectorAlmostEqual(model_part.GetNode(2).GetSolutionStepValue(KratosMultiphysics.MESH_DISPLACEMENT), [0.0, 2.0, 4.0])
        self.assertVectorAlmostEqual(model_part.GetNode(3).GetSolutionStepValue(KratosMultiphysics.MESH_DISPLACEMENT), [0.0, 2.0, 4.0])
        self.assertVectorAlmostEqual(model_part.GetNode(4).GetSolutionStepValue(KratosMultiphysics.MESH_DISPLACEMENT), [-1.0, 2.0, 5.0])
        self.assertVectorAlmostEqual(model_part.GetNode(5).GetSolutionStepValue(KratosMultiphysics.MESH_DISPLACEMENT), [-1.0, 2.0, 5.0])
        self.assertVectorAlmostEqual(model_part.GetNode(6).GetSolutionStepValue(KratosMultiphysics.MESH_DISPLACEMENT), [-1.0, 2.0, 5.0])
        self.assertVectorAlmostEqual(model_part.GetNode(7).GetSolutionStepValue(KratosMultiphysics.MESH_DISPLACEMENT), [-2.0, 2.0, 6.0])
        self.assertVectorAlmostEqual(model_part.GetNode(8).GetSolutionStepValue(KratosMultiphysics.MESH_DISPLACEMENT), [-2.0, 2.0, 6.0])
        self.assertVectorAlmostEqual(model_part.GetNode(9).GetSolutionStepValue(KratosMultiphysics.MESH_DISPLACEMENT), [-2.0, 2.0, 6.0])

    def test_axis_and_angle(self):
        model, model_part = GenerateModel()

        parameters = KratosMultiphysics.Parameters("""{
            "model_part_name"     : "Main",
            "rotation_definition" : "rotation_axis",
            "rotation_axis"       : [0, 1, 0],
            "reference_point"     : [-1, 0, 0],
            "rotation_angle"      : -1.57079632679,
            "translation_vector"  : [1, 2, 3]
        }""")

        ImposeMeshMotionProcess(
            model,
            parameters
        ).ExecuteInitializeSolutionStep()

        self.CheckNodes(model_part)

    def test_euler_angles(self):
        model, model_part = GenerateModel()

        parameters = KratosMultiphysics.Parameters("""{
            "model_part_name"     : "Main",
            "rotation_definition" : "euler_angles",
            "euler_angles"        : [-1.57079632679, -1.57079632679, 1.57079632679],
            "reference_point"     : [-1, 0, 0],
            "translation_vector"  : [1, 2, 3]
        }""")

        ImposeMeshMotionProcess(
            model,
            parameters
        ).ExecuteInitializeSolutionStep()

        self.CheckNodes(model_part)

    def test_Parametric(self):
        model, model_part = GenerateModel()

        parameters = KratosMultiphysics.Parameters("""{
            "model_part_name" : "Main",
            "rotation_definition" : "rotation_axis",
            "rotation_axis" : ["1.0-t", "t", 0.0],
            "rotation_angle" : "-1.57079632679 * t",
            "reference_point" : ["-t", "z", 0.0],
            "translation_vector" : [1, 2, 3]
        }""")

        model_part.ProcessInfo[KratosMultiphysics.TIME] = 1.0
        ImposeMeshMotionProcess(
            model,
            parameters
        ).ExecuteInitializeSolutionStep()

        self.CheckNodes(model_part)


if __name__ == "__main__":
    KratosUnittest.main()