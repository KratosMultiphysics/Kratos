import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.MeshioPlusPlusApplication.modelers.meshio_interpolate_modeler import Factory as MeshioInterpolateModelerFactory


def _CreateCubeOfTetrahedra(model_part):
    """A unit cube split into six linear tetrahedra."""
    model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
    model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
    model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
    model_part.CreateNewNode(4, 0.0, 1.0, 0.0)
    model_part.CreateNewNode(5, 0.0, 0.0, 1.0)
    model_part.CreateNewNode(6, 1.0, 0.0, 1.0)
    model_part.CreateNewNode(7, 1.0, 1.0, 1.0)
    model_part.CreateNewNode(8, 0.0, 1.0, 1.0)

    properties = model_part.CreateNewProperties(1)
    connectivities = [
        [1, 2, 3, 7], [1, 2, 7, 6], [1, 6, 7, 5],
        [1, 3, 4, 7], [1, 4, 8, 7], [1, 8, 5, 7],
    ]
    for i, nodes in enumerate(connectivities):
        model_part.CreateNewElement("Element3D4N", i + 1, nodes, properties)


class TestMeshioInterpolateModeler(KratosUnittest.TestCase):
    """Tests for the meshio++ interpolate modeler."""

    def setUp(self):
        self.model = KratosMultiphysics.Model()
        self.source = self.model.CreateModelPart("Source")
        _CreateCubeOfTetrahedra(self.source)
        for node in self.source.Nodes:
            node.SetValue(KratosMultiphysics.TEMPERATURE, node.Z)
        self.target = self.model.CreateModelPart("Target")
        _CreateCubeOfTetrahedra(self.target)

    def testInterpolateOntoAnIdenticalTarget(self):
        # Source and target are geometrically identical, so nearest sampling is exact.
        settings = KratosMultiphysics.Parameters("""{
            "source_model_part_name" : "Source",
            "target_model_part_name" : "Target",
            "output_model_part_name" : "Interpolated",
            "operation_settings" : {
                "method" : "nearest",
                "on_conflict" : "overwrite",
                "nodal_data_value_variables" : ["TEMPERATURE"]
            }
        }""")
        modeler = MeshioInterpolateModelerFactory(self.model, settings)

        # The destination model part must exist right after construction so solvers can
        # add nodal variables before the mesh is filled
        self.assertTrue(self.model.HasModelPart("Interpolated"))
        self.assertEqual(self.model["Interpolated"].NumberOfNodes(), 0)

        modeler.SetupGeometryModel()
        modeler.PrepareGeometryModel()
        modeler.SetupModelPart()

        destination = self.model["Interpolated"]
        self.assertEqual(destination.NumberOfNodes(), self.target.NumberOfNodes())
        for node in destination.Nodes:
            self.assertAlmostEqual(node.GetValue(KratosMultiphysics.TEMPERATURE), node.Z, places=10)

    def testMissingSourceModelPartNameRaises(self):
        settings = KratosMultiphysics.Parameters("""{
            "target_model_part_name" : "Target",
            "output_model_part_name" : "Interpolated"
        }""")
        with self.assertRaisesRegex(Exception, "source_model_part_name"):
            MeshioInterpolateModelerFactory(self.model, settings)

    def testMissingTargetModelPartNameRaises(self):
        settings = KratosMultiphysics.Parameters("""{
            "source_model_part_name" : "Source",
            "output_model_part_name" : "Interpolated"
        }""")
        with self.assertRaisesRegex(Exception, "target_model_part_name"):
            MeshioInterpolateModelerFactory(self.model, settings)

    def testModelerIsRegistered(self):
        self.assertTrue(KratosMultiphysics.Registry.HasItem(
            "Modelers.KratosMultiphysics.MeshioPlusPlusApplication.MeshioInterpolateModeler"))


if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
