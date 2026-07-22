import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest


def _CreateBoxSourceModelPart(model, name, x_min, y_min, z_min, x_max, y_max, z_max):
    """Create a model part containing only the 8 corner nodes of an axis-aligned box."""
    mp = model.CreateModelPart(name)
    mp.CreateNewNode(1, x_min, y_min, z_min)
    mp.CreateNewNode(2, x_max, y_min, z_min)
    mp.CreateNewNode(3, x_max, y_max, z_min)
    mp.CreateNewNode(4, x_min, y_max, z_min)
    mp.CreateNewNode(5, x_min, y_min, z_max)
    mp.CreateNewNode(6, x_max, y_min, z_max)
    mp.CreateNewNode(7, x_max, y_max, z_max)
    mp.CreateNewNode(8, x_min, y_max, z_max)
    return mp


class TestCartesianMeshGeneratorModeler(KratosUnittest.TestCase):

    # ------------------------------------------------------------------
    # Direct API (ModelPart + element_size constructor)
    # ------------------------------------------------------------------

    def test_3d_node_count(self):
        """1x1x1 box, element_size=0.5 → 3^3=27 nodes, 2^3*6=48 tets."""
        model = KratosMultiphysics.Model()
        source = _CreateBoxSourceModelPart(model, "source", 0.0, 0.0, 0.0, 1.0, 1.0, 1.0)
        output = model.CreateModelPart("output")
        output.CreateNewProperties(0)

        modeler = KratosMultiphysics.CartesianMeshGeneratorModeler(source, 0.5)
        modeler.GenerateMesh(output, "Element3D4N")

        self.assertEqual(output.NumberOfNodes(), 27)
        self.assertEqual(output.NumberOfElements(), 48)

    def test_3d_nodes_inside_bounding_box(self):
        """All generated nodes must lie inside (or on the boundary of) the source BB."""
        model = KratosMultiphysics.Model()
        source = _CreateBoxSourceModelPart(model, "source", -1.0, -1.0, -1.0, 1.0, 1.0, 1.0)
        output = model.CreateModelPart("output")
        output.CreateNewProperties(0)

        modeler = KratosMultiphysics.CartesianMeshGeneratorModeler(source, 0.5)
        modeler.GenerateMesh(output, "Element3D4N")

        tol = 1e-12
        for node in output.Nodes:
            self.assertGreaterEqual(node.X, -1.0 - tol)
            self.assertLessEqual(node.X,     1.0 + tol)
            self.assertGreaterEqual(node.Y, -1.0 - tol)
            self.assertLessEqual(node.Y,     1.0 + tol)
            self.assertGreaterEqual(node.Z, -1.0 - tol)
            self.assertLessEqual(node.Z,     1.0 + tol)

    def test_3d_element_connectivity(self):
        """Each element must have 4 nodes that exist in the output model part."""
        model = KratosMultiphysics.Model()
        source = _CreateBoxSourceModelPart(model, "source", 0.0, 0.0, 0.0, 1.0, 1.0, 1.0)
        output = model.CreateModelPart("output")
        output.CreateNewProperties(0)

        # element_size = 1.0 → 1 cell per direction → 6 tetrahedra
        modeler = KratosMultiphysics.CartesianMeshGeneratorModeler(source, 1.0)
        modeler.GenerateMesh(output, "Element3D4N")

        self.assertEqual(output.NumberOfElements(), 6)
        node_ids = {node.Id for node in output.Nodes}
        for elem in output.Elements:
            self.assertEqual(len(elem.GetNodes()), 4)
            for node in elem.GetNodes():
                self.assertIn(node.Id, node_ids)

    def test_3d_non_unit_box(self):
        """3x2x1 box, element_size=0.5 → 7*5*3=105 nodes, 6*4*2*6=288 tets."""
        model = KratosMultiphysics.Model()
        source = _CreateBoxSourceModelPart(model, "source", 0.0, 0.0, 0.0, 3.0, 2.0, 1.0)
        output = model.CreateModelPart("output")
        output.CreateNewProperties(0)

        modeler = KratosMultiphysics.CartesianMeshGeneratorModeler(source, 0.5)
        modeler.GenerateMesh(output, "Element3D4N")

        self.assertEqual(output.NumberOfNodes(), 105)
        self.assertEqual(output.NumberOfElements(), 288)

    def test_3d_uniform_spacing(self):
        """Node coordinates should be uniformly spaced by element_size."""
        model = KratosMultiphysics.Model()
        source = _CreateBoxSourceModelPart(model, "source", 0.0, 0.0, 0.0, 1.0, 1.0, 1.0)
        output = model.CreateModelPart("output")
        output.CreateNewProperties(0)

        element_size = 0.5
        modeler = KratosMultiphysics.CartesianMeshGeneratorModeler(source, element_size)
        modeler.GenerateMesh(output, "Element3D4N")

        # Collect unique X values and check spacing
        xs = sorted({round(node.X, 10) for node in output.Nodes})
        for i in range(1, len(xs)):
            self.assertAlmostEqual(xs[i] - xs[i-1], element_size, places=10)

    # ------------------------------------------------------------------
    # Factory / registry constructor (Model + Parameters)
    # ------------------------------------------------------------------

    def test_setup_model_part(self):
        """SetupModelPart() via factory constructor must produce the same mesh."""
        model = KratosMultiphysics.Model()
        _CreateBoxSourceModelPart(model, "source", 0.0, 0.0, 0.0, 2.0, 2.0, 2.0)

        params = KratosMultiphysics.Parameters("""{
            "input_model_part_name"  : "source",
            "output_model_part_name" : "output",
            "element_name"           : "Element3D4N",
            "element_size"           : 1.0
        }""")

        modeler = KratosMultiphysics.CartesianMeshGeneratorModeler(model, params)
        modeler.SetupModelPart()

        output = model.GetModelPart("output")
        # 2x2x2 box, element_size=1 → 3^3=27 nodes, 2^3*6=48 tets
        self.assertEqual(output.NumberOfNodes(), 27)
        self.assertEqual(output.NumberOfElements(), 48)

    def test_setup_model_part_creates_output(self):
        """SetupModelPart() must create the output model part if it does not exist."""
        model = KratosMultiphysics.Model()
        _CreateBoxSourceModelPart(model, "box", 0.0, 0.0, 0.0, 1.0, 1.0, 1.0)

        params = KratosMultiphysics.Parameters("""{
            "input_model_part_name"  : "box",
            "output_model_part_name" : "mesh",
            "element_name"           : "Element3D4N",
            "element_size"           : 1.0
        }""")

        KratosMultiphysics.CartesianMeshGeneratorModeler(model, params).SetupModelPart()
        self.assertTrue(model.HasModelPart("mesh"))

    def test_missing_input_name_raises(self):
        """Empty input_model_part_name must raise an error."""
        model = KratosMultiphysics.Model()
        model.CreateModelPart("source")

        params = KratosMultiphysics.Parameters("""{
            "input_model_part_name"  : "",
            "output_model_part_name" : "output",
            "element_name"           : "Element3D4N",
            "element_size"           : 1.0
        }""")

        modeler = KratosMultiphysics.CartesianMeshGeneratorModeler(model, params)
        with self.assertRaises(RuntimeError):
            modeler.SetupModelPart()

    def test_missing_output_name_raises(self):
        """Empty output_model_part_name must raise an error."""
        model = KratosMultiphysics.Model()
        _CreateBoxSourceModelPart(model, "source", 0.0, 0.0, 0.0, 1.0, 1.0, 1.0)

        params = KratosMultiphysics.Parameters("""{
            "input_model_part_name"  : "source",
            "output_model_part_name" : "",
            "element_name"           : "Element3D4N",
            "element_size"           : 1.0
        }""")

        modeler = KratosMultiphysics.CartesianMeshGeneratorModeler(model, params)
        with self.assertRaises(RuntimeError):
            modeler.SetupModelPart()

    # ------------------------------------------------------------------
    # Info / string interface
    # ------------------------------------------------------------------

    def test_info(self):
        """Info() must return the class name."""
        model = KratosMultiphysics.Model()
        source = model.CreateModelPart("source")
        modeler = KratosMultiphysics.CartesianMeshGeneratorModeler(source, 1.0)
        self.assertEqual(str(modeler).strip(), "CartesianMeshGeneratorModeler")

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
