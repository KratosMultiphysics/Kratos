# Import necessary modules from KratosMultiphysics
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.modelers.import_obj_modeler import ImportOBJModeler

class TestImportOBJModeler(KratosUnittest.TestCase):
    """
    Unit tests for the ImportOBJModeler class.
    """

    def setUp(self):
        """
        Sets up the test environment.
        """
        self.model = KratosMultiphysics.Model()
        self.settings = KratosMultiphysics.Parameters('''{
            "input_filename"  : "auxiliar_files_for_python_unittest/obj_files/cube.obj",
            "model_part_name" : "Main"
        }''')

    def test_ImportOBJModeler(self):
        """
        Tests the SetupModelPart method.
        """
        modeler = ImportOBJModeler(self.model, self.settings)
        modeler.SetupModelPart()

        # Retrieve model part
        model_part_name = self.settings["model_part_name"].GetString()
        self.assertTrue(self.model.HasModelPart(model_part_name))
        self.model_part = self.model.GetModelPart(model_part_name)

        # Assert that the model part has the correct number of nodes and elements
        self.assertEqual(self.model_part.NumberOfNodes(), 8)
        self.assertEqual(self.model_part.NumberOfElements(), 6)

        # Assert that the normals are the same as the nodes coordinates
        for node in self.model_part.Nodes:
            normal = node.GetValue(KratosMultiphysics.NORMAL)
            self.assertAlmostEqual(normal[0], node.X)
            self.assertAlmostEqual(normal[1], node.Y)
            self.assertAlmostEqual(normal[2], node.Z)

if __name__ == '__main__':
    # Set the logger severity level and run the KratosUnittest
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()