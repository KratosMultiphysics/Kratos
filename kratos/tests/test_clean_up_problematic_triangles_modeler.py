# Import necessary modules from KratosMultiphysics
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestCleanUpProblematicTrianglesModeler(KratosUnittest.TestCase):
    """
    Define a test class for testing the CleanUpProblematicTrianglesModeler
    """
    @classmethod
    def setUpClass(cls):
        """
        Set up the test class before running tests.
        """
        pass

    @classmethod
    def tearDownClass(cls):
        """
        Clean up after running tests.
        """
        pass

    def test_CleanUpProblematicTrianglesModeler(self):
        """
        Test the CleanUpProblematicTrianglesModeler
        """
        # Create a model part and set the domain size
        self.current_model = KratosMultiphysics.Model()
        self.model_part = self.current_model.CreateModelPart("Main")
        self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 3

        # Nodes (ID, X, Y, Z)
        nodes = [
            (-0.5, -0.5, 0.5),
            (0.5, -0.5, 0.5),
            (0.5, 0.5, 0.5),
            (-0.5, 0.5, 0.5),
            (-0.5, -0.5, -0.5),
            (0.5, -0.5, -0.5),
            (0.5, 0.5, -0.5),
            (-0.5, 0.5, -0.5),
            (-0.5, -0.5, 0.5),  # Extra vertex from the face definitions
            (-0.5, 0.5, -0.5),  # Extra vertex from the face definitions
        ]

        # Create nodes
        id = 1
        for node in nodes:
            self.model_part.CreateNewNode(id, node[0], node[1], node[2])
            id += 1

        # Elements (quadrilateral and triangular surfaces)
        # Assuming all elements are surfaces with the 'element' ID and nodal connectivity.
        elements = [
            # Quadrilateral faces
            (1, 2, 3, 4),  # Face 1 (Vertices 1, 2, 3, 4)
            (8, 7, 6, 5),  # Face 2 (Vertices 8, 7, 6, 5)
            (4, 8, 5, 1),  # Face 3 (Vertices 4, 8, 5, 1)
            (2, 6, 7, 3),  # Face 4 (Vertices 2, 6, 7, 3)
            (5, 6, 2, 1),  # Face 5 (Vertices 5, 6, 2, 1)
            (4, 3, 7, 8),  # Face 6 (Vertices 4, 3, 7, 8)

            # Triangular faces (two extra triangles)
            (1, 2, 9),     # Face 7 (Vertices 1, 2, 9)
            (8, 7, 10)     # Face 8 (Vertices 8, 7, 10)
        ]

        # Create elements
        id = 1
        for element in elements:
            if len(element) == 3:
                self.model_part.CreateNewElement("SurfaceElement3D3N", id, element, self.model_part.GetProperties()[1])
            elif len(element) == 4:
                self.model_part.CreateNewElement("SurfaceElement3D4N", id, element, self.model_part.GetProperties()[1])
            id += 1

        clean_settings = KratosMultiphysics.Parameters("""{
            "model_part_name"    : "Main"
        }""")
        KratosMultiphysics.CleanUpProblematicTrianglesModeler(self.current_model, clean_settings).SetupModelPart()

        # Assert that the model part has the correct number of nodes and elements
        self.assertEqual(self.model_part.NumberOfNodes(), 8)
        self.assertEqual(self.model_part.NumberOfElements(), 6)

if __name__ == '__main__':
    # Set the logger severity level and run the KratosUnittest
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()