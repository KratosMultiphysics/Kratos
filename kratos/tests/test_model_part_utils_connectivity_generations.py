import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import os

def GetFilePath(fileName):
    """
    Constructs the absolute path for a given file name, based on the directory of the current script.

    Parameters:
    - fileName (str): The name of the file for which the absolute path is to be constructed.

    Returns:
    - str: The absolute path to the given file.
    """
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestModelPartUtilsConnectivityGenerations(KratosUnittest.TestCase):
    """
    A test class for verifying the functionalities of model part utilities in KratosMultiphysics.

    Attributes:
    - model (KratosMultiphysics.Model): The Kratos model object.
    - main_model_part (KratosMultiphysics.ModelPart): The main model part used for testing.
    - connectivities (list): A list to hold the connectivities of conditions for verification.
    """

    def setUp(self):
        """
        Set up the test environment before each test case.

        Initializes the model, main model part, and constructs a basic model part structure with nodes and conditions.

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        self.model = KratosMultiphysics.Model()
        self.main_model_part = self.model.CreateModelPart("Main", 2)

        # Model part structure:
        # 7 - 8 - 9
        # | \ | \ |
        # 4 - 5 - 6
        # | \ | \ |
        # 1 - 2 - 3

        # Populate model part with nodes
        for id, (x, y, z) in enumerate([(0,0,0), (1,0,0), (2,0,0), (0,1,0), (1,1,0), (2,1,0), (0,2,0), (1,2,0), (2,2,0)], start=1):
            self.main_model_part.CreateNewNode(id, x, y, z)

        # Create inlet conditions with properties
        prop_0 = self.main_model_part.GetProperties()[0]
        connectivities = [
            (1, 4, 2), (4, 5, 2), (2, 5, 3), (5, 6, 3),
            (4, 7, 5), (7, 8, 5), (6, 5, 8), (8, 9, 6)
        ]
        for i, nodes_ids in enumerate(connectivities, start=1):
            self.main_model_part.CreateNewCondition("SurfaceCondition3D3N", i, list(nodes_ids), prop_0)

        self.connectivities = [list(nodes_ids) for nodes_ids in connectivities]

    def test_ModelPartUtils_FromConnectivityGenerateElements(self):
        """
        Test the generation of elements from connectivities using ModelPartUtils.

        Verifies if elements are correctly generated and their connectivity matches the predefined connectivities.

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        # Generate elements
        KratosMultiphysics.ModelPartUtils.FromConnectivityGenerateElements(
            "Element3D3N",
            self.connectivities,
            self.main_model_part.Nodes,
            self.main_model_part.Elements,
            None
        )

        # Check connectivity of elements
        for elem in self.main_model_part.Elements:
            ref = [node.Id for node in elem.GetGeometry()]
            self.assertTrue(ref in self.connectivities)

    def test_ModelPartUtils_FromConnectivityGenerateConditions(self):
        """
        Test the generation of conditions from connectivities using ModelPartUtils.

        Verifies if conditions are correctly generated and their connectivity matches the predefined connectivities.

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        # Generate conditions
        KratosMultiphysics.ModelPartUtils.FromConnectivityGenerateConditions(
            "SurfaceCondition3D3N",
            self.connectivities,
            self.main_model_part.Nodes,
            self.main_model_part.Conditions,
            None
        )

        # Check connectivity of conditions
        for cond in self.main_model_part.Conditions:
            ref = [node.Id for node in cond.GetGeometry()]
            self.assertTrue(ref in self.connectivities)

if __name__ == '__main__':
    # Set logging severity to warning to minimize log output during tests
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
