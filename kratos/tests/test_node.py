import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestNode(KratosUnittest.TestCase):
    """
    TestNode class extends KratosUnittest.TestCase to provide unit tests for the functionality
    of KratosMultiphysics nodes. It includes tests for checking the presence of nodal solution step variables
    and their correct assignment and retrieval.

    Methods:
        __SetUpTestModelPart(self, model): Sets up a test environment with a ModelPart containing nodes.
        testSolutionStepValue(self): Validates the assignment and retrieval of solution step values at nodes.
    """

    def __SetUpTestModelPart(self, model):
        """
        Private method to set up a ModelPart for testing. This method initializes a ModelPart
        with a specific configuration of nodes and assigns them solution step variables.

        Args:
            model (KratosMultiphysics.Model): The Kratos model used to create the ModelPart.

        Returns:
            KratosMultiphysics.ModelPart: A configured ModelPart with nodes and predefined solution step variables.
        """
        mp = model.CreateModelPart("Main")
        # Adding necessary nodal solution step variables to the model part.
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)

        # Creating nodes with specific IDs and coordinates.
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, 2.0, 2.0, 2.0)
        mp.CreateNewNode(3, 3.0, 3.0, 3.0)

        return mp

    def testSolutionStepValue(self):
        """
        Tests the functionality related to solution step values of nodes within a ModelPart.
        This includes checking the existence of specific solution step variables and
        verifying their correct assignment and retrieval.
        """
        current_model = KratosMultiphysics.Model()
        test_model_part = self.__SetUpTestModelPart(current_model)

        # Verify the presence of solution step variables for each node.
        for node in test_model_part.Nodes:
            self.assertTrue(node.HasSolutionStepValue(KratosMultiphysics.DISPLACEMENT))
            self.assertTrue(node.HasSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X))
            self.assertTrue(node.HasSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y))
            self.assertTrue(node.HasSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z))
            self.assertTrue(node.HasSolutionStepValue(KratosMultiphysics.REACTION))
            self.assertTrue(node.HasSolutionStepValue(KratosMultiphysics.REACTION_X))
            self.assertTrue(node.HasSolutionStepValue(KratosMultiphysics.REACTION_Y))
            self.assertTrue(node.HasSolutionStepValue(KratosMultiphysics.REACTION_Z))
            self.assertTrue(node.HasSolutionStepValue(KratosMultiphysics.PRESSURE))
            self.assertFalse(node.HasSolutionStepValue(KratosMultiphysics.TEMPERATURE))

        # Assign PRESSURE solution step value based on node X coordinate and verify.
        for node in test_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.PRESSURE, node.X)
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.PRESSURE), node.X)

if __name__ == '__main__':
    # Set logging severity to WARNING to reduce clutter during test execution.
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
