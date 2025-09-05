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

    def testSetValues(self):
        """
        Tests the functionality related to solution step values and node values within a ModelPart.
        This test covers:
        - Checking the existence of specific solution step variables.
        - Verifying the correct assignment and retrieval of these variables.

        The test assigns PRESSURE and TEMPERATURE based on the X coordinate of each node,
        then verifies the presence and values of solution step variables and node values.
        """
        # Initialize the Kratos Multiphysics model
        current_model = KratosMultiphysics.Model()
        test_model_part = self.__SetUpTestModelPart(current_model)

        # Assign PRESSURE as a solution step value and TEMPERATURE as a node value
        # based on the node's X coordinate. Verify these assignments later.
        for node in test_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.PRESSURE, node.X)
            node.SetValue(KratosMultiphysics.TEMPERATURE, node.X)

        # Verify the presence of solution step variables and their correct assignment for each node
        for node in test_model_part.Nodes:
            # Check for the presence of DISPLACEMENT and its components as solution step variables
            self.assertTrue(node.HasSolutionStepValue(KratosMultiphysics.DISPLACEMENT))
            self.assertFalse(node.Has(KratosMultiphysics.DISPLACEMENT))  # Verify DISPLACEMENT is not a regular node value
            self.assertTrue(node.HasSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X))
            self.assertFalse(node.Has(KratosMultiphysics.DISPLACEMENT_X))
            self.assertTrue(node.HasSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y))
            self.assertFalse(node.Has(KratosMultiphysics.DISPLACEMENT_Y))
            self.assertTrue(node.HasSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z))
            self.assertFalse(node.Has(KratosMultiphysics.DISPLACEMENT_Z))

            # Check for the presence of REACTION and its components as solution step variables
            self.assertTrue(node.HasSolutionStepValue(KratosMultiphysics.REACTION))
            self.assertFalse(node.Has(KratosMultiphysics.REACTION))  # Verify REACTION is not a regular node value
            self.assertTrue(node.HasSolutionStepValue(KratosMultiphysics.REACTION_X))
            self.assertFalse(node.Has(KratosMultiphysics.REACTION_X))
            self.assertTrue(node.HasSolutionStepValue(KratosMultiphysics.REACTION_Y))
            self.assertFalse(node.Has(KratosMultiphysics.REACTION_Y))
            self.assertTrue(node.HasSolutionStepValue(KratosMultiphysics.REACTION_Z))
            self.assertFalse(node.Has(KratosMultiphysics.REACTION_Z))

            # Verify PRESSURE is correctly assigned as a solution step variable
            self.assertTrue(node.HasSolutionStepValue(KratosMultiphysics.PRESSURE))
            self.assertFalse(node.Has(KratosMultiphysics.PRESSURE))  # Verify PRESSURE is not a regular node value

            # Verify TEMPERATURE is not set as a solution step variable but as a regular node value
            self.assertFalse(node.HasSolutionStepValue(KratosMultiphysics.TEMPERATURE))
            self.assertTrue(node.Has(KratosMultiphysics.TEMPERATURE))

        # Verify the assigned values of PRESSURE and TEMPERATURE for each node
        for node in test_model_part.Nodes:
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.PRESSURE), node.X, msg="PRESSURE value mismatch")
            self.assertAlmostEqual(node.GetValue(KratosMultiphysics.TEMPERATURE), node.X, msg="TEMPERATURE value mismatch")

if __name__ == '__main__':
    # Set logging severity to WARNING to reduce clutter during test execution.
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
