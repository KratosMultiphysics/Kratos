import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Import assign_nodal_elements_to_nodes_process
import KratosMultiphysics.StructuralMechanicsApplication.assign_nodal_elements_to_nodes_process as assign_nodal_elements_to_nodes_process

class TestAssignNodalElementsToNodesProcess(KratosUnittest.TestCase):
    """
    This test suite checks the correct work of the AssignNodalElementsToNodesProcess process.
    """

    def setUp(self):
        """
        Code here will be placed BEFORE every test in this TestCase.
        """
        self.model = KratosMultiphysics.Model()
        self.model_part = self.model.CreateModelPart("Main")

    def AssignNodalElementsToNodesProcessCreateModelPart(self, model_part):
        """
        Helper function to create the model part with nodes.
        """
        # Create nodes
        model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
        model_part.CreateNewNode(4, 0.0, 1.0, 0.0)
        model_part.CreateNewNode(5, 2.0, 0.0, 0.0)
        model_part.CreateNewNode(6, 2.0, 1.0, 0.0)

    def test_AssignNodalElementsToNodesProcess1(self):
        """
        Test 1. No rotation
        """
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        self.AssignNodalElementsToNodesProcessCreateModelPart(self.model_part)

        parameters = KratosMultiphysics.Parameters("""
        {
            "model_part_name" : "Main",
            "nodal_mass"      : 1.0,
            "nodal_stiffness" : [1.0, null, 0.0]
        }""")

        process = StructuralMechanicsApplication.AssignNodalElementsToNodesProcess(self.model, parameters)
        process.Execute()

        self.assertEqual(len(self.model_part.Elements), len(self.model_part.Nodes))

        for elem in self.model_part.Elements:
            self.assertEqual(elem.Properties[KratosMultiphysics.NODAL_MASS], 1.0)
            self.assertEqual(elem.Properties[StructuralMechanicsApplication.NODAL_DISPLACEMENT_STIFFNESS][0], 1.0)
            self.assertEqual(elem.Properties[StructuralMechanicsApplication.NODAL_DISPLACEMENT_STIFFNESS][1], 0.0)
            self.assertEqual(elem.Properties[StructuralMechanicsApplication.NODAL_DISPLACEMENT_STIFFNESS][2], 0.0)

    def test_AssignNodalElementsToNodesProcess2(self):
        """
        Test 2. With rotation
        """
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ANGULAR_ACCELERATION)
        self.AssignNodalElementsToNodesProcessCreateModelPart(self.model_part)

        parameters = KratosMultiphysics.Parameters("""
        {
            "model_part_name"            : "Main",
            "nodal_mass"                 : 1.0,
            "nodal_inertia"              : [1.0, 1.0, 1.0],
            "nodal_stiffness"            : [1.0, null, 0.0],
            "nodal_rotational_stiffness" : [1.0, null, 1.0]
        }""")

        process = StructuralMechanicsApplication.AssignNodalElementsToNodesProcess(self.model, parameters)
        process.Execute()

        self.assertEqual(len(self.model_part.Elements), len(self.model_part.Nodes))

        for elem in self.model_part.Elements:
            self.assertEqual(elem.Properties[KratosMultiphysics.NODAL_MASS], 1.0)
            self.assertEqual(elem.Properties[StructuralMechanicsApplication.NODAL_DISPLACEMENT_STIFFNESS][0], 1.0)
            self.assertEqual(elem.Properties[StructuralMechanicsApplication.NODAL_DISPLACEMENT_STIFFNESS][1], 0.0)
            self.assertEqual(elem.Properties[StructuralMechanicsApplication.NODAL_DISPLACEMENT_STIFFNESS][2], 0.0)
            self.assertEqual(elem.Properties[StructuralMechanicsApplication.NODAL_INERTIA][0], 1.0)
            self.assertEqual(elem.Properties[StructuralMechanicsApplication.NODAL_INERTIA][1], 1.0)
            self.assertEqual(elem.Properties[StructuralMechanicsApplication.NODAL_INERTIA][2], 1.0)
            self.assertEqual(elem.Properties[StructuralMechanicsApplication.NODAL_ROTATIONAL_STIFFNESS][0], 1.0)
            self.assertEqual(elem.Properties[StructuralMechanicsApplication.NODAL_ROTATIONAL_STIFFNESS][1], 0.0)
            self.assertEqual(elem.Properties[StructuralMechanicsApplication.NODAL_ROTATIONAL_STIFFNESS][2], 1.0)

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()