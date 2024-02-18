import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics

class TestDofs(KratosUnittest.TestCase):
    """
    This module contains unit tests for KratosMultiphysics DOF-related functionality.
    It covers DOF creation, manipulation, and DOF lists.
    """

    def __SetUpTestModelPart(self, model):
        """
        Set up a test ModelPart with nodes and DOFs.

        Args:
            model (KratosMultiphysics.Model): The Kratos model to create the ModelPart in.

        Returns:
            KratosMultiphysics.ModelPart: The created ModelPart with nodes and DOFs.
        """
        mp = model.CreateModelPart("Main")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)

        n1 = mp.CreateNewNode(1,0.0,0.0,0.0)
        n1.AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X)
        n1.AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y)
        n1.AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z)
        n1.AddDof(KratosMultiphysics.PRESSURE)

        n2 = mp.CreateNewNode(2,2.0,2.0,2.0)
        n2.AddDof(KratosMultiphysics.PRESSURE)

        n3 = mp.CreateNewNode(3,3.0,3.0,3.0)
        n3.AddDof(KratosMultiphysics.PRESSURE)

        return mp

    def testDofExport(self):
        """
        Test the export of DOF information.

        This test creates a ModelPart with nodes and DOFs, exports and manipulates DOFs,
        and performs various assertions.

        """
        current_model = KratosMultiphysics.Model()
        test_model_part = self.__SetUpTestModelPart(current_model)

        n = test_model_part.GetNode(1)
        n2 = test_model_part.GetNode(2)
        p = n.GetDof(KratosMultiphysics.PRESSURE)
        dx = n.GetDof(KratosMultiphysics.DISPLACEMENT_X)
        dy = n.GetDof(KratosMultiphysics.DISPLACEMENT_Y)
        dz = n.GetDof(KratosMultiphysics.DISPLACEMENT_Z)
        p2 = n2.GetDof(KratosMultiphysics.PRESSURE)

        # Fixing only displacement in X
        n.Fix(KratosMultiphysics.DISPLACEMENT_X)
        self.assertTrue(n.IsFixed(KratosMultiphysics.DISPLACEMENT_X))
        self.assertFalse(n.IsFixed(KratosMultiphysics.DISPLACEMENT_Y))
        self.assertFalse(n.IsFixed(KratosMultiphysics.DISPLACEMENT_Z))
        self.assertFalse(n.IsFixed(KratosMultiphysics.PRESSURE))

        # Release the displacement in X and checking again
        n.Free(KratosMultiphysics.DISPLACEMENT_X)
        self.assertFalse(n.IsFixed(KratosMultiphysics.DISPLACEMENT_X))
        self.assertFalse(n.IsFixed(KratosMultiphysics.DISPLACEMENT_Y))
        self.assertFalse(n.IsFixed(KratosMultiphysics.DISPLACEMENT_Z))
        self.assertFalse(n.IsFixed(KratosMultiphysics.PRESSURE))

        # Fixing again dof
        n.Fix(KratosMultiphysics.DISPLACEMENT_X)

        # Assign equation Id
        dx.EquationId = 5
        dy.EquationId = 6
        dz.EquationId = 7

        # Checks
        self.assertEqual(p.GetVariable(), KratosMultiphysics.PRESSURE)
        self.assertEqual(dx.GetVariable(), KratosMultiphysics.DISPLACEMENT_X)
        self.assertEqual(dx.GetReaction(), KratosMultiphysics.REACTION_X)
        self.assertEqual(dx.EquationId, 5)

        self.assertLess(p,p2)
        self.assertGreater(dy,dx)
        self.assertGreater(p2,dz)

    def testDofListValues(self):
        """
        Test DOF list values and manipulation.

        This test creates a ModelPart with nodes and DOFs, creates a DOF list, sets and
        retrieves values, and performs assertions.

        """
        current_model = KratosMultiphysics.Model()
        test_model_part = self.__SetUpTestModelPart(current_model)

        n = test_model_part.GetNode(1)
        p = n.GetDof(KratosMultiphysics.PRESSURE)
        dx = n.GetDof(KratosMultiphysics.DISPLACEMENT_X)
        dy = n.GetDof(KratosMultiphysics.DISPLACEMENT_Y)
        dz = n.GetDof(KratosMultiphysics.DISPLACEMENT_Z)

        ##dof list
        dofs = KratosMultiphysics.DofsArrayType()
        dofs.append(p)
        dofs.append(dx)
        dofs.append(dy)
        dofs.append(dz)

        counter = 0
        for dof in dofs:
            dof.EquationId = counter
            dof.SetSolutionStepValue(counter)
            counter += 1

        values = dofs.GetValues()
        counter = 0
        for value in values:
            self.assertEqual(value, counter)
            counter += 1

        values[0] = 14.0 #this is the pressure
        dofs.SetValues(values)
        self.assertEqual(p.GetSolutionStepValue(), 14.0)

    def testDofListAppend(self):
        """
        Test appending DOFs to a list.

        This test creates a ModelPart with nodes and DOFs and appends DOFs to a list,
        asserting the length of the list.

        """
        current_model = KratosMultiphysics.Model()
        test_model_part = self.__SetUpTestModelPart(current_model)

        dofs_vector = KratosMultiphysics.DofsArrayType()
        for node in test_model_part.Nodes:
            dofs_vector.append(node.GetDof(KratosMultiphysics.PRESSURE))

        self.assertEqual(len(dofs_vector), test_model_part.NumberOfNodes())

    def testDofListUnique(self):
        """
        Test making a DOF list unique.

        This test creates a ModelPart with nodes and DOFs, appends multiple copies of
        the same DOF to a list, makes the list unique, and asserts its length.

        """
        current_model = KratosMultiphysics.Model()
        test_model_part = self.__SetUpTestModelPart(current_model)

        dofs_vector = KratosMultiphysics.DofsArrayType()
        for node in test_model_part.Nodes:
            dofs_vector.append(node.GetDof(KratosMultiphysics.PRESSURE))
            dofs_vector.append(node.GetDof(KratosMultiphysics.PRESSURE))
            dofs_vector.append(node.GetDof(KratosMultiphysics.PRESSURE))

        self.assertEqual(len(dofs_vector), test_model_part.NumberOfNodes())

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
