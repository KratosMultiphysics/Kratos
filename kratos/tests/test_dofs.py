import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics

class TestDofs(KratosUnittest.TestCase):

    def __SetUpTestModelPart(self, model):
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
        current_model = KratosMultiphysics.Model()
        test_model_part = self.__SetUpTestModelPart(current_model)

        n = test_model_part.GetNode(1)
        n2 = test_model_part.GetNode(2)
        p = n.GetDof(KratosMultiphysics.PRESSURE)
        dx = n.GetDof(KratosMultiphysics.DISPLACEMENT_X)
        dy = n.GetDof(KratosMultiphysics.DISPLACEMENT_Y)
        dz = n.GetDof(KratosMultiphysics.DISPLACEMENT_Z)
        p2 = n2.GetDof(KratosMultiphysics.PRESSURE)

        n.Fix(KratosMultiphysics.DISPLACEMENT_X)
        dx.EquationId = 5
        dy.EquationId = 6
        dz.EquationId = 7

        self.assertEqual(p.GetVariable(), KratosMultiphysics.PRESSURE)
        self.assertEqual(dx.GetVariable(), KratosMultiphysics.DISPLACEMENT_X)
        self.assertEqual(dx.GetReaction(), KratosMultiphysics.REACTION_X)
        self.assertEqual(dx.EquationId, 5)

        self.assertLess(p,p2)
        self.assertGreater(dy,dx)
        self.assertGreater(p2,dz)

    def testDofListValues(self):
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
        current_model = KratosMultiphysics.Model()
        test_model_part = self.__SetUpTestModelPart(current_model)

        dofs_vector = KratosMultiphysics.DofsArrayType()
        for node in test_model_part.Nodes:
            dofs_vector.append(node.GetDof(KratosMultiphysics.PRESSURE))

        self.assertEqual(len(dofs_vector), test_model_part.NumberOfNodes())

    def testDofListUnique(self):
        current_model = KratosMultiphysics.Model()
        test_model_part = self.__SetUpTestModelPart(current_model)

        dofs_vector = KratosMultiphysics.DofsArrayType()
        for node in test_model_part.Nodes:
            dofs_vector.append(node.GetDof(KratosMultiphysics.PRESSURE))
            dofs_vector.append(node.GetDof(KratosMultiphysics.PRESSURE))
            dofs_vector.append(node.GetDof(KratosMultiphysics.PRESSURE))

        self.assertEqual(len(dofs_vector), 3*test_model_part.NumberOfNodes())
        dofs_vector.unique()
        self.assertEqual(len(dofs_vector), test_model_part.NumberOfNodes())

if __name__ == '__main__':
    KratosUnittest.main()
