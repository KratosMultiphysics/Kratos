from __future__ import print_function, absolute_import, division

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics

class TestDofs(KratosUnittest.TestCase):
    
    def test_model_part_sub_model_parts(self):
        current_model = KratosMultiphysics.Model()

        mp = current_model.CreateModelPart("Main")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)

        n = mp.CreateNewNode(1,0.0,0.0,0.0)
        n.AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X)
        n.AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y)
        n.AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z)
        n.AddDof(KratosMultiphysics.PRESSURE)
        n2 = mp.CreateNewNode(2,2.0,2.0,2.0)
        n2.AddDof(KratosMultiphysics.PRESSURE)

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

        #comparison
        self.assertTrue(p<p2)
        self.assertTrue(dy>dx)
        self.assertTrue(p2>dz)



if __name__ == '__main__':
    KratosUnittest.main()
