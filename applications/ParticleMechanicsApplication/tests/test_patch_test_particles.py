from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.ParticleMechanicsApplication as KratosParticle
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestLoadingConditionsPoint(KratosUnittest.TestCase):
    def _execute_point_load_condition_test(self, Dimension):
        mp = KratosMultiphysics.ModelPart("solid_part")
        mp.SetBufferSize(2)

        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_LOAD)

        # create node
        node = mp.CreateNewNode(1,0.0,0.0,0.0)

        # ensure that the property 1 is created
        mp.GetProperties()[1]

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)

        if Dimension == 2:
            cond = mp.CreateNewCondition("PointLoadCondition2D1N", 1, [1], mp.GetProperties()[1])
        elif Dimension == 3:
            cond = mp.CreateNewCondition("PointLoadCondition3D1N", 1, [1], mp.GetProperties()[1])
        else:
            raise RuntimeError("Wrong Dimension")

        lhs = KratosMultiphysics.Matrix(0,0)
        rhs = KratosMultiphysics.Vector(0)

        # first we apply a load to the condition
        load_on_cond = KratosMultiphysics.Vector(3)
        load_on_cond[0] = 1.8
        load_on_cond[1] = 2.6
        load_on_cond[2] = -11.47

        cond.SetValue(StructuralMechanicsApplication.POINT_LOAD, load_on_cond)

        cond.CalculateLocalSystem(lhs,rhs,mp.ProcessInfo)

        self.assertEqual(rhs[0], load_on_cond[0])
        self.assertEqual(rhs[1], load_on_cond[1])
        if Dimension == 3:
            self.assertEqual(rhs[2], load_on_cond[2])

        # now we apply a load to the node
        nodal_load = KratosMultiphysics.Vector(3)
        nodal_load[0] = -5.5
        nodal_load[1] = 1.2
        nodal_load[2] = 9.3

        node.SetSolutionStepValue(StructuralMechanicsApplication.POINT_LOAD, nodal_load)

        cond.CalculateLocalSystem(lhs,rhs,mp.ProcessInfo)

        self.assertEqual(rhs[0], load_on_cond[0] + nodal_load[0])
        self.assertEqual(rhs[1], load_on_cond[1] + nodal_load[1])
        if Dimension == 3:
            self.assertEqual(rhs[2], load_on_cond[2] + nodal_load[2])

    def _execute_point_moment_condition_test(self):
        mp = KratosMultiphysics.ModelPart("solid_part")
        mp.SetBufferSize(2)

        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_MOMENT)
        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_MOMENT)

        # create node
        node = mp.CreateNewNode(1,0.0,0.0,0.0)

        # ensure that the property 1 is created
        mp.GetProperties()[1]

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_X, KratosMultiphysics.REACTION_MOMENT_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Y, KratosMultiphysics.REACTION_MOMENT_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Z, KratosMultiphysics.REACTION_MOMENT_Z,mp)

        cond = mp.CreateNewCondition("PointMomentCondition3D1N", 1, [1], mp.GetProperties()[1])

        lhs = KratosMultiphysics.Matrix(0,0)
        rhs = KratosMultiphysics.Vector(0)

        # first we apply a load to the condition
        load_on_cond = KratosMultiphysics.Vector(3)
        load_on_cond[0] = 1.8
        load_on_cond[1] = 2.6
        load_on_cond[2] = -11.47

        cond.SetValue(StructuralMechanicsApplication.POINT_MOMENT, load_on_cond)

        cond.CalculateLocalSystem(lhs,rhs,mp.ProcessInfo)

        self.assertEqual(rhs[0], load_on_cond[0])
        self.assertEqual(rhs[1], load_on_cond[1])
        self.assertEqual(rhs[2], load_on_cond[2])

        # now we apply a load to the node
        nodal_load = KratosMultiphysics.Vector(3)
        nodal_load[0] = -5.5
        nodal_load[1] = 1.2
        nodal_load[2] = 9.3

        node.SetSolutionStepValue(StructuralMechanicsApplication.POINT_MOMENT, nodal_load)

        cond.CalculateLocalSystem(lhs,rhs,mp.ProcessInfo)

        self.assertEqual(rhs[0], load_on_cond[0] + nodal_load[0])
        self.assertEqual(rhs[1], load_on_cond[1] + nodal_load[1])
        self.assertEqual(rhs[2], load_on_cond[2] + nodal_load[2])


    def test_PointLoadCondition2D1N(self):
        pass
        #self._execute_point_load_condition_test(Dimension=2)

    def test_PointLoadCondition3D1N(self):
        pass
        #self._execute_point_load_condition_test(Dimension=3)

    def test_PointMomentCondition3D1N(self):
        pass
        #self._execute_point_moment_condition_test()

    def test_PointMomentCondition3D1N444(self):
        pass
        #self._execute_point_moment_condition_test()


if __name__ == '__main__':
    KratosUnittest.main()
