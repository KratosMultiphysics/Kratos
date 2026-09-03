import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestLoadingConditionsPoint(KratosUnittest.TestCase):

    def _execute_point_load_condition_test(self, current_model, Dimension):
        mp: KratosMultiphysics.ModelPart = current_model.CreateModelPart("solid_part")
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

        cond: KratosMultiphysics.Condition
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

        # Check influencing variables.
        influencing_variables: list[KratosMultiphysics.IAdjoint.DynamicVariable] = cond.GetInfluencingVariables(
            KratosMultiphysics.IAdjoint.ResidualTerm.Mass,
            mp.ProcessInfo)
        self.assertEqual(
            len(influencing_variables),
            0)

        influencing_variables = cond.GetInfluencingVariables(
            KratosMultiphysics.IAdjoint.ResidualTerm.Damping,
            mp.ProcessInfo)
        self.assertEqual(
            len(influencing_variables),
            0)

        influencing_variables = cond.GetInfluencingVariables(
            KratosMultiphysics.IAdjoint.ResidualTerm.Stiffness,
            mp.ProcessInfo)
        self.assertEqual(
            len(influencing_variables),
            0)

        influencing_variables = cond.GetInfluencingVariables(
            KratosMultiphysics.IAdjoint.ResidualTerm.Load,
            mp.ProcessInfo)
        reference_influencing_variables = [KratosMultiphysics.IAdjoint.DynamicVariable(variable, 0) for variable in [
            KratosMultiphysics.StructuralMechanicsApplication.POINT_LOAD_X,
            KratosMultiphysics.StructuralMechanicsApplication.POINT_LOAD_Y,
            KratosMultiphysics.StructuralMechanicsApplication.POINT_LOAD_Z
        ]][:Dimension]
        for reference_variable in reference_influencing_variables:
            self.assertIn(
                reference_variable,
                influencing_variables)

        # Check derivatives (requested in a different order than the condition provided them).
        influencing_variables = influencing_variables[::-1]
        derivatives: KratosMultiphysics.Matrix = cond.ComputeDerivative(
            KratosMultiphysics.IAdjoint.ResidualTerm.Load,
            influencing_variables,
            mp.ProcessInfo,
            0)

        for i_row in range(derivatives.Size1()):
            for i_variable, variable in enumerate(influencing_variables):
                i_component = reference_influencing_variables.index(variable)
                self.assertEqual(
                    1.0 if i_row == i_component else 0.0,
                    derivatives[i_row, i_variable],
                    msg = f"{i_row} {i_component} : {variable.Name()} {i_component}")

    def _execute_point_moment_condition_test(self, current_model):
        mp = current_model.CreateModelPart("solid_part")
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
        current_model = KratosMultiphysics.Model()
        self._execute_point_load_condition_test(current_model, Dimension=2)

    def test_PointLoadCondition3D1N(self):
        current_model = KratosMultiphysics.Model()
        self._execute_point_load_condition_test(current_model, Dimension=3)

    def test_PointMomentCondition3D1N(self):
        current_model = KratosMultiphysics.Model()
        self._execute_point_moment_condition_test(current_model)


if __name__ == '__main__':
    KratosUnittest.main()
