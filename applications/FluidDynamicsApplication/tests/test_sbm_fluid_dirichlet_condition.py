import math

import KratosMultiphysics as KM
import KratosMultiphysics.FluidDynamicsApplication as CFD
import KratosMultiphysics.KratosUnittest as KratosUnittest


class SbmFluidDirichletConditionTest(KratosUnittest.TestCase):

    def test_face_condition_created_by_surrogate_boundary_modeler(self):
        model = KM.Model()
        skin = model.CreateModelPart("Skin")
        skin_properties = skin.CreateNewProperties(1)
        skin.CreateNewNode(1, 0.75, 0.75, 0.0)
        skin.CreateNewNode(2, 1.25, 0.75, 0.0)
        skin.CreateNewNode(3, 1.25, 1.25, 0.0)
        skin.CreateNewNode(4, 0.75, 1.25, 0.0)
        skin.CreateNewCondition("LineCondition2D2N", 1, [1, 2], skin_properties)
        skin.CreateNewCondition("LineCondition2D2N", 2, [2, 3], skin_properties)
        skin.CreateNewCondition("LineCondition2D2N", 3, [3, 4], skin_properties)
        skin.CreateNewCondition("LineCondition2D2N", 4, [4, 1], skin_properties)

        background = model.CreateModelPart("Background")
        background.AddNodalSolutionStepVariable(KM.VELOCITY)
        background.AddNodalSolutionStepVariable(KM.ACCELERATION)
        background.AddNodalSolutionStepVariable(KM.PRESSURE)
        properties = background.CreateNewProperties(1)
        properties.SetValue(KM.DYNAMIC_VISCOSITY, 0.8)
        properties.SetValue(KM.CONSTITUTIVE_LAW, CFD.Newtonian2DLaw())
        properties.SetValue(KM.PENALTY_COEFFICIENT, 10.0)

        settings = KM.Parameters(r'''{
            "input_model_part_name" : "Skin",
            "output_model_part_name" : "Background",
            "lower_point" : [0.0, 0.0, 0.0],
            "upper_point" : [2.0, 2.0, 0.0],
            "number_of_elements" : [8, 8, 1],
            "physical_domain" : "outside",
            "lambda" : 1.0,
            "condition_name" : "SbmFluidDirichletCondition2D4N"
        }''')
        modeler = KM.CreateModeler("SurrogateBoundaryModeler", model, settings)
        modeler.SetupModelPart()

        surrogate_conditions = list(
            background.GetSubModelPart("SurrogateBoundary").Conditions
        )
        face_keys = set()
        for surrogate_condition in surrogate_conditions:
            face = surrogate_condition.GetValue(
                KM.SURROGATE_BOUNDARY_FACE_COORDINATES
            )
            face_keys.add(
                tuple(
                    sorted(
                        (
                            (round(face[0, 0], 12), round(face[0, 1], 12)),
                            (round(face[1, 0], 12), round(face[1, 1], 12)),
                        )
                    )
                )
            )
        self.assertEqual(len(face_keys), len(surrogate_conditions))

        for node in background.Nodes:
            node.AddDof(KM.VELOCITY_X)
            node.AddDof(KM.VELOCITY_Y)
            node.AddDof(KM.PRESSURE)
            velocity = node.GetSolutionStepValue(KM.VELOCITY)
            velocity[0] = 0.05 + 0.02 * node.X
            velocity[1] = -0.03 + 0.01 * node.Y
            node.SetSolutionStepValue(KM.PRESSURE, 1.0 + 0.1 * node.X)

        condition = surrogate_conditions[0]
        prescribed_velocities = KM.Matrix(2, 3)
        prescribed_velocities[0, 0] = 0.2
        prescribed_velocities[0, 1] = -0.1
        prescribed_velocities[1, 0] = 0.25
        prescribed_velocities[1, 1] = -0.08
        condition.SetValue(CFD.SBM_BOUNDARY_VELOCITIES, prescribed_velocities)
        condition.Initialize(background.ProcessInfo)

        lhs_penalty = KM.Matrix()
        rhs_penalty = KM.Vector()
        condition.CalculateLocalVelocityContribution(
            lhs_penalty, rhs_penalty, background.ProcessInfo
        )
        self.assertEqual(lhs_penalty.Size1(), 12)
        self.assertEqual(lhs_penalty.Size2(), 12)
        self.assertEqual(rhs_penalty.Size(), 12)
        self.assertTrue(
            all(
                math.isfinite(lhs_penalty[i, j])
                for i in range(lhs_penalty.Size1())
                for j in range(lhs_penalty.Size2())
            )
        )
        self.assertTrue(all(math.isfinite(value) for value in rhs_penalty))

        properties.SetValue(KM.PENALTY_COEFFICIENT, 0.0)
        lhs_penalty_free = KM.Matrix()
        rhs_penalty_free = KM.Vector()
        condition.CalculateLocalVelocityContribution(
            lhs_penalty_free, rhs_penalty_free, background.ProcessInfo
        )
        maximum_lhs_difference = max(
            abs(lhs_penalty[i, j] - lhs_penalty_free[i, j])
            for i in range(lhs_penalty.Size1())
            for j in range(lhs_penalty.Size2())
        )
        maximum_rhs_difference = max(
            abs(rhs_penalty[i] - rhs_penalty_free[i])
            for i in range(rhs_penalty.Size())
        )
        self.assertGreater(maximum_lhs_difference, 1.0e-8)
        self.assertGreater(maximum_rhs_difference, 1.0e-8)


if __name__ == "__main__":
    KratosUnittest.main()
