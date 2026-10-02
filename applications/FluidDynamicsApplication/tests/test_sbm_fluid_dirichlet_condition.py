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

    def test_gap_modeler_creates_exact_boundary_and_interface_entities(self):
        model = KM.Model()
        skin = model.CreateModelPart("Skin")
        skin_properties = skin.CreateNewProperties(1)
        coordinates = (
            (0.73, 0.67),
            (1.31, 0.67),
            (1.31, 1.37),
            (0.73, 1.37),
        )
        for node_id, (x, y) in enumerate(coordinates, start=1):
            skin.CreateNewNode(node_id, x, y, 0.0)
        for condition_id in range(4):
            skin.CreateNewCondition(
                "LineCondition2D2N",
                condition_id + 1,
                [condition_id + 1, (condition_id + 1) % 4 + 1],
                skin_properties,
            )

        background = model.CreateModelPart("Background")
        for variable in (
            KM.VELOCITY,
            KM.ACCELERATION,
            KM.PRESSURE,
            KM.BODY_FORCE,
        ):
            background.AddNodalSolutionStepVariable(variable)
        properties = background.CreateNewProperties(1)
        properties.SetValue(KM.DENSITY, 1.0)
        properties.SetValue(KM.DYNAMIC_VISCOSITY, 0.8)
        properties.SetValue(KM.CONSTITUTIVE_LAW, CFD.Newtonian2DLaw())
        properties.SetValue(KM.PENALTY_COEFFICIENT, 0.0)

        settings = KM.Parameters(r'''{
            "input_model_part_name" : "Skin",
            "output_model_part_name" : "Background",
            "lower_point" : [0.0, 0.0, 0.0],
            "upper_point" : [2.0, 2.0, 0.0],
            "number_of_elements" : [12, 12, 1],
            "physical_domain" : "outside",
            "lambda" : 1.0,
            "element_name" : "SymbolicStokes2D4N",
            "condition_name" : "SbmFluidDirichletCondition2D4N",
            "gap_element_name" : "SymbolicStokes2D4N",
            "gap_interface_condition_name" : "SbmFluidGapInterfaceCondition2D"
        }''')
        modeler = KM.CreateModeler(
            "GapSurrogateBoundaryModeler", model, settings
        )
        modeler.SetupModelPart()
        background.SetBufferSize(3)
        background.ProcessInfo[KM.DELTA_TIME] = 1.0
        background.ProcessInfo[KM.DYNAMIC_TAU] = 0.0
        bdf_coefficients = KM.Vector(3)
        bdf_coefficients[0] = 0.0
        bdf_coefficients[1] = 0.0
        bdf_coefficients[2] = 0.0
        background.ProcessInfo[KM.BDF_COEFFICIENTS] = bdf_coefficients

        surrogate = background.GetSubModelPart("SurrogateBoundary")
        gap_elements = background.GetSubModelPart("GapElements")
        gap_interfaces = background.GetSubModelPart("GapInterfaces")
        self.assertGreater(surrogate.NumberOfConditions(), 0)
        self.assertEqual(
            gap_elements.NumberOfElements(), surrogate.NumberOfConditions()
        )
        self.assertGreater(gap_interfaces.NumberOfConditions(), 0)

        for node in background.Nodes:
            node.AddDof(KM.VELOCITY_X)
            node.AddDof(KM.VELOCITY_Y)
            node.AddDof(KM.PRESSURE)
            node.SetSolutionStepValue(
                KM.VELOCITY,
                [0.1 + 0.03 * node.X, -0.02 + 0.04 * node.Y, 0.0],
            )
            node.SetSolutionStepValue(KM.PRESSURE, 0.2 * node.X - 0.1 * node.Y)

        exact_condition = next(iter(surrogate.Conditions))
        boundary_quadrature = exact_condition.GetValue(
            KM.SURROGATE_BOUNDARY_PROJECTION
        )
        self.assertEqual(boundary_quadrature.Size1(), 2)
        self.assertEqual(boundary_quadrature.Size2(), 7)
        exact_condition.SetValue(CFD.SBM_BOUNDARY_VELOCITIES, KM.Matrix(2, 3))
        exact_condition.Initialize(background.ProcessInfo)
        boundary_lhs = KM.Matrix()
        boundary_rhs = KM.Vector()
        exact_condition.CalculateLocalVelocityContribution(
            boundary_lhs, boundary_rhs, background.ProcessInfo
        )
        self.assertEqual(boundary_lhs.Size1(), 12)
        self.assertTrue(all(math.isfinite(value) for value in boundary_rhs))

        # The gap entity is an unmodified registered Stokes element whose
        # QuadraturePointGeometry stores the extended parent shape functions.
        # Exercising its local system checks that no gap nodes or DOFs are
        # required and that the custom 2x2 quadrature reaches the fluid kernel.
        gap_element = next(iter(gap_elements.Elements))
        gap_element.Initialize(background.ProcessInfo)
        gap_lhs = KM.Matrix()
        gap_rhs = KM.Vector()
        gap_element.CalculateLocalSystem(
            gap_lhs, gap_rhs, background.ProcessInfo
        )
        self.assertEqual(gap_lhs.Size1(), 12)
        self.assertEqual(gap_lhs.Size2(), 12)
        self.assertEqual(gap_rhs.Size(), 12)
        self.assertTrue(
            all(
                math.isfinite(gap_lhs[i, j])
                for i in range(gap_lhs.Size1())
                for j in range(gap_lhs.Size2())
            )
        )
        self.assertTrue(all(math.isfinite(value) for value in gap_rhs))

        interface_condition = next(iter(gap_interfaces.Conditions))
        interface_condition.Initialize(background.ProcessInfo)
        interface_lhs = KM.Matrix()
        interface_rhs = KM.Vector()
        interface_condition.CalculateLocalVelocityContribution(
            interface_lhs, interface_rhs, background.ProcessInfo
        )
        expected_size = 3 * interface_condition.GetGeometry().PointsNumber()
        self.assertEqual(interface_lhs.Size1(), expected_size)
        self.assertEqual(interface_lhs.Size2(), expected_size)
        self.assertEqual(interface_rhs.Size(), expected_size)
        self.assertTrue(
            all(
                math.isfinite(interface_lhs[i, j])
                for i in range(interface_lhs.Size1())
                for j in range(interface_lhs.Size2())
            )
        )
        self.assertTrue(all(math.isfinite(value) for value in interface_rhs))

    def test_gap_modeler_accepts_triangular_patch_at_skin_corner(self):
        model = KM.Model()
        skin = model.CreateModelPart("Skin")
        skin_properties = skin.CreateNewProperties(1)
        coordinates = (
            (0.50, 0.25),
            (0.75, 0.50),
            (0.50, 0.75),
            (0.25, 0.50),
        )
        for node_id, (x, y) in enumerate(coordinates, start=1):
            skin.CreateNewNode(node_id, x, y, 0.0)
        for condition_id in range(4):
            skin.CreateNewCondition(
                "LineCondition2D2N",
                condition_id + 1,
                [condition_id + 1, (condition_id + 1) % 4 + 1],
                skin_properties,
            )

        background = model.CreateModelPart("Background")
        for variable in (
            KM.VELOCITY,
            KM.ACCELERATION,
            KM.PRESSURE,
            KM.BODY_FORCE,
        ):
            background.AddNodalSolutionStepVariable(variable)
        properties = background.CreateNewProperties(1)
        properties.SetValue(KM.DENSITY, 1.0)
        properties.SetValue(KM.DYNAMIC_VISCOSITY, 0.8)
        properties.SetValue(KM.CONSTITUTIVE_LAW, CFD.Newtonian2DLaw())
        properties.SetValue(KM.PENALTY_COEFFICIENT, 0.0)

        settings = KM.Parameters(r'''{
            "input_model_part_name" : "Skin",
            "output_model_part_name" : "Background",
            "lower_point" : [0.0, 0.0, 0.0],
            "upper_point" : [1.0, 1.0, 0.0],
            "number_of_elements" : [16, 16, 1],
            "physical_domain" : "outside",
            "lambda" : 1.0,
            "element_name" : "SymbolicStokes2D4N",
            "condition_name" : "SbmFluidDirichletCondition2D4N",
            "gap_element_name" : "SymbolicStokes2D4N",
            "gap_interface_condition_name" : "SbmFluidGapInterfaceCondition2D"
        }''')
        modeler = KM.CreateModeler(
            "GapSurrogateBoundaryModeler", model, settings
        )
        modeler.SetupModelPart()

        background.SetBufferSize(3)
        background.ProcessInfo[KM.DELTA_TIME] = 1.0
        background.ProcessInfo[KM.DYNAMIC_TAU] = 0.0
        bdf_coefficients = KM.Vector(3)
        background.ProcessInfo[KM.BDF_COEFFICIENTS] = bdf_coefficients
        for node in background.Nodes:
            node.AddDof(KM.VELOCITY_X)
            node.AddDof(KM.VELOCITY_Y)
            node.AddDof(KM.PRESSURE)

        surrogate = background.GetSubModelPart("SurrogateBoundary")
        gap_elements = background.GetSubModelPart("GapElements")
        self.assertEqual(
            gap_elements.NumberOfElements(), surrogate.NumberOfConditions()
        )

        collapsed_conditions = []
        for condition in surrogate.Conditions:
            boundary_quadrature = condition.GetValue(
                KM.SURROGATE_BOUNDARY_PROJECTION
            )
            if all(
                abs(boundary_quadrature[g, 6]) < 1.0e-15
                for g in range(boundary_quadrature.Size1())
            ):
                collapsed_conditions.append(condition)
        self.assertGreater(len(collapsed_conditions), 0)

        for condition in collapsed_conditions:
            condition.SetValue(
                CFD.SBM_BOUNDARY_VELOCITIES, KM.Matrix(2, 3)
            )
            condition.Initialize(background.ProcessInfo)
            lhs = KM.Matrix()
            rhs = KM.Vector()
            condition.CalculateLocalVelocityContribution(
                lhs, rhs, background.ProcessInfo
            )
            self.assertEqual(lhs.Size1(), 12)
            self.assertEqual(lhs.Size2(), 12)
            self.assertEqual(rhs.Size(), 12)
            self.assertTrue(
                all(
                    abs(lhs[i, j]) < 1.0e-15
                    for i in range(lhs.Size1())
                    for j in range(lhs.Size2())
                )
            )
            self.assertTrue(all(abs(value) < 1.0e-15 for value in rhs))

    def test_gap_neumann_condition_integrates_exact_boundary_traction(self):
        model = KM.Model()
        skin = model.CreateModelPart("Skin")
        skin_properties = skin.CreateNewProperties(1)
        coordinates = (
            (0.73, 0.67),
            (1.31, 0.67),
            (1.31, 1.37),
            (0.73, 1.37),
        )
        for node_id, (x, y) in enumerate(coordinates, start=1):
            skin.CreateNewNode(node_id, x, y, 0.0)
        for condition_id in range(4):
            skin.CreateNewCondition(
                "LineCondition2D2N",
                condition_id + 1,
                [condition_id + 1, (condition_id + 1) % 4 + 1],
                skin_properties,
            )

        background = model.CreateModelPart("Background")
        for variable in (KM.VELOCITY, KM.ACCELERATION, KM.PRESSURE):
            background.AddNodalSolutionStepVariable(variable)
        background.CreateNewProperties(1)

        settings = KM.Parameters(r'''{
            "input_model_part_name" : "Skin",
            "output_model_part_name" : "Background",
            "lower_point" : [0.0, 0.0, 0.0],
            "upper_point" : [2.0, 2.0, 0.0],
            "number_of_elements" : [12, 12, 1],
            "physical_domain" : "outside",
            "lambda" : 1.0,
            "element_name" : "SymbolicStokes2D4N",
            "condition_name" : "SbmFluidNeumannCondition2D4N",
            "gap_element_name" : "SymbolicStokes2D4N",
            "gap_interface_condition_name" : "SbmFluidGapInterfaceCondition2D"
        }''')
        modeler = KM.CreateModeler(
            "GapSurrogateBoundaryModeler", model, settings
        )
        modeler.SetupModelPart()

        background.ProcessInfo[KM.BDF_COEFFICIENTS] = KM.Vector(3)
        for node in background.Nodes:
            node.AddDof(KM.VELOCITY_X)
            node.AddDof(KM.VELOCITY_Y)
            node.AddDof(KM.PRESSURE)

        condition = next(
            iter(
                background.GetSubModelPart(
                    "SurrogateBoundary"
                ).Conditions
            )
        )
        quadrature = condition.GetValue(
            KM.SURROGATE_BOUNDARY_PROJECTION
        )
        prescribed_stress = KM.Matrix(2, 4)
        for point_index in range(2):
            prescribed_stress[point_index, 0] = 1.0
            prescribed_stress[point_index, 1] = 0.0
            prescribed_stress[point_index, 2] = 0.0
            prescribed_stress[point_index, 3] = 1.0
        condition.SetValue(KM.CAUCHY_STRESS_TENSOR, prescribed_stress)
        self.assertEqual(condition.Check(background.ProcessInfo), 0)

        lhs = KM.Matrix()
        rhs = KM.Vector()
        condition.CalculateLocalSystem(lhs, rhs, background.ProcessInfo)
        self.assertEqual(lhs.Size1(), 12)
        self.assertEqual(lhs.Size2(), 12)
        self.assertEqual(rhs.Size(), 12)
        self.assertTrue(
            all(
                abs(lhs[i, j]) < 1.0e-15
                for i in range(lhs.Size1())
                for j in range(lhs.Size2())
            )
        )
        self.assertTrue(
            all(
                abs(rhs[3 * node_index + 2]) < 1.0e-15
                for node_index in range(4)
            )
        )

        expected_x = sum(
            quadrature[g, 6] * quadrature[g, 3] for g in range(2)
        )
        expected_y = sum(
            quadrature[g, 6] * quadrature[g, 4] for g in range(2)
        )
        integrated_x = sum(rhs[3 * node_index] for node_index in range(4))
        integrated_y = sum(rhs[3 * node_index + 1] for node_index in range(4))
        self.assertAlmostEqual(integrated_x, expected_x, places=12)
        self.assertAlmostEqual(integrated_y, expected_y, places=12)


if __name__ == "__main__":
    KratosUnittest.main()
