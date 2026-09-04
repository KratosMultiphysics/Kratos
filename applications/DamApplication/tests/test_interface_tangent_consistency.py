import math
import os

import KratosMultiphysics as KM
import KratosMultiphysics.DamApplication
import KratosMultiphysics.PoromechanicsApplication as KratosPoro
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.DamApplication.dam_analysis import DamAnalysis


class TestSmallDisplacementInterfaceTangent(KratosUnittest.TestCase):
    """Directional finite-difference diagnostics for ElasticCohesive2DLaw.

    CalculateLocalSystem returns the positive stiffness matrix and the negative
    internal-force residual.  Consequently this test compares d(RHS)/du to
    ``-LHS * p``.  States that straddle zero normal jump are intentionally
    recorded but not asserted: the elastic cohesive law is discontinuously
    differentiable at that point.
    """

    _eps_scales = (1.0e-4, 1.0e-3, 1.0e-2)

    def _make_element(self, distort_current_coordinates_before_initialize=False):
        model = KM.Model()
        model_part = model.CreateModelPart("interface", 2)
        model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KM.VOLUME_ACCELERATION)
        model_part.AddNodalSolutionStepVariable(KratosPoro.INITIAL_STRESS_TENSOR)
        model_part.AddNodalSolutionStepVariable(KratosPoro.NODAL_JOINT_WIDTH)
        model_part.AddNodalSolutionStepVariable(KratosPoro.NODAL_JOINT_AREA)

        # A unit, zero-thickness vertical interface. Nodes 1--2 and 3--4 are
        # its two faces; the local normal is x and the local tangent is y.
        for node_id, x, y in ((1, 0.0, 0.0), (2, 0.0, 1.0),
                              (3, 0.0, 1.0), (4, 0.0, 0.0)):
            node = model_part.CreateNewNode(node_id, x, y, 0.0)
            node.AddDof(KM.DISPLACEMENT_X)
            node.AddDof(KM.DISPLACEMENT_Y)
            node.SetSolutionStepValue(KM.DISPLACEMENT, KM.Array3())
            node.SetSolutionStepValue(KM.VOLUME_ACCELERATION, KM.Array3())
            initial_stress = KM.Matrix(2, 2)
            for i in range(2):
                for j in range(2):
                    initial_stress[i, j] = 0.0
            node.SetSolutionStepValue(KratosPoro.INITIAL_STRESS_TENSOR, initial_stress)

        properties = model_part.CreateNewProperties(1)
        properties.SetValue(KM.DENSITY, 0.0)
        properties.SetValue(KM.THICKNESS, 1.0)
        # These legacy Poromechanics variables are registered globally rather
        # than exported as Python module attributes.
        properties.SetValue(KM.KratosGlobals.GetVariable("INITIAL_JOINT_WIDTH"), 0.0)
        properties.SetValue(KM.KratosGlobals.GetVariable("NORMAL_STIFFNESS"), 3.0e5)
        properties.SetValue(KM.KratosGlobals.GetVariable("SHEAR_STIFFNESS"), 3.0e5)
        properties.SetValue(KM.KratosGlobals.GetVariable("PENALTY_STIFFNESS"), 1.0e3)
        properties.SetValue(KM.CONSTITUTIVE_LAW, KratosPoro.ElasticCohesive2DLaw())
        element = model_part.CreateNewElement(
            "SmallDisplacementInterfaceElement2D4N", 1, [1, 2, 3, 4], properties)
        model_part.ProcessInfo.SetValue(KratosPoro.IS_CONVERGED, True)
        if distort_current_coordinates_before_initialize:
            model_part.GetNode(3).X += 1.0e-2
            model_part.GetNode(4).X += 1.0e-2
        element.Initialize(model_part.ProcessInfo)
        return model_part, element

    @staticmethod
    def _set_displacement(model_part, direction, scale):
        for i, node in enumerate(model_part.Nodes):
            displacement = KM.Array3()
            displacement[0] = scale * direction[2 * i]
            displacement[1] = scale * direction[2 * i + 1]
            node.SetSolutionStepValue(KM.DISPLACEMENT, displacement)

    @staticmethod
    def _local_normal_jump(element, process_info):
        jumps = element.CalculateOnIntegrationPoints(
            KratosPoro.LOCAL_RELATIVE_DISPLACEMENT_VECTOR, process_info)
        return sum(value[1] for value in jumps) / len(jumps)

    @staticmethod
    def _norm(vector):
        return math.sqrt(sum(value * value for value in vector))

    def _rhs(self, element, process_info):
        lhs = KM.Matrix()
        rhs = KM.Vector()
        element.CalculateLocalSystem(lhs, rhs, process_info)
        return rhs

    @staticmethod
    def _local_system(element, process_info):
        lhs = KM.Matrix()
        rhs = KM.Vector()
        element.CalculateLocalSystem(lhs, rhs, process_info)
        return lhs, rhs

    def _assert_matrix_near(self, first, second, tolerance=1.0e-12):
        self.assertEqual(first.Size1(), second.Size1())
        self.assertEqual(first.Size2(), second.Size2())
        for i in range(first.Size1()):
            for j in range(first.Size2()):
                self.assertAlmostEqual(first[i, j], second[i, j], delta=tolerance)

    def _assert_vector_near(self, first, second, tolerance=1.0e-12):
        self.assertEqual(len(first), len(second))
        for first_i, second_i in zip(first, second):
            self.assertAlmostEqual(first_i, second_i, delta=tolerance)

    def _error(self, model_part, element, direction, normal_jump, epsilon):
        self._set_displacement(model_part, direction, normal_jump + epsilon)
        rhs_plus = self._rhs(element, model_part.ProcessInfo)
        self._set_displacement(model_part, direction, normal_jump - epsilon)
        rhs_minus = self._rhs(element, model_part.ProcessInfo)
        self._set_displacement(model_part, direction, normal_jump)
        lhs = KM.Matrix()
        rhs = KM.Vector()
        element.CalculateLocalSystem(lhs, rhs, model_part.ProcessInfo)

        finite_difference = [(rhs_plus[i] - rhs_minus[i]) / (2.0 * epsilon)
                             for i in range(len(rhs))]
        analytical = [-sum(lhs[i, j] * direction[j] for j in range(len(rhs)))
                      for i in range(len(rhs))]
        difference = [a - b for a, b in zip(finite_difference, analytical)]
        return self._norm(difference) / max(self._norm(finite_difference), self._norm(analytical), 1.0)

    def test_tangent_is_consistent_inside_opening_and_compression_branches(self):
        model_part, element = self._make_element()

        # Apply opposite face displacements. Flip the direction once, if
        # necessary, so positive scale is definitively the opening branch.
        direction = [-1.0, 0.0, -1.0, 0.0, 1.0, 0.0, 1.0, 0.0]
        self._set_displacement(model_part, direction, 1.0)
        if self._local_normal_jump(element, model_part.ProcessInfo) < 0.0:
            direction = [-value for value in direction]

        # The small states and epsilons stay in their branch. The zero state
        # is not asserted because its central perturbation crosses the switch.
        for state in (1.0e-3, 1.0e-5, -1.0e-5, -1.0e-3):
            branch = "opening" if state > 0.0 else "compression"
            for epsilon_scale in self._eps_scales:
                epsilon = max(abs(state), 1.0e-5) * epsilon_scale
                error = self._error(model_part, element, direction, state, epsilon)
                KM.Logger.PrintInfo("InterfaceTangentDiagnostic",
                    "state={} branch={} epsilon={} crossing_zero=false relative_error={}".format(
                        state, branch, epsilon, error))
                self.assertLess(error, 1.0e-6,
                    "Tangent inconsistency at normal jump {} and epsilon {}".format(state, epsilon))

        for epsilon_scale in self._eps_scales:
            epsilon = 1.0e-5 * epsilon_scale
            crossing_error = self._error(model_part, element, direction, 0.0, epsilon)
            KM.Logger.PrintInfo("InterfaceTangentDiagnostic",
                "state=0.0 branch=boundary epsilon={} crossing_zero=true relative_error={}".format(
                    epsilon, crossing_error))
            self.assertGreater(crossing_error, 1.0e-3,
                "The zero-jump diagnostic must continue to cross the constitutive switch")

    def test_displacement_controlled_cycle_exercises_both_branches(self):
        test_directory = os.path.dirname(os.path.realpath(__file__))
        parameter_file = os.path.join(test_directory,
            "joint_elastic_cohesive_2d_normal",
            "joint_elastic_cohesive_2d_normal_parameters.json")
        with open(parameter_file, "r") as parameters_stream:
            parameters = KM.Parameters(parameters_stream.read())

        history = (2.0e-4, 5.0e-4, 2.0e-5, -2.0e-5,
                   -2.0e-4, -2.0e-5, 2.0e-5, 2.0e-4)
        parameters["problem_data"]["end_time"].SetDouble(float(len(history)))
        parameters["problem_data"]["time_step"].SetDouble(1.0)
        parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString(
            os.path.join(test_directory, "joint_elastic_cohesive_2d_normal",
                         "joint_elastic_cohesive_2d_normal"))
        parameters["loads_process_list"][0]["Parameters"]["modulus"].SetString("0.0")

        class DisplacementControlledAnalysis(DamAnalysis):
            boundary_node_ids = (7, 11, 12, 16, 20, 24)

            def Initialize(self):
                super().Initialize()
                self.normal_jumps = []
                for node_id in self.boundary_node_ids:
                    self.main_model_part.GetNode(node_id).Fix(KM.DISPLACEMENT_X)

            def RunMainTemporalLoop(self):
                for step, prescribed_displacement in enumerate(history, 1):
                    self.time += self.delta_time
                    self.main_model_part.CloneTimeStep(self.time)
                    self.main_model_part.ProcessInfo[KM.STEP] = step
                    for node_id in self.boundary_node_ids:
                        self.main_model_part.GetNode(node_id).SetSolutionStepValue(
                            KM.DISPLACEMENT_X, prescribed_displacement)
                    for process in self.list_of_processes:
                        process.ExecuteInitializeSolutionStep()
                    self.solver.Solve()
                    interface = self.main_model_part.GetElement(21)
                    local_jumps = interface.CalculateOnIntegrationPoints(
                        KratosPoro.LOCAL_RELATIVE_DISPLACEMENT_VECTOR,
                        self.main_model_part.ProcessInfo)
                    self.normal_jumps.append(local_jumps[0][1])
                    for process in self.list_of_processes:
                        process.ExecuteFinalizeSolutionStep()

        analysis = DisplacementControlledAnalysis(KM.Model(), parameters)
        analysis.Run()
        branches = [jump > 0.0 for jump in analysis.normal_jumps]
        self.assertTrue(any(branches), "Cycle did not exercise opening")
        self.assertTrue(any(not branch for branch in branches), "Cycle did not exercise compression")
        self.assertTrue(any(branches[i-1] and not branches[i] for i in range(1, len(branches))))
        self.assertTrue(any(not branches[i-1] and branches[i] for i in range(1, len(branches))))

    def test_reference_geometry_is_invariant_while_opening_remains_live(self):
        model_part, element = self._make_element()

        # Besides opening the two faces, distort the moved mid-plane enough
        # that a current-geometry rotation and measure would visibly change.
        opening_displacements = (
            (-1.0e-3, 0.0), (9.0e-2, 2.0e-1),
            (9.2e-2, 2.0e-1), (1.0e-3, 0.0))
        for node, displacement in zip(model_part.Nodes, opening_displacements):
            value = KM.Array3()
            value[0], value[1] = displacement
            node.SetSolutionStepValue(KM.DISPLACEMENT, value)

        trial_jump = element.CalculateOnIntegrationPoints(
            KratosPoro.LOCAL_RELATIVE_DISPLACEMENT_VECTOR,
            model_part.ProcessInfo)
        if sum(value[1] for value in trial_jump) < 0.0:
            for node in model_part.Nodes:
                displacement = node.GetSolutionStepValue(KM.DISPLACEMENT)
                displacement[0] *= -1.0
                displacement[1] *= -1.0
                node.SetSolutionStepValue(KM.DISPLACEMENT, displacement)

        reference_lhs, reference_rhs = self._local_system(
            element, model_part.ProcessInfo)
        reference_jump = element.CalculateOnIntegrationPoints(
            KratosPoro.LOCAL_RELATIVE_DISPLACEMENT_VECTOR,
            model_part.ProcessInfo)

        for node in model_part.Nodes:
            displacement = node.GetSolutionStepValue(KM.DISPLACEMENT)
            node.X = node.X0 + displacement[0]
            node.Y = node.Y0 + displacement[1]

        moved_lhs, moved_rhs = self._local_system(element, model_part.ProcessInfo)
        moved_jump = element.CalculateOnIntegrationPoints(
            KratosPoro.LOCAL_RELATIVE_DISPLACEMENT_VECTOR,
            model_part.ProcessInfo)
        opening_stress = element.CalculateOnIntegrationPoints(
            KratosPoro.LOCAL_STRESS_VECTOR, model_part.ProcessInfo)

        self._assert_matrix_near(reference_lhs, moved_lhs)
        self._assert_vector_near(reference_rhs, moved_rhs)
        for reference_value, moved_value in zip(reference_jump, moved_jump):
            self._assert_vector_near(reference_value, moved_value)

        opening_normal_jump = sum(value[1] for value in moved_jump) / len(moved_jump)
        self.assertGreater(opening_normal_jump, 0.0)
        opening_joint_width = max(opening_normal_jump, 0.0)

        for node in model_part.Nodes:
            displacement = node.GetSolutionStepValue(KM.DISPLACEMENT)
            displacement[0] *= -1.0
            displacement[1] *= -1.0
            node.SetSolutionStepValue(KM.DISPLACEMENT, displacement)
            node.X = node.X0 + displacement[0]
            node.Y = node.Y0 + displacement[1]

        compressed_jump = element.CalculateOnIntegrationPoints(
            KratosPoro.LOCAL_RELATIVE_DISPLACEMENT_VECTOR,
            model_part.ProcessInfo)
        compressed_stress = element.CalculateOnIntegrationPoints(
            KratosPoro.LOCAL_STRESS_VECTOR, model_part.ProcessInfo)
        compressed_normal_jump = sum(value[1] for value in compressed_jump) / len(compressed_jump)
        compressed_joint_width = max(compressed_normal_jump, 0.0)
        opening_normal_stress = sum(value[1] for value in opening_stress) / len(opening_stress)
        compressed_normal_stress = sum(value[1] for value in compressed_stress) / len(compressed_stress)

        self.assertLess(compressed_normal_jump, 0.0)
        self.assertGreater(opening_joint_width, compressed_joint_width)
        self.assertEqual(compressed_joint_width, 0.0)
        self.assertGreater(opening_normal_stress, 0.0)
        self.assertLess(compressed_normal_stress, 0.0)
        self.assertGreater(
            abs(compressed_normal_stress), 100.0 * abs(opening_normal_stress))

    def test_initial_gap_uses_reference_geometry(self):
        reference_model_part, reference_element = self._make_element()
        moved_model_part, moved_element = self._make_element(True)

        direction = [-1.0, 0.0, -1.0, 0.0, 1.0, 0.0, 1.0, 0.0]
        self._set_displacement(reference_model_part, direction, 1.0e-4)
        self._set_displacement(moved_model_part, direction, 1.0e-4)

        reference_lhs, reference_rhs = self._local_system(
            reference_element, reference_model_part.ProcessInfo)
        moved_lhs, moved_rhs = self._local_system(
            moved_element, moved_model_part.ProcessInfo)
        self._assert_matrix_near(reference_lhs, moved_lhs)
        self._assert_vector_near(reference_rhs, moved_rhs)

    def test_nonzero_initial_stress_uses_reference_frame(self):
        model_part, element = self._make_element()
        initial_stress = KM.Matrix(2, 2)
        initial_stress[0, 0] = 4.0
        initial_stress[0, 1] = 1.0
        initial_stress[1, 0] = 1.0
        initial_stress[1, 1] = 2.0
        for node in model_part.Nodes:
            node.SetSolutionStepValue(KratosPoro.INITIAL_STRESS_TENSOR, initial_stress)

        reference_lhs, reference_rhs = self._local_system(
            element, model_part.ProcessInfo)
        reference_stress = element.CalculateOnIntegrationPoints(
            KratosPoro.LOCAL_STRESS_VECTOR, model_part.ProcessInfo)

        model_part.GetNode(2).X += 2.0e-1
        model_part.GetNode(3).X += 2.0e-1
        moved_lhs, moved_rhs = self._local_system(element, model_part.ProcessInfo)
        moved_stress = element.CalculateOnIntegrationPoints(
            KratosPoro.LOCAL_STRESS_VECTOR, model_part.ProcessInfo)

        self._assert_matrix_near(reference_lhs, moved_lhs)
        self._assert_vector_near(reference_rhs, moved_rhs)
        for reference_value, moved_value in zip(reference_stress, moved_stress):
            self._assert_vector_near(reference_value, moved_value)
