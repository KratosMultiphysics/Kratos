import os
import sys

import KratosMultiphysics as KM
import KratosMultiphysics.DamApplication as KratosDam
import KratosMultiphysics.PoromechanicsApplication as KratosPoro
import KratosMultiphysics.KratosUnittest as KratosUnittest

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
if TESTS_DIR not in sys.path:
    sys.path.insert(0, TESTS_DIR)

import test_interface_tangent_consistency as interface_test


class TestDamInitialStressHandling(KratosUnittest.TestCase):
    """Regression tests for INITIAL_STRESS_TENSOR storage in the Dam schemes."""

    @staticmethod
    def _make_model_part():
        model = KM.Model()
        model_part = model.CreateModelPart("initial_stress_test")
        model_part.AddNodalSolutionStepVariable(KratosPoro.INITIAL_STRESS_TENSOR)
        model_part.ProcessInfo.SetValue(KM.DELTA_TIME, 0.1)
        model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        return model_part

    @staticmethod
    def _schemes():
        return [
            KratosDam.IncrementalUpdateStaticSmoothingScheme(),
            KratosDam.IncrementalUpdateStaticDampedSmoothingScheme(0.0, 0.0),
            KratosDam.BossakDisplacementSmoothingScheme(0.0, 0.0, 0.0),
            KratosDam.DamUPScheme(0.25, 0.5, 0.0, 0.0),
            KratosDam.DamPScheme(0.25, 0.5),
        ]

    def _assert_tensor(self, tensor, expected):
        self.assertEqual(tensor.Size1(), 3)
        self.assertEqual(tensor.Size2(), 3)
        for i in range(3):
            for j in range(3):
                self.assertAlmostEqual(tensor[i, j], expected[i][j], places=14)

    def test_schemes_preserve_preloaded_3x3_initial_stress(self):
        expected = [[1.0e6, 2.0e6, 0.0], [2.0e6, 3.0e6, 0.0], [0.0, 0.0, 4.0e6]]
        for scheme in self._schemes():
            with self.subTest(scheme=type(scheme).__name__):
                model_part = self._make_model_part()
                for node in model_part.Nodes:
                    node.SetSolutionStepValue(KratosPoro.INITIAL_STRESS_TENSOR,
                                              KM.Matrix(expected))
                scheme.Initialize(model_part)
                for node in model_part.Nodes:
                    self._assert_tensor(
                        node.GetSolutionStepValue(KratosPoro.INITIAL_STRESS_TENSOR),
                        expected)

    def test_schemes_initialize_empty_tensor_to_zero_3x3(self):
        for scheme in self._schemes():
            with self.subTest(scheme=type(scheme).__name__):
                model_part = self._make_model_part()
                scheme.Initialize(model_part)
                for node in model_part.Nodes:
                    self._assert_tensor(
                        node.GetSolutionStepValue(KratosPoro.INITIAL_STRESS_TENSOR),
                        [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])

    def test_schemes_convert_2x2_initial_stress_to_3x3(self):
        expected = [[1.0e6, 2.0e6, 0.0], [2.0e6, 3.0e6, 0.0], [0.0, 0.0, 0.0]]
        for scheme in self._schemes():
            with self.subTest(scheme=type(scheme).__name__):
                model_part = self._make_model_part()
                initial_stress = KM.Matrix(2, 2)
                initial_stress[0, 0] = 1.0e6
                initial_stress[0, 1] = 2.0e6
                initial_stress[1, 0] = 2.0e6
                initial_stress[1, 1] = 3.0e6
                for node in model_part.Nodes:
                    node.SetSolutionStepValue(KratosPoro.INITIAL_STRESS_TENSOR,
                                              initial_stress)
                scheme.Initialize(model_part)
                scheme.Initialize(model_part)  # idempotence
                for node in model_part.Nodes:
                    self._assert_tensor(
                        node.GetSolutionStepValue(KratosPoro.INITIAL_STRESS_TENSOR),
                        expected)

    def test_schemes_reject_unsupported_tensor_sizes(self):
        for unsupported_size in ((1, 1), (4, 4)):
            for scheme in self._schemes():
                with self.subTest(scheme=type(scheme).__name__,
                                  size=unsupported_size):
                    model_part = self._make_model_part()
                    matrix = KM.Matrix(unsupported_size[0], unsupported_size[1])
                    for node in model_part.Nodes:
                        node.SetSolutionStepValue(KratosPoro.INITIAL_STRESS_TENSOR,
                                                  matrix)
                    with self.assertRaises(RuntimeError):
                        scheme.Initialize(model_part)

    def test_interface_initial_stress_survives_scheme_initialization(self):
        test_case = interface_test.TestSmallDisplacementInterfaceTangent(
            "test_initial_gap_uses_reference_geometry")

        initial_stress = KM.Matrix(3, 3)
        for i in range(3):
            for j in range(3):
                initial_stress[i, j] = 1.0e6 if i < 2 and j < 2 else 0.0

        model_part, element = test_case._make_element()
        for node in model_part.Nodes:
            node.SetSolutionStepValue(KratosPoro.INITIAL_STRESS_TENSOR,
                                      initial_stress)
        _, rhs_before = test_case._local_system(element, model_part.ProcessInfo)

        KratosDam.IncrementalUpdateStaticSmoothingScheme().Initialize(model_part)

        for node in model_part.Nodes:
            self._assert_tensor(
                node.GetSolutionStepValue(KratosPoro.INITIAL_STRESS_TENSOR),
                [[1.0e6, 1.0e6, 0.0], [1.0e6, 1.0e6, 0.0], [0.0, 0.0, 0.0]])

        _, rhs_after = test_case._local_system(element, model_part.ProcessInfo)
        test_case._assert_vector_near(rhs_before, rhs_after)

        zero_model_part, zero_element = test_case._make_element()
        _, rhs_zero = test_case._local_system(
            zero_element, zero_model_part.ProcessInfo)
        self.assertGreater(
            max(abs(a - b) for a, b in zip(rhs_before, rhs_zero)), 1.0)


if __name__ == "__main__":
    KratosUnittest.main()
