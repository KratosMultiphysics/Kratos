# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.StructuralMechanicsApplication.step_controller import DefaultStepController, GeometricStepController

class TestStepControllers(KratosUnittest.TestCase):
    def test_DefaultStepController(self):
        params = KratosMultiphysics.Parameters("""{}""")
        step_controller = DefaultStepController(params)

        step_controller.Initialize(10, 15)
        self.assertTrue(step_controller.IsCompleted(100, False))
        self.assertTrue(step_controller.IsCompleted(100, True))
        self.assertTrue(step_controller.IsCompleted(12, False))
        self.assertTrue(step_controller.IsCompleted(12, True))
        self.assertTrue(step_controller.IsCompleted(5, False))
        self.assertTrue(step_controller.IsCompleted(5, True))

        self.assertEqual(step_controller.GetNextStep(14, False), 15)
        self.assertEqual(step_controller.GetNextStep(14, True), 15)

    def test_GeometricStepController1(self):
        params = KratosMultiphysics.Parameters("""{
            "type"                                       : "geometric_step_controller",
            "divergence_factor"                          : 0.25,
            "convergence_factor"                         : 1.5,
            "delta_t_init"                               : 1.1,
            "delta_t_min"                                : 0.12,
            "delta_t_max"                                : 3.1,
            "max_number_of_sub_steps"                    : 25,
            "number_of_successful_attempts_for_increment": 2,
            "number_of_failed_attempts_for_termination"  : 5
        }""")

        step_controller = GeometricStepController(params)
        step_controller.Initialize(0.2, 10.0)

        self.assertFalse(step_controller.IsCompleted(0.1, True))
        self.assertFalse(step_controller.IsCompleted(0.2, True))
        self.assertFalse(step_controller.IsCompleted(5, True))
        self.assertFalse(step_controller.IsCompleted(10, False))
        self.assertTrue(step_controller.IsCompleted(10, True))

        self.assertAlmostEqual(step_controller.GetNextStep(0.3, True), 0.3 + 1.1)
        self.assertAlmostEqual(step_controller.GetNextStep(0.3, True), 0.3 + 1.1)
        self.assertAlmostEqual(step_controller.GetNextStep(0.3, True), 0.3 + 1.1 * 1.5)
        self.assertAlmostEqual(step_controller.GetNextStep(0.3, True), 0.3 + 1.1 * 1.5 * 1.5)
        self.assertAlmostEqual(step_controller.GetNextStep(0.3, True), 0.3 + 3.1)

        self.assertAlmostEqual(step_controller.GetNextStep(0.3, False), 0.3 + 3.1 * 0.25)
        self.assertAlmostEqual(step_controller.GetNextStep(0.3, False), 0.3 + 3.1 * 0.25 * 0.25)

        with self.assertRaises(RuntimeError):
            self.assertAlmostEqual(step_controller.GetNextStep(0.3, False), 0.3 + 3.1 * 0.25 * 0.25 * 0.25)

        self.assertAlmostEqual(step_controller.GetNextStep(0.3, True), 0.3 + 3.1 * 0.25 * 0.25 * 0.25)
        self.assertAlmostEqual(step_controller.GetNextStep(0.3, True), 0.3 + 3.1 * 0.25 * 0.25 * 0.25)
        self.assertAlmostEqual(step_controller.GetNextStep(0.3, True), 0.3 + 3.1 * 0.25 * 0.25 * 0.25 * 1.5)
        self.assertAlmostEqual(step_controller.GetNextStep(0.3, True), 0.3 + 3.1 * 0.25 * 0.25 * 0.25 * 1.5 * 1.5)
        self.assertAlmostEqual(step_controller.GetNextStep(0.3, True), 0.3 + 3.1 * 0.25 * 0.25 * 0.25 * 1.5 * 1.5 * 1.5)
        self.assertAlmostEqual(step_controller.GetNextStep(0.3, True), 0.3 + 3.1 * 0.25 * 0.25 * 0.25 * 1.5 * 1.5 * 1.5 * 1.5)
        self.assertAlmostEqual(step_controller.GetNextStep(0.3, True), 0.3 + 3.1 * 0.25 * 0.25 * 0.25 * 1.5 * 1.5 * 1.5 * 1.5 * 1.5)
        self.assertAlmostEqual(step_controller.GetNextStep(0.3, True), 0.3 + 3.1 * 0.25 * 0.25 * 0.25 * 1.5 * 1.5 * 1.5 * 1.5 * 1.5 * 1.5)
        self.assertAlmostEqual(step_controller.GetNextStep(0.3, True), 0.3 + 3.1 * 0.25 * 0.25 * 0.25 * 1.5 * 1.5 * 1.5 * 1.5 * 1.5 * 1.5 * 1.5)
        self.assertAlmostEqual(step_controller.GetNextStep(0.3, True), 0.3 + 3.1 * 0.25 * 0.25 * 0.25 * 1.5 * 1.5 * 1.5 * 1.5 * 1.5 * 1.5 * 1.5 * 1.5)
        self.assertAlmostEqual(step_controller.GetNextStep(0.3, True), 0.3 + 3.1 * 0.25 * 0.25 * 0.25 * 1.5 * 1.5 * 1.5 * 1.5 * 1.5 * 1.5 * 1.5 * 1.5 * 1.5)
        self.assertAlmostEqual(step_controller.GetNextStep(0.3, True), 0.3 + 3.1 * 0.25 * 0.25 * 0.25 * 1.5 * 1.5 * 1.5 * 1.5 * 1.5 * 1.5 * 1.5 * 1.5 * 1.5 * 1.5)
        self.assertAlmostEqual(step_controller.GetNextStep(0.3, True), 0.3 + 3.1)

    def test_GeometricStepController2(self):
        params = KratosMultiphysics.Parameters("""{
            "type"                                       : "geometric_step_controller",
            "divergence_factor"                          : 0.25,
            "convergence_factor"                         : 1.5,
            "delta_t_init"                               : 1.1,
            "delta_t_min"                                : 0.001,
            "delta_t_max"                                : 3.1,
            "max_number_of_sub_steps"                    : 25,
            "number_of_successful_attempts_for_increment": 2,
            "number_of_failed_attempts_for_termination"  : 2
        }""")

        step_controller = GeometricStepController(params)
        step_controller.Initialize(0.2, 10.0)

        step_controller.GetNextStep(0.3, False)
        step_controller.GetNextStep(0.3, False)
        with self.assertRaises(RuntimeError):
            step_controller.GetNextStep(0.3, False)

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()