
import KratosMultiphysics as Kratos

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.OptimizationApplication.utilities.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.utilities.communicators.optimization_component_communicator import OptimizationComponentCommunicator
from KratosMultiphysics.OptimizationApplication.responses.mass_response_function import MassResponseFunction

class TestOptimizationInfo(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.optimization_info = OptimizationInfo()

        cls.response_function = MassResponseFunction(cls.model, Kratos.Parameters("""{"evaluated_model_part_names": ["test"]}"""), cls.optimization_info)
        cls.communicator = OptimizationComponentCommunicator(cls.optimization_info)
        cls.communicator.AddResponseFunction("mass", cls.response_function)

    def test_AdvanceStep(self):
        problem_data = self.communicator.GetProblemData()

        problem_data["mass_1"] = OptimizationInfo(3)
        initial_step = self.communicator.GetStep()

        problem_data["mass_1/int"] = 1
        self.assertEqual(self.communicator.GetStep(), initial_step)

        self.communicator.AdvanceStep()
        problem_data["mass_1/int"] = 2
        self.assertEqual(problem_data["mass_1/int", 1], 1)
        self.assertEqual(self.communicator.GetStep(), initial_step + 1)
        self.assertEqual(self.communicator.GetResponseFunction("mass"), self.response_function)

        self.communicator.AdvanceStep()
        problem_data["mass_1/int"] = 3
        self.assertEqual(problem_data["mass_1/int", 1], 2)
        self.assertEqual(problem_data["mass_1/int", 2], 1)
        self.assertEqual(self.communicator.GetStep(), initial_step + 2)
        self.assertEqual(self.communicator.GetResponseFunction("mass"), self.response_function)

        self.communicator.AdvanceStep()
        problem_data["mass_1/int"] = 4
        self.assertEqual(problem_data["mass_1/int", 1], 3)
        self.assertEqual(problem_data["mass_1/int", 2], 2)
        self.assertEqual(self.communicator.GetStep(), initial_step + 3)
        self.assertEqual(self.communicator.GetResponseFunction("mass"), self.response_function)

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()