
import KratosMultiphysics as Kratos

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.buffered_dict import BufferedDict
from KratosMultiphysics.OptimizationApplication.responses.mass_response_function import MassResponseFunction
class TestOptimizationInfo(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.optimization_info = OptimizationProblem()

        cls.response_function = MassResponseFunction(cls.model, Kratos.Parameters("""{"evaluated_model_part_names": ["test"]}"""), cls.optimization_info)
        cls.optimization_info.AddResponse("mass", cls.response_function)

    def test_GetReponseData(self):
        mass_problem_data = self.optimization_info.GetReponseData("mass")
        self.assertEqual(mass_problem_data, self.optimization_info.GetReponseData("mass"))
        self.assertEqual(self.optimization_info.GetProblemDataContainer(), self.optimization_info.GetReponseData("mass").GetRoot())

    def test_AdvanceStep(self):
        mass_problem_data = self.optimization_info.GetReponseData("mass")
        mass_problem_data["mass_1"] = BufferedDict(3)

        initial_step = self.optimization_info.GetStep()

        mass_problem_data["mass_1/int"] = 1
        self.assertEqual(self.optimization_info.GetStep(), initial_step)

        self.optimization_info.AdvanceStep()
        mass_problem_data["mass_1/int"] = 2
        self.assertEqual(mass_problem_data["mass_1/int", 1], 1)
        self.assertEqual(self.optimization_info.GetStep(), initial_step + 1)
        self.assertEqual(self.optimization_info.GetResponse("mass"), self.response_function)

        self.optimization_info.AdvanceStep()
        mass_problem_data["mass_1/int"] = 3
        self.assertEqual(mass_problem_data["mass_1/int", 1], 2)
        self.assertEqual(mass_problem_data["mass_1/int", 2], 1)
        self.assertEqual(self.optimization_info.GetStep(), initial_step + 2)
        self.assertEqual(self.optimization_info.GetResponse("mass"), self.response_function)

        self.optimization_info.AdvanceStep()
        mass_problem_data["mass_1/int"] = 4
        self.assertEqual(mass_problem_data["mass_1/int", 1], 3)
        self.assertEqual(mass_problem_data["mass_1/int", 2], 2)
        self.assertEqual(self.optimization_info.GetStep(), initial_step + 3)
        self.assertEqual(self.optimization_info.GetResponse("mass"), self.response_function)

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()