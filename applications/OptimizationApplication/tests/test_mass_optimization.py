
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.OptimizationApplication.optimization_analysis import OptimizationAnalysis

class TestMassOptimization(kratos_unittest.TestCase):
    def test_MassOptimization(self):
        model = Kratos.Model()
        with kratos_unittest.WorkFolderScope("mass_optimization", __file__):
            with open("optimization_parameters.json", "r") as file_input:
                parameters = Kratos.Parameters(file_input.read())

            analysis = OptimizationAnalysis(model, parameters)
            analysis.Run()

            self.assertAlmostEqual(analysis.optimization_info.GetValue("problem_data/response_data/Structure/mass_Structure/value"), 13446.777989000835, 10)
            self.assertAlmostEqual(analysis.optimization_info.GetValue("problem_data/response_data/Structure/strain_energy/value"), 0.12021307691786137, 10)
            self.assertAlmostEqual(analysis.optimization_info.GetValue("problem_data/algorithm_data/info/projection_norm"), 0.7009336912136516, 10)
            self.assertAlmostEqual(analysis.optimization_info.GetValue("problem_data/algorithm_data/info/correction_norm"), 0.011216149756261142, 10)
            self.assertAlmostEqual(analysis.optimization_info.GetValue("problem_data/algorithm_data/info/search_direction_norm"), 0.12000523834406328, 10)

    @classmethod
    def tearDownClass(cls):
        with kratos_unittest.WorkFolderScope("mass_optimization", __file__):
            DeleteFileIfExisting("Structure.time")

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()