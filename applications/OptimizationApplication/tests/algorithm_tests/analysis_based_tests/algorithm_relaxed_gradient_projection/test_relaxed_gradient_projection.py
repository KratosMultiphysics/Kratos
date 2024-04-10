import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.OptimizationApplication.optimization_analysis import OptimizationAnalysis
from KratosMultiphysics.compare_two_files_check_process import CompareTwoFilesCheckProcess

# Temporaly failing
# @KratosUnittest.skipIf(True, "Temporaly Failing")
class TestRelaxedGradientProjectionAnalysis(KratosUnittest.TestCase):
    def test_gradient_projection_analysis(self):
        with KratosUnittest.WorkFolderScope(".", __file__):
            with open("optimization_parameters.json", "r") as file_input:
                parameters = Kratos.Parameters(file_input.read())

            model = Kratos.Model()
            analysis = OptimizationAnalysis(model, parameters)
            analysis.Run()

            CompareTwoFilesCheckProcess(Kratos.Parameters("""
            {
                "reference_file_name"   : "summary_orig.csv",
                "output_file_name"      : "summary.csv",
                "remove_output_file"    : true,
                "comparison_type"       : "csv_file",
                "tolerance"             : 1e-6,
                "relative_tolerance"    : 1e-6
            }""")).Execute()

    @classmethod
    def tearDownClass(cls) -> None:
        with KratosUnittest.WorkFolderScope(".", __file__):
            DeleteFileIfExisting("Structure.time")

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    KratosUnittest.main()

