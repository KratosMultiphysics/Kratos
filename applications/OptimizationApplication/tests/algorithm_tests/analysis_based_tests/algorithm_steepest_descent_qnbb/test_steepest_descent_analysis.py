import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as kratos_unittest

from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.OptimizationApplication.optimization_analysis import OptimizationAnalysis
from KratosMultiphysics.compare_two_files_check_process import CompareTwoFilesCheckProcess

class TestQNBBSteepestDescentAnalysis(kratos_unittest.TestCase):
    def test_steepest_descent_analysis(self):
        with kratos_unittest.WorkFolderScope(".", __file__):
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
                "comparison_type"       : "deterministic"
            }""")).Execute()

    @classmethod
    def tearDownClass(cls) -> None:
        with kratos_unittest.WorkFolderScope(".", __file__):
            DeleteFileIfExisting("Structure.time")

if __name__ == "__main__":
    kratos_unittest.main()

