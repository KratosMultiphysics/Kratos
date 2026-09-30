import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.OptimizationApplication.optimization_analysis import OptimizationAnalysis
from KratosMultiphysics.compare_two_files_check_process import CompareTwoFilesCheckProcess


class TestMMAAnalysis(KratosUnittest.TestCase):
    def setUp(self):
        # MMA's asymptote/move-limit heuristics react very sensitively to this problem's
        # design space (a shell thickness control whose raw, unfiltered representation is
        # poorly scaled), so tiny floating-point differences from OpenMP's thread-order-
        # dependent reduction sums get amplified across outer iterations. Pin to a single
        # thread so the reference CSV comparison is reproducible regardless of the machine's
        # core count.
        self.previous_threads = Kratos.ParallelUtilities.GetNumThreads()
        Kratos.ParallelUtilities.SetNumThreads(1)

    def tearDown(self):
        Kratos.ParallelUtilities.SetNumThreads(self.previous_threads)

    def test_mma_analysis(self):
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
    KratosUnittest.main()
