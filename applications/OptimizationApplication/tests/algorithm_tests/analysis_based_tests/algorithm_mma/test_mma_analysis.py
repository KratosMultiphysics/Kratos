import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.OptimizationApplication.optimization_analysis import OptimizationAnalysis
from KratosMultiphysics.compare_two_files_check_process import CompareTwoFilesCheckProcess


@KratosUnittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
class TestMMAAnalysis(KratosUnittest.TestCase):
    """Minimises the mass of a shell plate (THICKNESS control, exact gradients) s.t. a strain energy limit."""

    def _RunAndCompare(self, variant: str) -> None:
        with KratosUnittest.WorkFolderScope(".", __file__):
            with open(f"optimization_parameters_{variant}.json", "r") as file_input:
                parameters = Kratos.Parameters(file_input.read())

            model = Kratos.Model()
            analysis = OptimizationAnalysis(model, parameters)
            analysis.Run()

            # Linear solver (pardiso/sparse_lu/skyline) and thread count change the MMA slack variable by up to 2e-5 and
            # the near-zero std constraint values (strain energy ~4e4 minus its limit) by ~1e-3 absolute.
            CompareTwoFilesCheckProcess(Kratos.Parameters("""
            {
                "reference_file_name"   : "summary_orig_""" + variant + """.csv",
                "output_file_name"      : "summary_""" + variant + """.csv",
                "remove_output_file"    : true,
                "comparison_type"       : "csv_file",
                "tolerance"             : 1e-2,
                "relative_tolerance"    : 1e-5
            }""")).Execute()

    def test_mma_analysis(self):
        self._RunAndCompare("mma")

    def test_gcmma_analysis(self):
        self._RunAndCompare("gcmma")

    @classmethod
    def tearDownClass(cls) -> None:
        with KratosUnittest.WorkFolderScope(".", __file__):
            DeleteFileIfExisting("shell_thin.time")
            DeleteFileIfExisting("Structure.time")


if __name__ == "__main__":
    KratosUnittest.main()
