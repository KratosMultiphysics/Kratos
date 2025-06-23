import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as kratos_unittest

from KratosMultiphysics.OptimizationApplication.optimization_analysis import OptimizationAnalysis
from KratosMultiphysics.compare_two_files_check_process import CompareTwoFilesCheckProcess
from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting


import csv, os

try:
    import nlopt
except ImportError:
    nlopt = None

@kratos_unittest.skipIf(nlopt==None, "Missing nlopt python libraries ")
class TestNLOPTOptimizer(kratos_unittest.TestCase):
    def test_MMA_optimizer(self):
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
