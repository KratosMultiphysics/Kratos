import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as kratos_unittest

from KratosMultiphysics.OptimizationApplication.optimization_analysis import OptimizationAnalysis
from KratosMultiphysics.compare_two_files_check_process import CompareTwoFilesCheckProcess

import csv, os

try:
    import scipy
except ImportError:
    scipy = None

@kratos_unittest.skipIf(scipy==None, "Missing scipy python libraries ")
class TestSciPyOptimizer(kratos_unittest.TestCase):
    def test_MMA_optimizer(self):
        with kratos_unittest.WorkFolderScope(".", __file__):
            with open("optimization_parameters.json", "r") as file_input:
                parameters = Kratos.Parameters(file_input.read())
            model = Kratos.Model()
            analysis = OptimizationAnalysis(model, parameters)
            analysis.Run()
            CompareTwoFilesCheckProcess(Kratos.Parameters("""
            {
                "reference_file_name"   : "summary_ref.csv",
                "output_file_name"      : "summary.csv",
                "remove_output_file"    : true,
                "comparison_type"       : "csv_file",
                "tolerance"             : 1e-6,
                "relative_tolerance"    : 1e-6
            }""")).Execute()

if __name__ == "__main__":
    kratos_unittest.main()
