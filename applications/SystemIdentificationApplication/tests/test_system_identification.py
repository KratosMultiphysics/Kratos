import subprocess
import os
import sys
import numpy
import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as UnitTest
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
from KratosMultiphysics.OptimizationApplication.optimization_analysis import OptimizationAnalysis
from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.compare_two_files_check_process import CompareTwoFilesCheckProcess

class TestSystemIdentification(UnitTest.TestCase):
    def __RunSingleThreadedProcess(self, script_path):
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = "1"
        env["OPENBLAS_NUM_THREADS"] = "1"
        env["NUMEXPR_NUM_THREADS"] = "1"
        if any(k.startswith("MKL") for k in env):
            env["MKL_NUM_THREADS"] = "1"

        subprocess.run(
            [sys.executable, script_path],
            env=env,
            check=True
        )

    def test_DamagedSystem(self):
        self.addCleanup(DeleteFileIfExisting, "auxiliary_files/damaged_problem/measured_data.csv")

        self._run_single_threaded_process("auxiliary_files/damaged_problem/MainKratos.py")

        data = numpy.loadtxt("auxiliary_files/damaged_problem/measured_data.csv", comments="#", usecols=[0,3,4,5,6], delimiter=",")
        ref_data = numpy.loadtxt("auxiliary_files/damaged_problem/measured_data_ref.csv", comments="#", usecols=[0,3,4,5,6], delimiter=",")
        numpy.testing.assert_allclose(data, ref_data, rtol=1e-13, atol=1e-13, verbose=True)

    def test_SystemIdentification(self):
        self._run_single_threaded_process("auxiliary_files/system_identification/MainKratos.py")

        params = Kratos.Parameters("""{
            "reference_file_name"   : "auxiliary_files/system_identification_summary_ref.csv",
            "output_file_name"      : "auxiliary_files/summary.csv",
            "remove_output_file"    : true,
            "comparison_type"       : "csv_file",
            "tolerance"             : 1e-11,
            "relative_tolerance"    : 1e-11,
            "dimension"             : 3
        }""")
        CompareTwoFilesCheckProcess(params).Execute()

    def test_SystemIdentificationPNorm(self):
        self._run_single_threaded_process("auxiliary_files/system_identification_p_norm/MainKratos.py")

        params = Kratos.Parameters("""{
            "reference_file_name"   : "auxiliary_files/system_identification_p_norm_summary_ref.csv",
            "output_file_name"      : "auxiliary_files/summary_p_norm.csv",
            "remove_output_file"    : true,
            "comparison_type"       : "csv_file",
            "tolerance"             : 1e-11,
            "relative_tolerance"    : 1e-11,
            "dimension"             : 3
        }""")
        CompareTwoFilesCheckProcess(params).Execute()

if __name__ == '__main__':
    UnitTest.main()

