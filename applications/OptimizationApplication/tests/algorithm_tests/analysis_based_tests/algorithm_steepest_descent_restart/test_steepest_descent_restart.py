import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as kratos_unittest

from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.kratos_utilities import DeleteDirectoryIfExisting
from KratosMultiphysics.OptimizationApplication.optimization_analysis import OptimizationAnalysis

class TestSteepestDescentRestart(kratos_unittest.TestCase):
    restart_files_path = "Optimization_Restart_test"

    def _ReadParameters(self) -> Kratos.Parameters:
        with open("optimization_parameters.json", "r") as file_input:
            return Kratos.Parameters(file_input.read())

    def _SetMaxIter(self, parameters: Kratos.Parameters, max_iter: int) -> None:
        parameters["algorithm_settings"]["settings"]["conv_settings"]["max_iter"].SetInt(max_iter)

    def _AddProcess(self, parameters: Kratos.Parameters, process_category: str, process_settings: Kratos.Parameters) -> None:
        optimization_data_processes = parameters["processes"]["optimization_data_processes"]
        if not optimization_data_processes.Has(process_category):
            optimization_data_processes.AddEmptyArray(process_category)
        optimization_data_processes[process_category].Append(process_settings)

    def _RunToConvergence(self, parameters: Kratos.Parameters) -> OptimizationAnalysis:
        model = Kratos.Model()
        analysis = OptimizationAnalysis(model, parameters)
        analysis.Run()
        return analysis

    def test_steepest_descent_restart(self):
        with kratos_unittest.WorkFolderScope(".", __file__):
            # reference run: solve iterations 0..6 in one uninterrupted process.
            reference_parameters = self._ReadParameters()
            self._SetMaxIter(reference_parameters, 6)
            reference_analysis = self._RunToConvergence(reference_parameters)
            reference_algorithm = reference_analysis.GetAlgorithm()
            reference_obj_value = reference_algorithm.GetOptimizedObjectiveValue()
            reference_control_field = list(reference_algorithm.GetCurrentControlField().data)

            # checkpoint run: solve iterations 0..5 only, writing a restart checkpoint every step.
            checkpoint_parameters = self._ReadParameters()
            self._SetMaxIter(checkpoint_parameters, 5)
            self._AddProcess(checkpoint_parameters, "output_processes", Kratos.Parameters("""{
                "type"    : "optimization_problem_restart_output_process",
                "module"  : "KratosMultiphysics.OptimizationApplication.processes",
                "settings": {
                    "restart_files_path"    : \"""" + self.restart_files_path + """\",
                    "restart_save_frequency": 1
                }
            }"""))
            self._RunToConvergence(checkpoint_parameters)

            # resume run: restore the checkpoint from iteration 5 and run one more iteration to 6.
            resume_parameters = self._ReadParameters()
            self._SetMaxIter(resume_parameters, 6)
            self._AddProcess(resume_parameters, "auxiliary_processes", Kratos.Parameters("""{
                "type"    : "optimization_problem_restart_input_process",
                "module"  : "KratosMultiphysics.OptimizationApplication.processes",
                "settings": {
                    "restart_files_path": \"""" + self.restart_files_path + """\",
                    "restart_load_step" : 5
                }
            }"""))
            resume_analysis = self._RunToConvergence(resume_parameters)
            resume_algorithm = resume_analysis.GetAlgorithm()
            resume_obj_value = resume_algorithm.GetOptimizedObjectiveValue()
            resume_control_field = list(resume_algorithm.GetCurrentControlField().data)

            self.assertAlmostEqual(reference_obj_value, resume_obj_value, places=9)
            # the control field holds YOUNG_MODULUS-scale values (~1e10), so an absolute
            # places-based comparison would demand ~19-20 significant digits -- well beyond
            # double precision. Compare with a delta scaled to the field's own magnitude instead,
            # loose enough to absorb floating-point noise (e.g. parallel reduction order) while
            # still catching a materially wrong restart restore.
            max_abs_reference_value = max(abs(value) for value in reference_control_field)
            self.assertVectorAlmostEqual(reference_control_field, resume_control_field, places=None, delta=max_abs_reference_value * 1e-9)

    @classmethod
    def tearDownClass(cls) -> None:
        with kratos_unittest.WorkFolderScope(".", __file__):
            DeleteFileIfExisting("Structure.time")
            DeleteDirectoryIfExisting(cls.restart_files_path)

if __name__ == "__main__":
    kratos_unittest.main()
