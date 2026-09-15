import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as kratos_unittest

from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.kratos_utilities import DeleteDirectoryIfExisting
from KratosMultiphysics.OptimizationApplication.optimization_analysis import OptimizationAnalysis

class TestSystemIdentificationRestart(kratos_unittest.TestCase):
    """Restart test for a full system-identification setup (sensors + a filtered/phi-space
    material control), mirroring OptimizationApplication's algorithm_steepest_descent_restart
    test but exercising two things that test cannot, because it uses neither sensors nor a
    filtered control:

    1. Sensors are never written to the restart checkpoint (they are re-derived from
       "list_of_sensors"/measurement-file settings on every run, restart or not -- see
       SensorModelPartController.ImportModelPart) -- optimization_problem_restart_output_process
       must skip them rather than crash (see ConvertLeafToRestartData), and the resumed run's
       damage_response must still match the reference run, proving the re-derived sensors are
       equivalent to the reference run's.
    2. MaterialPropertiesControl (SI) is Initialize()'d twice on a resumed run -- once explicitly
       by optimization_problem_restart_input_process before restoring tensor data, once more via
       Algorithm.Initialize() -- and its "element_specific_properties_created" ModelPart status
       flag (material_properties_control.py) must make the second call a no-op; otherwise entity
       specific properties would be created twice and the resumed run would diverge (or crash).
       "output_all_fields": true is set below so the control's ComponentDataView also holds
       several diagnostic-only tensors (not just the design), stressing the restart process's
       "which tensor is actually the design" gating (see RestoreBufferedDict).
    """
    restart_files_path = "Optimization_Restart_test"

    def _ReadParameters(self) -> Kratos.Parameters:
        with open("auxiliary_files/system_identification_restart/optimization_parameters.json", "r") as file_input:
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

    def test_system_identification_restart(self):
        with kratos_unittest.WorkFolderScope(".", __file__):
            # reference run: solve iterations 0..2 in one uninterrupted process.
            reference_parameters = self._ReadParameters()
            self._SetMaxIter(reference_parameters, 3)
            reference_analysis = self._RunToConvergence(reference_parameters)
            reference_algorithm = reference_analysis.GetAlgorithm()
            reference_obj_value = reference_algorithm.GetOptimizedObjectiveValue()
            reference_control_field = list(reference_algorithm.GetCurrentControlField().data)

            # checkpoint run: solve iterations 0..1 only, writing a restart checkpoint every step.
            checkpoint_parameters = self._ReadParameters()
            self._SetMaxIter(checkpoint_parameters, 2)
            self._AddProcess(checkpoint_parameters, "output_processes", Kratos.Parameters("""{
                "type"    : "optimization_problem_restart_output_process",
                "module"  : "KratosMultiphysics.OptimizationApplication.processes",
                "settings": {
                    "restart_files_path"    : \"""" + self.restart_files_path + """\",
                    "restart_save_frequency": 1
                }
            }"""))
            self._RunToConvergence(checkpoint_parameters)

            # resume run: restore the checkpoint from iteration 2 and run one more iteration to 3.
            # This is also where MaterialPropertiesControl.Initialize() runs twice (see class
            # docstring): once from optimization_problem_restart_input_process.ExecuteInitialize,
            # once more from Algorithm.Initialize() right after.
            resume_parameters = self._ReadParameters()
            self._SetMaxIter(resume_parameters, 3)
            self._AddProcess(resume_parameters, "auxiliary_processes", Kratos.Parameters("""{
                "type"    : "optimization_problem_restart_input_process",
                "module"  : "KratosMultiphysics.OptimizationApplication.processes",
                "settings": {
                    "restart_files_path": \"""" + self.restart_files_path + """\",
                    "restart_load_step" : 2
                }
            }"""))
            resume_analysis = self._RunToConvergence(resume_parameters)
            resume_algorithm = resume_analysis.GetAlgorithm()
            resume_obj_value = resume_algorithm.GetOptimizedObjectiveValue()
            resume_control_field = list(resume_algorithm.GetCurrentControlField().data)

            self.assertAlmostEqual(reference_obj_value, resume_obj_value, places=9)
            # GetCurrentControlField() here is the *phi-space* (filtered/clamped) field, not the
            # raw physical one, so a places-based comparison isn't meaningful either way -- compare
            # with a delta scaled to the field's own magnitude instead. Unlike
            # algorithm_steepest_descent_restart's linear, direct (unfiltered, unbounded) primal
            # solve, this response runs an adjoint FEM solve through an iterative (amgcl) linear
            # solver each iteration, so run-to-run reproducibility has a higher noise floor (order
            # 1e-8 relative here, vs. 1e-9 there) -- still tight enough to catch a materially wrong
            # restart restore, e.g. sensors/measurements not reproduced identically on resume, or
            # entity specific properties created twice (this delta failed at ~3e-1 relative, not
            # 1e-8, when this test first caught the control's Initialize()-is-idempotent bug fixed
            # alongside it -- see MaterialPropertiesControl.Initialize()).
            max_abs_reference_value = max(abs(value) for value in reference_control_field)
            self.assertVectorAlmostEqual(reference_control_field, resume_control_field, places=None, delta=max_abs_reference_value * 1e-6)

    @classmethod
    def tearDownClass(cls) -> None:
        with kratos_unittest.WorkFolderScope(".", __file__):
            DeleteFileIfExisting("Structure.time")
            DeleteFileIfExisting("AdjointStructure.time")
            DeleteDirectoryIfExisting(cls.restart_files_path)

if __name__ == "__main__":
    kratos_unittest.main()
