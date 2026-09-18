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
       SensorModelPartController.ImportModelPart) -- OptimizationProblemRestartOutput must skip
       them rather than crash (they have no container to launder through), and the resumed run's
       damage_response must still match the reference run, proving the re-derived sensors are
       equivalent to the reference run's.
    2. Restart now saves/restores the whole Kratos.Model in one shot, so "AdjointStructure"
       (independently mdpa-imported) and "Structure" (derived from it via
       connectivity_preserving_model_part_controller, so it has no independent import of its own)
       come back together, both carrying the MODEL_PART_STATUS they had at save time. This is what
       actually makes MaterialPropertiesControl.Initialize()'s primal/adjoint consistency check
       (material_properties_control.py) pass on resume -- restoring only one side used to trip it.
       "output_all_fields": true is set below so the control's ComponentDataView also holds
       several diagnostic-only tensors (not just the design), stressing the restart process's
       "which tensor is actually the design" gating.

    """
    restart_files_path = "Optimization_Restart_test"

    def _ReadParameters(self) -> Kratos.Parameters:
        with open("auxiliary_files/system_identification_restart/optimization_parameters.json", "r") as file_input:
            return Kratos.Parameters(file_input.read())

    def _SetMaxIter(self, parameters: Kratos.Parameters, max_iter: int) -> None:
        parameters["algorithm_settings"]["settings"]["conv_settings"]["max_iter"].SetInt(max_iter)

    def _RunToConvergence(self, parameters: Kratos.Parameters) -> OptimizationAnalysis:
        model = Kratos.Model()
        analysis = OptimizationAnalysis(model, parameters)
        analysis.Run()
        return analysis

    def _GetPhysicalField(self, analysis: OptimizationAnalysis) -> 'list[float]':
        """The identified material field as it sits on the mesh -- the restart-invariant quantity.

        GetCurrentControlField() is deliberately not used: with an explicit filter the control
        (phi-space) field is a pure delta accumulator -- Update() applies only
        new_control_field - control_phi_field, and ExplicitFilter.UnfilterField() seeds it at zero
        on every Initialize() -- so a resumed run legitimately carries a different absolute phi
        value for the same physical field.
        """
        control, = analysis.GetAlgorithm().master_control.GetListOfControls()
        return list(control.GetPhysicalField().data)

    def test_system_identification_restart(self):
        with kratos_unittest.WorkFolderScope(".", __file__):
            # reference run: solve iterations 0..2 in one uninterrupted process.
            reference_parameters = self._ReadParameters()
            self._SetMaxIter(reference_parameters, 3)
            reference_analysis = self._RunToConvergence(reference_parameters)
            reference_algorithm = reference_analysis.GetAlgorithm()
            reference_obj_value = reference_algorithm.GetOptimizedObjectiveValue()
            reference_physical_field = self._GetPhysicalField(reference_analysis)

            # checkpoint run: solve iterations 0..1 only, writing a restart checkpoint every step.
            checkpoint_parameters = self._ReadParameters()
            self._SetMaxIter(checkpoint_parameters, 2)
            checkpoint_parameters.AddValue("restart_settings", Kratos.Parameters("""{
                "save_restart"           : true,
                "restart_file_name"      : \"""" + self.restart_files_path + """/restart_<step>.pkl",
                "restart_save_frequency" : 1
            }"""))
            self._RunToConvergence(checkpoint_parameters)

            # resume run: restore the checkpoint from iteration 2 and run one more iteration to 3.
            resume_parameters = self._ReadParameters()
            self._SetMaxIter(resume_parameters, 3)
            resume_parameters.AddValue("restart_settings", Kratos.Parameters("""{
                "load_restart"       : true,
                "restart_file_name"  : \"""" + self.restart_files_path + """/restart_<step>.pkl",
                "restart_load_step"  : 2
            }"""))
            resume_analysis = self._RunToConvergence(resume_parameters)
            resume_algorithm = resume_analysis.GetAlgorithm()
            resume_obj_value = resume_algorithm.GetOptimizedObjectiveValue()
            resume_physical_field = self._GetPhysicalField(resume_analysis)

            self.assertAlmostEqual(reference_obj_value, resume_obj_value, places=9)
            # This response runs an adjoint FEM solve through an iterative (amgcl) linear solver
            # each iteration, so run-to-run reproducibility has a higher noise floor than
            # algorithm_steepest_descent_restart's linear, direct solve -- compare with a delta
            # scaled to the field's own magnitude instead of an absolute places=. Still tight
            # enough to catch a materially wrong restore, e.g. sensors/measurements not reproduced
            # identically on resume, or entity specific properties created twice (this delta failed
            # at ~3e-1 relative when this test first caught the control's Initialize()-is-idempotent
            # bug fixed alongside it -- see MaterialPropertiesControl.Initialize()).
            max_abs_reference_value = max(abs(value) for value in reference_physical_field)
            self.assertVectorAlmostEqual(reference_physical_field, resume_physical_field, places=None, delta=max_abs_reference_value * 1e-6)

    @classmethod
    def tearDownClass(cls) -> None:
        with kratos_unittest.WorkFolderScope(".", __file__):
            DeleteFileIfExisting("Structure.time")
            DeleteFileIfExisting("AdjointStructure.time")
            DeleteDirectoryIfExisting(cls.restart_files_path)

if __name__ == "__main__":
    kratos_unittest.main()
