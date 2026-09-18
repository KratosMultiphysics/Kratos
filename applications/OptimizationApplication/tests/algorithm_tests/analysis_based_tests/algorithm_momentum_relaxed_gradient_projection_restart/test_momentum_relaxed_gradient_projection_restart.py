import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as kratos_unittest

from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.kratos_utilities import DeleteDirectoryIfExisting
from KratosMultiphysics.OptimizationApplication.optimization_analysis import OptimizationAnalysis

class TestMomentumRelaxedGradientProjectionRestart(kratos_unittest.TestCase):
    """Restart test for the momentum-based relaxed gradient projection algorithm, mirroring
    algorithm_nesterov_accelerated_gradient_restart.

    AlgorithmMomentumRelaxedGradientProjection.ComputeControlUpdate implements the same
    Nesterov-style momentum recurrence as AlgorithmNesterovAcceleratedGradient (see that
    algorithm's docstring/comments), and had the identical bug: momentum was kept in a plain
    Python attribute (self.prev_update) invisible to the restart snapshot, so a resumed run
    silently lost it. Fixed the same way -- bookkept via ComponentDataView under "momentum"
    instead -- and verified here the same way.

    Restart now always saves/restores the whole Kratos.Model (no per-model-part toggle), so any
    NEIGHBOUR_ELEMENTS/NEIGHBOUR_CONDITIONS left on the nodes gets serialized too, which used to
    segfault; OptAppModelPartUtils.Clear/RestoreNeighbourEntitiesData (used internally by
    OptimizationProblemRestartOutput) clears it around the save.

    The control filters explicitly on purpose -- do not switch it back to "implicit_filter".
    The implicit (Helmholtz/Sobolev) filter builds its own auxiliary root model part of
    HelmholtzElements, which whole-Model save now serializes like any other, and
    Serializer::Create<HelmholtzElement<...>> crashes on load inside
    EntityCalculationUtils::CreateSolidGeometry -- a pre-existing C++ deserialization bug in the
    element, out of scope here. So implicit filtering has no restart coverage for now; nothing this
    test checks (momentum surviving a checkpoint/resume) depends on which filter is used.
    """
    restart_files_path = "Optimization_Restart_test"

    def _ReadParameters(self) -> Kratos.Parameters:
        with open("optimization_parameters.json", "r") as file_input:
            return Kratos.Parameters(file_input.read())

    def _SetMaxIter(self, parameters: Kratos.Parameters, max_iter: int) -> None:
        parameters["algorithm_settings"]["settings"]["conv_settings"]["max_iter"].SetInt(max_iter)

    def _RunToConvergence(self, parameters: Kratos.Parameters) -> OptimizationAnalysis:
        model = Kratos.Model()
        analysis = OptimizationAnalysis(model, parameters)
        analysis.Run()
        return analysis

    def _GetPhysicalField(self, analysis: OptimizationAnalysis) -> 'list[float]':
        """The design as it sits on the mesh (element THICKNESS) -- the restart-invariant quantity.

        GetCurrentControlField() is deliberately not used: with an explicit filter the control
        (phi-space) field is a pure delta accumulator -- Update() applies only
        new_control_field - control_phi_field, and ExplicitFilter.UnfilterField() seeds it at zero
        on every Initialize() -- so a resumed run legitimately carries a different absolute phi
        value for the same physical design.
        """
        control, = analysis.GetAlgorithm().master_control.GetListOfControls()
        return list(control.GetPhysicalField().data)

    def test_momentum_relaxed_gradient_projection_restart(self):
        with kratos_unittest.WorkFolderScope(".", __file__):
            # reference run: solve iterations 0..3 in one uninterrupted process.
            reference_parameters = self._ReadParameters()
            self._SetMaxIter(reference_parameters, 4)
            reference_analysis = self._RunToConvergence(reference_parameters)
            reference_algorithm = reference_analysis.GetAlgorithm()
            reference_obj_value = reference_algorithm.GetOptimizedObjectiveValue()
            reference_physical_field = self._GetPhysicalField(reference_analysis)

            # checkpoint run: solve iterations 0..2 only (momentum is already non-trivial by
            # then), writing a restart checkpoint every step.
            checkpoint_parameters = self._ReadParameters()
            self._SetMaxIter(checkpoint_parameters, 3)
            checkpoint_parameters.AddValue("restart_settings", Kratos.Parameters("""{
                "save_restart"           : true,
                "restart_file_name"      : \"""" + self.restart_files_path + """/restart_<step>.pkl",
                "restart_save_frequency" : 1
            }"""))
            self._RunToConvergence(checkpoint_parameters)

            # resume run: restore the checkpoint from iteration 3 and run one more iteration to 4.
            resume_parameters = self._ReadParameters()
            self._SetMaxIter(resume_parameters, 4)
            resume_parameters.AddValue("restart_settings", Kratos.Parameters("""{
                "load_restart"       : true,
                "restart_file_name"  : \"""" + self.restart_files_path + """/restart_<step>.pkl",
                "restart_load_step"  : 3
            }"""))
            resume_analysis = self._RunToConvergence(resume_parameters)
            resume_algorithm = resume_analysis.GetAlgorithm()
            resume_obj_value = resume_algorithm.GetOptimizedObjectiveValue()
            resume_physical_field = self._GetPhysicalField(resume_analysis)

            # This control is filtered and passed through a nonlinear sigmoidal projection, and the
            # constraint runs an FEM solve each iteration, so reproducibility has a small but
            # non-zero noise floor (observed ~1e-9 relative on both quantities below) -- compare
            # relative to each field's own magnitude rather than with an absolute places=.
            self.assertAlmostEqual(reference_obj_value, resume_obj_value, delta=abs(reference_obj_value) * 1e-7)
            max_abs_reference_value = max(abs(value) for value in reference_physical_field)
            self.assertVectorAlmostEqual(reference_physical_field, resume_physical_field, places=None, delta=max_abs_reference_value * 1e-6)

    @classmethod
    def tearDownClass(cls) -> None:
        with kratos_unittest.WorkFolderScope(".", __file__):
            DeleteFileIfExisting("Structure.time")
            DeleteDirectoryIfExisting(cls.restart_files_path)

if __name__ == "__main__":
    kratos_unittest.main()
