import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as kratos_unittest

from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.kratos_utilities import DeleteDirectoryIfExisting
from KratosMultiphysics.OptimizationApplication.optimization_analysis import OptimizationAnalysis

class TestSystemIdentificationRestart(kratos_unittest.TestCase):
    """Restart with sensors and a filtered SysId control: sensors are re-derived (not checkpointed),
    and the primal/adjoint model parts come back together with the whole Model."""
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
        """Material field read off the mesh. The explicit filter's phi field is a delta accumulator
        reset to zero on Initialize(), so unlike the physical field it is not restart-invariant."""
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
            # the adjoint solve uses an iterative (amgcl) solver: compare relative to the field magnitude
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
