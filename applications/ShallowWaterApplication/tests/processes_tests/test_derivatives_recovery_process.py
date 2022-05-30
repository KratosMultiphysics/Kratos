import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.ShallowWaterApplication.derivatives_recovery_process as derivatives_recovery_process
from unittest.mock import patch


class TestDerivativesRecoveryProcess(KratosUnittest.TestCase):

    def setUp(self):
        self.model = KM.Model()
        self.model_part = self.model.CreateModelPart('test_model_part')
        self.model_part.ProcessInfo[KM.DOMAIN_SIZE] = 2
        self.patcher1 = patch('KratosMultiphysics.ShallowWaterApplication.DerivativesRecoveryUtility2D.Check', autospec=True)
        self.patcher2 = patch('KratosMultiphysics.ShallowWaterApplication.DerivativesRecoveryUtility2D.RecoverGradient', autospec=True)
        self.patcher3 = patch('KratosMultiphysics.ShallowWaterApplication.DerivativesRecoveryUtility2D.RecoverDivergence', autospec=True)
        self.patcher4 = patch('KratosMultiphysics.ShallowWaterApplication.DerivativesRecoveryUtility2D.RecoverLaplacian', autospec=True)
        self.RecoveryCheck = self.patcher1.start()
        self.RecoverGradient = self.patcher2.start()
        self.RecoverDivergence = self.patcher3.start()
        self.RecoverLaplacian = self.patcher4.start()

    def tearDown(self):
        self.patcher1.stop()
        self.patcher2.stop()
        self.patcher3.stop()
        self.patcher4.stop()

    def test_DerivativesRecoveryProcess(self):
        settings = KM.Parameters("""
        {
            "Parameters": {
                "model_part_name": "test_model_part",
                "list_of_operations"       : [{
                    "operation"           : "gradient",
                    "primitive_variable"  : "DISTANCE",
                    "derivative_variable" : "DISTANCE_GRADIENT",
                    "buffer_step"         : 0,
                    "process_step"        : "ExecuteInitializeSolutionStep"
                },{
                    "operation"           : "divergence",
                    "primitive_variable"  : "DISPLACEMENT",
                    "derivative_variable" : "DETERMINANT_F",
                    "buffer_step"         : 0,
                    "process_step"        : "ExecuteFinalizeSolutionStep"
                },{
                    "operation"           : "laplacian",
                    "primitive_variable"  : "VELOCITY",
                    "derivative_variable" : "VELOCITY_LAPLACIAN",
                    "buffer_step"         : 0,
                    "process_step"        : "ExecuteFinalizeSolutionStep"
                }],
                "compute_neighbors"        : false,
                "update_mesh_topology"     : false
            }
        }""")
        process = derivatives_recovery_process.Factory(settings, self.model)
        process.Check()
        self.assertEqual(self.RecoveryCheck.call_count, 1)
        self.assertEqual(self.RecoverGradient.call_count, 0)
        self.assertEqual(self.RecoverDivergence.call_count, 0)
        self.assertEqual(self.RecoverLaplacian.call_count, 0)
        process.ExecuteInitializeSolutionStep()
        self.assertEqual(self.RecoveryCheck.call_count, 1)
        self.assertEqual(self.RecoverGradient.call_count, 1)
        self.assertEqual(self.RecoverDivergence.call_count, 0)
        self.assertEqual(self.RecoverLaplacian.call_count, 0)
        process.ExecuteFinalizeSolutionStep()
        self.assertEqual(self.RecoveryCheck.call_count, 1)
        self.assertEqual(self.RecoverGradient.call_count, 1)
        self.assertEqual(self.RecoverDivergence.call_count, 1)
        self.assertEqual(self.RecoverLaplacian.call_count, 1)

if __name__ == "__main__":
    KratosUnittest.main()
