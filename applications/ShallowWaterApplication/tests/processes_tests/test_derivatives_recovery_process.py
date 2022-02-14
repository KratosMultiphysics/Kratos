import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.ShallowWaterApplication.derivatives_recovery_process as derivatives_recovery_process
import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestDerivativesRecoveryProcess(KratosUnittest.TestCase):

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
                "compute_neighbors"        : true,
                "update_mesh_topology"     : false
            }
        }""")
        self.model = KM.Model()
        self.model_part = self.model.CreateModelPart('test_model_part')
        self.model_part.AddNodalSolutionStepVariable(SW.FIRST_DERIVATIVE_WEIGHTS)
        self.model_part.AddNodalSolutionStepVariable(SW.SECOND_DERIVATIVE_WEIGHTS)
        self.model_part.AddNodalSolutionStepVariable(KM.DISTANCE)
        self.model_part.AddNodalSolutionStepVariable(KM.DISTANCE_GRADIENT)
        self.model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        self.model_part.AddNodalSolutionStepVariable(KM.DETERMINANT_F)
        self.model_part.AddNodalSolutionStepVariable(KM.VELOCITY)
        self.model_part.AddNodalSolutionStepVariable(KM.VELOCITY_LAPLACIAN)
        self.model_part.ProcessInfo[KM.DOMAIN_SIZE] = 2
        KM.ModelPartIO(GetFilePath("model_part")).ReadModelPart(self.model_part)
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(KM.DISTANCE, node.X)
            node.SetSolutionStepValue(KM.DISPLACEMENT, [node.X, node.Y, 0])
            node.SetSolutionStepValue(KM.VELOCITY, [node.X**2, node.Y**2, 0])
        process = derivatives_recovery_process.Factory(settings, self.model)
        process.Check()
        process.ExecuteInitialize()
        first_node = self.model_part.Nodes.__iter__().__next__()
        self.assertVectorAlmostEqual(first_node.GetSolutionStepValue(KM.DISTANCE_GRADIENT), [0, 0, 0])
        self.assertVectorAlmostEqual(first_node.GetSolutionStepValue(KM.VELOCITY_LAPLACIAN), [0, 0, 0])
        self.assertAlmostEqual(first_node.GetSolutionStepValue(KM.DETERMINANT_F), 0)
        process.ExecuteInitializeSolutionStep()
        self.assertVectorAlmostEqual(first_node.GetSolutionStepValue(KM.DISTANCE_GRADIENT), [1, 0, 0])
        self.assertVectorAlmostEqual(first_node.GetSolutionStepValue(KM.VELOCITY_LAPLACIAN), [0, 0, 0])
        self.assertAlmostEqual(first_node.GetSolutionStepValue(KM.DETERMINANT_F), 0)
        process.ExecuteFinalizeSolutionStep()
        self.assertVectorAlmostEqual(first_node.GetSolutionStepValue(KM.DISTANCE_GRADIENT), [1, 0, 0])
        self.assertVectorAlmostEqual(first_node.GetSolutionStepValue(KM.VELOCITY_LAPLACIAN), [2, 2, 0])
        self.assertAlmostEqual(first_node.GetSolutionStepValue(KM.DETERMINANT_F), 2)

if __name__ == "__main__":
    KratosUnittest.main()
