import os
import math

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory

# from KratosMultiphysics.gid_output_process import GiDOutputProcess

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestEdgeBasedGradientRecoveryProcess(KratosUnittest.TestCase):

    def tearDown(self):
        # Remove the .time file
        try:
            os.remove('levelset_convection_process_mesh.time')
        except :
            pass

    def _SetUpTestModelPart(self, IsOriginHistorical, IsGradientHistorical, OriginVar = None, GradientVar = None):
        # Create the background mesh
        self.model = KratosMultiphysics.Model()
        model_part = self.model.CreateModelPart("TestModelPart")
        if IsOriginHistorical:
            model_part.AddNodalSolutionStepVariable(OriginVar)
        if IsGradientHistorical:
            model_part.AddNodalSolutionStepVariable(GradientVar)

        problem_domain = KratosMultiphysics.Quadrilateral2D4(
            KratosMultiphysics.Node(1, 0.0, 0.0, 0.0),
            KratosMultiphysics.Node(2, 0.0, 1.0, 0.0),
            KratosMultiphysics.Node(3, 1.0, 1.0, 0.0),
            KratosMultiphysics.Node(4, 1.0, 0.0, 0.0))
        parameters = KratosMultiphysics.Parameters("{}")
        parameters.AddEmptyValue("number_of_divisions").SetInt(20)
        parameters.AddEmptyValue("element_name").SetString("Element2D3N")
        parameters.AddEmptyValue("create_skin_sub_model_part").SetBool(False)
        # parameters.AddEmptyValue("condition_name").SetString("LineCondition2D2N")

        KratosMultiphysics.StructuredMeshGeneratorProcess(problem_domain, model_part, parameters).Execute()

    def testScalarHistoricalAndHistorical(self):
        # Set up the test model part
        self._SetUpTestModelPart(True, True, KratosMultiphysics.DISTANCE, KratosMultiphysics.DISTANCE_GRADIENT)

        # Set the origin variable field
        test_model_part = self.model.GetModelPart("TestModelPart")
        for node in test_model_part.Nodes:
            value = node.X**2 - math.sinh(4*node.X)/math.sinh(4)
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, value)

        # Set the gradient recovery process
        linear_solver_settings = KratosMultiphysics.Parameters("""{
            "solver_type" : "skyline_lu_factorization"
        }""")
        linear_solver = linear_solver_factory.ConstructSolver(linear_solver_settings)

        gradient_process_settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "TestModelPart",
            "origin_variable" : "DISTANCE",
            "gradient_variable" : "DISTANCE_GRADIENT",
            "gradient_penalty_coefficient" : 1.0e-6,
            "calculate_nodal_neighbours" : true,
            "is_historical_origin_variable" : true,
            "is_historical_gradient_variable": true
        }""")
        KratosMultiphysics.EdgeBasedGradientRecoveryProcessScalar(
            self.model,
            linear_solver,
            gradient_process_settings).Execute()

        # gid_output = GiDOutputProcess(model_part,
        #                            "levelset_test_2D_supg",
        #                            KratosMultiphysics.Parameters("""
        #                                {
        #                                    "result_file_configuration" : {
        #                                        "gidpost_flags": {
        #                                            "GiDPostMode": "GiD_PostBinary",
        #                                            "WriteDeformedMeshFlag": "WriteUndeformed",
        #                                            "WriteConditionsFlag": "WriteConditions",
        #                                            "MultiFileFlag": "SingleFile"
        #                                        },
        #                                        "nodal_results"       : ["DISTANCE","VELOCITY"]
        #                                    }
        #                                }
        #                                """)
        #                            )

        # gid_output.ExecuteInitialize()
        # gid_output.ExecuteBeforeSolutionLoop()
        # gid_output.ExecuteInitializeSolutionStep()
        # gid_output.PrintOutput()
        # gid_output.ExecuteFinalizeSolutionStep()
        # gid_output.ExecuteFinalize()

        for node in test_model_part.Nodes:
            grad_exact_x = 2*node.X - 4*math.cosh(4*node.X)/math.sinh(4)
            grad_exact_y = 0.0
            grad_obtained = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE_GRADIENT)
            self.assertAlmostEqual(distance_gradient[0], grad_exact_x)
            self.assertAlmostEqual(distance_gradient[1], grad_exact_y)

if __name__ == '__main__':
    KratosUnittest.main()