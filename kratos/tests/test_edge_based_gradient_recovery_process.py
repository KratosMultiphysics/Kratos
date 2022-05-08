import os
import math

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory

from KratosMultiphysics.gid_output_process import GiDOutputProcess

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestEdgeBasedGradientRecoveryProcess(KratosUnittest.TestCase):

    def _SetUpTestModelPart(self, OriginVar = None, GradientVar = None):
        # Create the background mesh
        self.model = KratosMultiphysics.Model()
        model_part = self.model.CreateModelPart("TestModelPart")
        if OriginVar:
            model_part.AddNodalSolutionStepVariable(OriginVar)
        if GradientVar:
            model_part.AddNodalSolutionStepVariable(GradientVar)

        model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 2

        problem_domain = KratosMultiphysics.Quadrilateral2D4(
            KratosMultiphysics.Node(1, 0.0, 0.0, 0.0),
            KratosMultiphysics.Node(2, 0.0, 1.0, 0.0),
            KratosMultiphysics.Node(3, 1.0, 1.0, 0.0),
            KratosMultiphysics.Node(4, 1.0, 0.0, 0.0))
        parameters = KratosMultiphysics.Parameters("{}")
        parameters.AddEmptyValue("number_of_divisions").SetInt(20)
        parameters.AddEmptyValue("element_name").SetString("Element2D3N")
        parameters.AddEmptyValue("create_skin_sub_model_part").SetBool(False)

        KratosMultiphysics.StructuredMeshGeneratorProcess(problem_domain, model_part, parameters).Execute()

    @classmethod
    def _SetUpLinearSolver(self):
        linear_solver_settings = KratosMultiphysics.Parameters("""{
            "solver_type" : "skyline_lu_factorization"
        }""")
        return linear_solver_factory.ConstructSolver(linear_solver_settings)

    def _CheckScalarResults(self, ValueGetter):
        err_norm_x = 0.0
        err_norm_y = 0.0
        test_model_part = self.model.GetModelPart("TestModelPart")
        KratosMultiphysics.CalculateNonHistoricalNodalAreaProcess(test_model_part,2).Execute()
        for node in test_model_part.Nodes:
            grad_exact_x = -25.0 * math.exp(-25*node.X)
            grad_exact_y = -25.0 * math.exp(-25*node.Y)
            nodal_area = node.GetValue(KratosMultiphysics.NODAL_AREA)
            grad_obtained = ValueGetter(node)
            err_norm_x += nodal_area * (grad_exact_x - grad_obtained[0])**2
            err_norm_y += nodal_area * (grad_exact_y - grad_obtained[1])**2
        self.assertAlmostEqual(math.sqrt(err_norm_x), 0.5935743750043295)
        self.assertAlmostEqual(math.sqrt(err_norm_y), 0.5935743752935422)

    def _CheckVectorResults(self, ValueGetter):
        err_norm = 0.0
        test_model_part = self.model.GetModelPart("TestModelPart")
        KratosMultiphysics.CalculateNonHistoricalNodalAreaProcess(test_model_part,2).Execute()
        for node in test_model_part.Nodes:
            nodal_area = node.GetValue(KratosMultiphysics.NODAL_AREA)
            grad_obtained = ValueGetter(node)
            grad_exact = [-25.0 * math.exp(-25*node.X), -25.0 * math.exp(-25*node.Y)]
            for i in range(2):
                for j in range(2):
                    err_norm += nodal_area * (grad_obtained[i,j] - grad_exact[j])**2
        self.assertAlmostEqual(math.sqrt(err_norm), 1.6713578322919727)

    def testScalarHistoricalAndHistorical(self):
        # Set up the test model part
        self._SetUpTestModelPart(KratosMultiphysics.DISTANCE, KratosMultiphysics.DISTANCE_GRADIENT)

        # Set the origin variable field
        test_model_part = self.model.GetModelPart("TestModelPart")
        for node in test_model_part.Nodes:
            value = math.exp(-25*node.X) + math.exp(-25*node.Y)
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, value)

        # Set the gradient recovery process
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
            self._SetUpLinearSolver(),
            gradient_process_settings).Execute()

        gid_output = GiDOutputProcess(
            test_model_part,
            "test_edge_based_gradient_recovery_process",
            KratosMultiphysics.Parameters("""{
                "result_file_configuration" : {
                    "gidpost_flags": {
                        "GiDPostMode": "GiD_PostBinary",
                        "WriteDeformedMeshFlag": "WriteUndeformed",
                        "WriteConditionsFlag": "WriteConditions",
                        "MultiFileFlag": "SingleFile"
                    },
                    "nodal_results" : ["DISTANCE","DISTANCE_GRADIENT"]
                }
            }""")
        )

        gid_output.ExecuteInitialize()
        gid_output.ExecuteBeforeSolutionLoop()
        gid_output.ExecuteInitializeSolutionStep()
        gid_output.PrintOutput()
        gid_output.ExecuteFinalizeSolutionStep()
        gid_output.ExecuteFinalize()

        # Check results with error norm
        self._CheckScalarResults(lambda node : node.GetSolutionStepValue(KratosMultiphysics.DISTANCE_GRADIENT))

    def testScalarNonHistoricalAndNonHistorical(self):
        # Set up the test model part
        self._SetUpTestModelPart()

        # Set the origin variable field
        test_model_part = self.model.GetModelPart("TestModelPart")
        for node in test_model_part.Nodes:
            value = math.exp(-25*node.X) + math.exp(-25*node.Y)
            node.SetValue(KratosMultiphysics.DISTANCE, value)

        # Set the gradient recovery process
        gradient_process_settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "TestModelPart",
            "origin_variable" : "DISTANCE",
            "gradient_variable" : "DISTANCE_GRADIENT",
            "gradient_penalty_coefficient" : 1.0e-6,
            "calculate_nodal_neighbours" : true,
            "is_historical_origin_variable" : false,
            "is_historical_gradient_variable": false
        }""")
        KratosMultiphysics.EdgeBasedGradientRecoveryProcessScalar(
            self.model,
            self._SetUpLinearSolver(),
            gradient_process_settings).Execute()

        # Check results with error norm
        self._CheckScalarResults(lambda node : node.GetValue(KratosMultiphysics.DISTANCE_GRADIENT))

    def testVectorHistoricalAndNonHistorical(self):
        # Set up the test model part
        self._SetUpTestModelPart(KratosMultiphysics.VELOCITY)

        # Set the origin variable field
        test_model_part = self.model.GetModelPart("TestModelPart")
        for node in test_model_part.Nodes:
            value = math.exp(-25*node.X) + math.exp(-25*node.Y)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, [value, value, 0.0])

        # Set the gradient recovery process
        gradient_process_settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "TestModelPart",
            "origin_variable" : "VELOCITY",
            "gradient_variable" : "CAUCHY_STRESS_TENSOR",
            "gradient_penalty_coefficient" : 1.0e-6,
            "calculate_nodal_neighbours" : true,
            "is_historical_origin_variable" : true,
            "is_historical_gradient_variable": false
        }""")
        KratosMultiphysics.EdgeBasedGradientRecoveryProcessArray(
            self.model,
            self._SetUpLinearSolver(),
            gradient_process_settings).Execute()

        # gid_output = GiDOutputProcess(
        #     test_model_part,
        #     "test_edge_based_gradient_recovery_process",
        #     KratosMultiphysics.Parameters("""{
        #         "result_file_configuration" : {
        #             "gidpost_flags": {
        #                 "GiDPostMode": "GiD_PostBinary",
        #                 "WriteDeformedMeshFlag": "WriteUndeformed",
        #                 "WriteConditionsFlag": "WriteConditions",
        #                 "MultiFileFlag": "SingleFile"
        #             },
        #             "nodal_results" : ["VELOCITY"],
        #             "nodal_nonhistorical_results" : ["CAUCHY_STRESS_TENSOR"]
        #         }
        #     }""")
        # )

        # gid_output.ExecuteInitialize()
        # gid_output.ExecuteBeforeSolutionLoop()
        # gid_output.ExecuteInitializeSolutionStep()
        # gid_output.PrintOutput()
        # gid_output.ExecuteFinalizeSolutionStep()
        # gid_output.ExecuteFinalize()

        # Check results with error norm
        self._CheckVectorResults(lambda node : node.GetValue(KratosMultiphysics.CAUCHY_STRESS_TENSOR))

if __name__ == '__main__':
    KratosUnittest.main()
