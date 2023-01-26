import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.gid_output_process as gid_output_process
from unittest.mock import patch


class TestGidOutputProcess(KratosUnittest.TestCase):

    def setUp(self):
        self.model = KM.Model()
        self.model_part = self.model.CreateModelPart('test_model_part')
        self.model_part.ProcessInfo[KM.STEP] = 0
        self.model_part.ProcessInfo[KM.TIME] = 0
        self.patcher1 = patch('KratosMultiphysics.GidIO.InitializeMesh', autospec=True)
        self.patcher2 = patch('KratosMultiphysics.GidIO.WriteMesh', autospec=True)
        self.patcher3 = patch('KratosMultiphysics.GidIO.WriteNodeMesh', autospec=True)
        self.patcher4 = patch('KratosMultiphysics.GidIO.FinalizeMesh', autospec=True)
        self.patcher5 = patch('KratosMultiphysics.GidIO.InitializeResults', autospec=True)
        self.patcher6 = patch('KratosMultiphysics.GidIO.WriteNodalResults', autospec=True)
        self.patcher7 = patch('KratosMultiphysics.GidIO.PrintOnGaussPoints', autospec=True)
        self.patcher8 = patch('KratosMultiphysics.GidIO.WriteNodalResultsNonHistorical', autospec=True)
        self.patcher9 = patch('KratosMultiphysics.GidIO.WriteNodalFlags', autospec=True)
        self.patcher10 = patch('KratosMultiphysics.GidIO.PrintFlagsOnGaussPoints', autospec=True)
        self.patcher11 = patch('KratosMultiphysics.GidIO.FinalizeResults', autospec=True)
        self.patcher12 = patch('KratosMultiphysics.CuttingUtility.FindSmallestEdge', autospec=True)
        self.patcher13 = patch('KratosMultiphysics.CuttingUtility.AddVariablesToCutModelPart', autospec=True)
        self.patcher14 = patch('KratosMultiphysics.CuttingUtility.AddSkinConditions', autospec=True)
        self.patcher15 = patch('KratosMultiphysics.CuttingUtility.GenerateCut', autospec=True)
        self.patcher16 = patch('KratosMultiphysics.CuttingUtility.UpdateCutData', autospec=True)
        self.GidIOInitializeMesh = self.patcher1.start()
        self.GidIOWriteMesh = self.patcher2.start()
        self.GidIOWriteNodeMesh = self.patcher3.start()
        self.GidIOFinalizeMesh = self.patcher4.start()
        self.GidIOInitializeResults = self.patcher5.start()
        self.GidIOWriteNodalResults = self.patcher6.start()
        self.GidIOPrintOnGaussPoints = self.patcher7.start()
        self.GidIOWriteNodalResultsNonHistorical = self.patcher8.start()
        self.GidIOWriteNodalFlags = self.patcher9.start()
        self.GidIOPrintFlagsOnGaussPoints = self.patcher10.start()
        self.GidIOFinalizeResults = self.patcher11.start()
        self.CuttingUtilityFindSmallestEdge = self.patcher12.start()
        self.CuttingUtilityAddVariablesToCutModelPart = self.patcher13.start()
        self.CuttingUtilityAddSkinConditions = self.patcher14.start()
        self.CuttingUtilityGenerateCut = self.patcher15.start()
        self.CuttingUtilityUpdateCutData = self.patcher16.start()


    def tearDown(self):
        self.patcher1.stop()
        self.patcher2.stop()
        self.patcher3.stop()
        self.patcher4.stop()
        self.patcher5.stop()
        self.patcher6.stop()
        self.patcher7.stop()
        self.patcher8.stop()
        self.patcher9.stop()
        self.patcher10.stop()
        self.patcher11.stop()
        self.patcher12.stop()
        self.patcher13.stop()
        self.patcher14.stop()
        self.patcher15.stop()
        self.patcher16.stop()


    def testGidOutputBodyIO(self):
        settings = KM.Parameters("""
        {
            "Parameters": {
                "model_part_name": "test_model_part",
                "output_name": "output_file",
                "postprocess_parameters" : {
                    "result_file_configuration": {
                        "gidpost_flags": {
                            "GiDPostMode": "GiD_PostBinary",
                            "WriteDeformedMeshFlag": "WriteUndeformed",
                            "WriteConditionsFlag": "WriteElementsOnly",
                            "MultiFileFlag": "SingleFile"
                        },
                        "file_label": "time",
                        "time_label_format": "{:.12f}",
                        "output_control_type": "step",
                        "output_interval": 1.0,
                        "flush_after_output": false,
                        "body_output": true,
                        "node_output": false,
                        "skin_output": false,
                        "plane_output": [],
                        "nodal_results": ["DISPLACEMENT"],
                        "nodal_nonhistorical_results": ["NODAL_AREA"],
                        "nodal_flags_results": ["STRUCTURE"],
                        "elemental_conditional_flags_results": ["FLUID"],
                        "gauss_point_results": ["STRESSES"],
                        "additional_list_files": []
                    }
                }
            }
        }""")
        process = gid_output_process.Factory(settings, self.model)
        process.ExecuteInitialize()
        process.ExecuteBeforeSolutionLoop()
        process.PrintOutput()
        self.assertEqual(self.GidIOInitializeMesh.call_count, 1)
        self.assertEqual(self.GidIOWriteMesh.call_count, 1)
        self.assertEqual(self.GidIOWriteNodeMesh.call_count, 0)
        self.assertEqual(self.GidIOFinalizeMesh.call_count, 1)
        self.assertEqual(self.GidIOInitializeResults.call_count, 1)
        self.assertEqual(self.GidIOWriteNodalResults.call_count, 1)
        self.assertEqual(self.GidIOPrintOnGaussPoints.call_count, 1)
        self.assertEqual(self.GidIOWriteNodalResultsNonHistorical.call_count, 1)
        self.assertEqual(self.GidIOWriteNodalFlags.call_count, 1)
        self.assertEqual(self.GidIOPrintFlagsOnGaussPoints.call_count, 1)
        process.ExecuteFinalize()
        self.assertEqual(self.GidIOFinalizeResults.call_count, 1)
        self.assertEqual(self.CuttingUtilityGenerateCut.call_count, 0)


    def testGidOutputNodeIO(self):
        settings = KM.Parameters("""
        {
            "Parameters": {
                "model_part_name": "test_model_part",
                "output_name": "output_file",
                "postprocess_parameters" : {
                    "result_file_configuration": {
                        "gidpost_flags": {
                            "GiDPostMode": "GiD_PostBinary",
                            "WriteDeformedMeshFlag": "WriteUndeformed",
                            "WriteConditionsFlag": "WriteElementsOnly",
                            "MultiFileFlag": "SingleFile"
                        },
                        "file_label": "time",
                        "time_label_format": "{:.12f}",
                        "output_control_type": "step",
                        "output_interval": 1.0,
                        "flush_after_output": false,
                        "body_output": false,
                        "node_output": true,
                        "skin_output": false,
                        "plane_output": [],
                        "nodal_results": ["DISPLACEMENT"],
                        "nodal_nonhistorical_results": ["NODAL_AREA"],
                        "nodal_flags_results": ["STRUCTURE"],
                        "elemental_conditional_flags_results": ["FLUID"],
                        "gauss_point_results": ["STRESSES"],
                        "additional_list_files": []
                    }
                }
            }
        }""")
        process = gid_output_process.Factory(settings, self.model)
        process.ExecuteInitialize()
        process.ExecuteBeforeSolutionLoop()
        process.PrintOutput()
        self.assertEqual(self.GidIOInitializeMesh.call_count, 1)
        self.assertEqual(self.GidIOWriteMesh.call_count, 0)
        self.assertEqual(self.GidIOWriteNodeMesh.call_count, 1)
        self.assertEqual(self.GidIOFinalizeMesh.call_count, 1)
        self.assertEqual(self.GidIOInitializeResults.call_count, 1)
        self.assertEqual(self.GidIOWriteNodalResults.call_count, 1)
        self.assertEqual(self.GidIOPrintOnGaussPoints.call_count, 0)
        self.assertEqual(self.GidIOWriteNodalResultsNonHistorical.call_count, 1)
        self.assertEqual(self.GidIOWriteNodalFlags.call_count, 1)
        self.assertEqual(self.GidIOPrintFlagsOnGaussPoints.call_count, 1)
        process.ExecuteFinalize()
        self.assertEqual(self.GidIOFinalizeResults.call_count, 1)
        self.assertEqual(self.CuttingUtilityGenerateCut.call_count, 0)


    def testGidOutputSkinIO(self):
        settings = KM.Parameters("""
        {
            "Parameters": {
                "model_part_name": "test_model_part",
                "output_name": "output_file",
                "postprocess_parameters" : {
                    "result_file_configuration": {
                        "gidpost_flags": {
                            "GiDPostMode": "GiD_PostBinary",
                            "WriteDeformedMeshFlag": "WriteUndeformed",
                            "WriteConditionsFlag": "WriteElementsOnly",
                            "MultiFileFlag": "SingleFile"
                        },
                        "file_label": "time",
                        "time_label_format": "{:.12f}",
                        "output_control_type": "step",
                        "output_interval": 1.0,
                        "flush_after_output": false,
                        "body_output": false,
                        "node_output": false,
                        "skin_output": true,
                        "plane_output": [],
                        "nodal_results": ["DISPLACEMENT"],
                        "nodal_nonhistorical_results": ["NODAL_AREA"],
                        "nodal_flags_results": ["STRUCTURE"],
                        "elemental_conditional_flags_results": ["FLUID"],
                        "gauss_point_results": ["STRESSES"],
                        "additional_list_files": []
                    }
                }
            }
        }""")
        process = gid_output_process.Factory(settings, self.model)
        process.ExecuteInitialize()
        process.ExecuteBeforeSolutionLoop()
        process.PrintOutput()
        self.assertEqual(self.GidIOInitializeMesh.call_count, 1)
        self.assertEqual(self.GidIOWriteMesh.call_count, 1)
        self.assertEqual(self.GidIOWriteNodeMesh.call_count, 0)
        self.assertEqual(self.GidIOFinalizeMesh.call_count, 1)
        self.assertEqual(self.GidIOInitializeResults.call_count, 1)
        self.assertEqual(self.GidIOWriteNodalResults.call_count, 1)
        self.assertEqual(self.GidIOPrintOnGaussPoints.call_count, 0)
        self.assertEqual(self.GidIOWriteNodalResultsNonHistorical.call_count, 1)
        self.assertEqual(self.GidIOWriteNodalFlags.call_count, 1)
        self.assertEqual(self.GidIOPrintFlagsOnGaussPoints.call_count, 1)
        process.ExecuteFinalize()
        self.assertEqual(self.GidIOFinalizeResults.call_count, 1)
        self.assertEqual(self.CuttingUtilityGenerateCut.call_count, 0)


    def testGidOutputPlaneIO(self):
        settings = KM.Parameters("""
        {
            "Parameters": {
                "model_part_name": "test_model_part",
                "output_name": "output_file",
                "postprocess_parameters" : {
                    "result_file_configuration": {
                        "gidpost_flags": {
                            "GiDPostMode": "GiD_PostBinary",
                            "WriteDeformedMeshFlag": "WriteUndeformed",
                            "WriteConditionsFlag": "WriteElementsOnly",
                            "MultiFileFlag": "SingleFile"
                        },
                        "file_label": "time",
                        "time_label_format": "{:.12f}",
                        "output_control_type": "step",
                        "output_interval": 1.0,
                        "flush_after_output": false,
                        "body_output": false,
                        "node_output": false,
                        "skin_output": false,
                        "plane_output": [{
                            "normal": [1.0, 0.0, 0.0],
                            "point" : [0.0, 0.0, 0.0]
                        }],
                        "nodal_results": ["DISPLACEMENT"],
                        "nodal_nonhistorical_results": ["NODAL_AREA"],
                        "nodal_flags_results": ["STRUCTURE"],
                        "elemental_conditional_flags_results": ["FLUID"],
                        "gauss_point_results": ["STRESSES"],
                        "additional_list_files": []
                    }
                }
            }
        }""")
        process = gid_output_process.Factory(settings, self.model)
        process.ExecuteInitialize()
        self.assertEqual(self.CuttingUtilityFindSmallestEdge.call_count, 1)
        self.assertEqual(self.CuttingUtilityAddVariablesToCutModelPart.call_count, 1)
        self.assertEqual(self.CuttingUtilityAddSkinConditions.call_count, 0)
        self.assertEqual(self.CuttingUtilityGenerateCut.call_count, 1)
        process.ExecuteBeforeSolutionLoop()
        process.PrintOutput()
        self.assertEqual(self.GidIOInitializeMesh.call_count, 1)
        self.assertEqual(self.GidIOWriteMesh.call_count, 1)
        self.assertEqual(self.GidIOWriteNodeMesh.call_count, 0)
        self.assertEqual(self.GidIOFinalizeMesh.call_count, 1)
        self.assertEqual(self.GidIOInitializeResults.call_count, 1)
        self.assertEqual(self.GidIOWriteNodalResults.call_count, 1)
        self.assertEqual(self.GidIOPrintOnGaussPoints.call_count, 0)
        self.assertEqual(self.GidIOWriteNodalResultsNonHistorical.call_count, 1)
        self.assertEqual(self.GidIOWriteNodalFlags.call_count, 1)
        self.assertEqual(self.GidIOPrintFlagsOnGaussPoints.call_count, 1)
        self.assertEqual(self.CuttingUtilityUpdateCutData.call_count, 4)
        process.ExecuteFinalize()
        self.assertEqual(self.GidIOFinalizeResults.call_count, 1)


if __name__ == "__main__":
    KratosUnittest.main()
