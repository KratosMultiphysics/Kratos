import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.HDF5Application.single_mesh_temporal_output_process as single_mesh_temporal_output_process
import KratosMultiphysics.HDF5Application.multiple_mesh_temporal_output_process as multiple_mesh_temporal_output_process
import KratosMultiphysics.HDF5Application.single_mesh_primal_output_process as single_mesh_primal_output_process
import KratosMultiphysics.HDF5Application.initialization_from_hdf5_process as initialization_from_hdf5_process
import KratosMultiphysics.HDF5Application.single_mesh_temporal_input_process as single_mesh_temporal_input_process
import KratosMultiphysics.HDF5Application.single_mesh_xdmf_output_process as single_mesh_xdmf_output_process
from unittest.mock import patch


class TestHDF5Processes(KratosUnittest.TestCase):

    def setUp(self):
        self.model = KratosMultiphysics.Model()
        self.model_part = self.model.CreateModelPart('test_model_part')
        self.model_part.ProcessInfo[KratosMultiphysics.TIME] = 0.0
        self.model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME] = 0.1
        self.patcher1 = patch(
            'KratosMultiphysics.HDF5Application.core.file_io.KratosHDF5.HDF5FileSerial', autospec=True)
        self.patcher2 = patch(
            'KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5ModelPartIO', autospec=True)
        self.patcher3 = patch(
            'KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5NodalSolutionStepDataIO', autospec=True)
        self.patcher4 = patch(
            'KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5NodalDataValueIO', autospec=True)
        self.patcher5 = patch(
            'KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5ElementDataValueIO', autospec=True)
        self.patcher6 = patch(
            'KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5NodalSolutionStepBossakIO', autospec=True)
        self.patcher7 = patch(
            'KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5ElementFlagValueIO', autospec=True)
        self.patcher8 = patch(
            'KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5NodalFlagValueIO', autospec=True)
        self.patcher9 = patch(
            'KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5ConditionFlagValueIO', autospec=True)
        self.patcher10 = patch(
            'KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5ConditionDataValueIO', autospec=True)
        self.HDF5FileSerial = self.patcher1.start()
        self.HDF5ModelPartIO = self.patcher2.start()
        self.HDF5NodalSolutionStepDataIO = self.patcher3.start()
        self.HDF5NodalDataValueIO = self.patcher4.start()
        self.HDF5ElementDataValueIO = self.patcher5.start()
        self.HDF5NodalSolutionStepBossakIO = self.patcher6.start()
        self.HDF5ElementFlagValueIO = self.patcher7.start()
        self.HDF5NodalFlagValueIO = self.patcher8.start()
        self.HDF5ConditionFlagValueIO = self.patcher9.start()
        self.HDF5ConditionDataValueIO = self.patcher10.start()

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

    def test_SingleMeshTemporalOutputProcess(self):
        settings = KratosMultiphysics.Parameters('''
            {
                "Parameters": {
                    "model_part_name": "test_model_part",
                    "file_settings": {
                        "file_access_mode": "truncate",
                        "echo_level": 1
                    },
                    "model_part_output_settings": {
                        "prefix": "/ModelData/<model_part_name>"
                    },
                    "nodal_solution_step_data_settings": {
                        "list_of_variables": ["DISPLACEMENT"],
                        "prefix": "/ResultsData/<model_part_name>/<time>",
                        "time_format": "0.2f"
                    },
                    "element_data_value_settings": {
                        "prefix": "/ResultsData/ElementDataValues"
                    },
                    "nodal_flag_value_settings": {
                        "prefix": "/ResultsData/NodalFlagValues"
                    },
                    "element_flag_value_settings": {
                        "prefix": "/ResultsData/ElementFlagValues"
                    },
                    "condition_data_value_settings": {
                        "prefix": "/ResultsData/ConditionDataValues"
                    },
                    "condition_flag_value_settings": {
                        "prefix": "/ResultsData/ConditionFlagValues"
                    },
                    "output_time_settings": {
                        "time_frequency": 0.2,
                        "step_frequency": 10
                    }
                }
            }
            ''')
        process = single_mesh_temporal_output_process.Factory(
            settings, self.model)
        process.ExecuteBeforeSolutionLoop()
        self.assertEqual(self.HDF5FileSerial.call_count, 1)
        self.assertEqual(
            self.HDF5FileSerial.call_args[0][0]['file_name'].GetString(), 'test_model_part.h5')
        for time in [0.09999999, 0.19999998]:
            self.model_part.CloneTimeStep(time)
            process.ExecuteFinalizeSolutionStep()
        self.assertEqual(self.HDF5FileSerial.call_count, 2)
        self.assertEqual(self.HDF5FileSerial.call_args[0][0]['file_name'].GetString(
        ), 'test_model_part-0.2000.h5')
        self.assertEqual(
            self.HDF5FileSerial.call_args[0][0]['file_access_mode'].GetString(), 'truncate')
        self.assertEqual(
            self.HDF5FileSerial.call_args[0][0]['echo_level'].GetInt(), 1)
        self.HDF5ModelPartIO.assert_called_once_with(
            self.HDF5FileSerial.return_value, '/ModelData/test_model_part')
        self.HDF5ModelPartIO.return_value.WriteModelPart.assert_called_once_with(
            self.model_part)
        self.assertEqual(self.HDF5NodalSolutionStepDataIO.call_count, 2)
        self.assertEqual(self.HDF5NodalSolutionStepDataIO.call_args[0][0]['prefix'].GetString(
        ), '/ResultsData/test_model_part/0.20')
        self.assertEqual(
            self.HDF5NodalSolutionStepDataIO.call_args[0][0]['list_of_variables'][0].GetString(), 'DISPLACEMENT')
        self.assertEqual(
            self.HDF5NodalSolutionStepDataIO.return_value.WriteNodalResults.call_count, 2)
        self.HDF5NodalSolutionStepDataIO.return_value.WriteNodalResults.assert_called_with(
            self.model_part.Nodes, 0)
        self.assertEqual(self.HDF5NodalDataValueIO.call_count, 2)
        self.assertEqual(
            self.HDF5NodalDataValueIO.call_args[0][0]['prefix'].GetString(), '/ResultsData')
        self.assertEqual(
            self.HDF5NodalDataValueIO.call_args[0][0]['list_of_variables'].size(), 0)
        self.assertEqual(
            self.HDF5NodalDataValueIO.return_value.WriteNodalResults.call_count, 2)
        self.HDF5NodalDataValueIO.return_value.WriteNodalResults.assert_called_with(
            self.model_part.Nodes)

        self.assertEqual(self.HDF5NodalFlagValueIO.call_count, 2)
        self.assertEqual(
            self.HDF5NodalFlagValueIO.call_args[0][0]['prefix'].GetString(), '/ResultsData/NodalFlagValues')
        self.assertEqual(
            self.HDF5NodalFlagValueIO.call_args[0][0]['list_of_variables'].size(), 0)
        self.assertEqual(
            self.HDF5NodalFlagValueIO.return_value.WriteNodalFlags.call_count, 2)
        self.HDF5NodalFlagValueIO.return_value.WriteNodalFlags.assert_called_with(
            self.model_part.Nodes)

        self.assertEqual(self.HDF5ElementDataValueIO.call_count, 2)
        self.assertEqual(self.HDF5ElementDataValueIO.call_args[0][0]['prefix'].GetString(
        ), '/ResultsData/ElementDataValues')
        self.assertEqual(
            self.HDF5ElementDataValueIO.call_args[0][0]['list_of_variables'].size(), 0)
        self.assertEqual(
            self.HDF5ElementDataValueIO.return_value.WriteElementResults.call_count, 2)
        self.HDF5ElementDataValueIO.return_value.WriteElementResults.assert_called_with(
            self.model_part.Elements)

        self.assertEqual(self.HDF5ElementFlagValueIO.call_count, 2)
        self.assertEqual(self.HDF5ElementFlagValueIO.call_args[0][0]['prefix'].GetString(
        ), '/ResultsData/ElementFlagValues')
        self.assertEqual(
            self.HDF5ElementFlagValueIO.call_args[0][0]['list_of_variables'].size(), 0)
        self.assertEqual(
            self.HDF5ElementFlagValueIO.return_value.WriteElementFlags.call_count, 2)
        self.HDF5ElementFlagValueIO.return_value.WriteElementFlags.assert_called_with(
            self.model_part.Elements)

        self.assertEqual(self.HDF5ConditionDataValueIO.call_count, 2)
        self.assertEqual(self.HDF5ConditionDataValueIO.call_args[0][0]['prefix'].GetString(
        ), '/ResultsData/ConditionDataValues')
        self.assertEqual(
            self.HDF5ConditionDataValueIO.call_args[0][0]['list_of_variables'].size(), 0)
        self.assertEqual(
            self.HDF5ConditionDataValueIO.return_value.WriteConditionResults.call_count, 2)
        self.HDF5ConditionDataValueIO.return_value.WriteConditionResults.assert_called_with(
            self.model_part.Conditions)

        self.assertEqual(self.HDF5ConditionFlagValueIO.call_count, 2)
        self.assertEqual(self.HDF5ConditionFlagValueIO.call_args[0][0]['prefix'].GetString(
        ), '/ResultsData/ConditionFlagValues')
        self.assertEqual(
            self.HDF5ConditionFlagValueIO.call_args[0][0]['list_of_variables'].size(), 0)
        self.assertEqual(
            self.HDF5ConditionFlagValueIO.return_value.WriteConditionFlags.call_count, 2)
        self.HDF5ConditionFlagValueIO.return_value.WriteConditionFlags.assert_called_with(
            self.model_part.Conditions)

    def test_MultipleMeshTemporalOutputProcess(self):
        settings = KratosMultiphysics.Parameters('''
            {
                "Parameters": {
                    "model_part_name": "test_model_part",
                    "file_settings": {
                        "file_name": "kratos-<time>.h5"
                    },
                    "output_time_settings": {
                        "time_frequency": 0.2,
                        "step_frequency": 1
                    }
                }
            }
            ''')
        process = multiple_mesh_temporal_output_process.Factory(
            settings, self.model)
        process.ExecuteBeforeSolutionLoop()
        self.assertEqual(self.HDF5FileSerial.call_count, 1)
        self.assertEqual(
            self.HDF5FileSerial.call_args[0][0]['file_name'].GetString(), 'kratos-0.0000.h5')
        for time in [0.09999999, 0.19999998]:
            self.model_part.CloneTimeStep(time)
            process.ExecuteFinalizeSolutionStep()
        self.assertEqual(self.HDF5FileSerial.call_count, 3)
        self.assertEqual(
            self.HDF5FileSerial.call_args[0][0]['file_name'].GetString(), 'kratos-0.2000.h5')
        self.assertEqual(
            self.HDF5FileSerial.call_args[0][0]['file_access_mode'].GetString(), 'exclusive')
        self.assertEqual(
            self.HDF5FileSerial.call_args[0][0]['echo_level'].GetInt(), 0)
        self.assertEqual(self.HDF5ModelPartIO.call_count, 3)
        self.HDF5ModelPartIO.assert_called_with(
            self.HDF5FileSerial.return_value, '/ModelData')
        self.assertEqual(
            self.HDF5ModelPartIO.return_value.WriteModelPart.call_count, 3)
        self.HDF5ModelPartIO.return_value.WriteModelPart.assert_called_with(
            self.model_part)

    def test_SingleMeshPrimalOutputProcess(self):
        settings = KratosMultiphysics.Parameters('''
            {
                "Parameters": {
                    "model_part_name": "test_model_part",
                    "nodal_solution_step_data_settings": {
                        "alpha_bossak": -0.2
                    },
                    "output_time_settings": {
                        "step_frequency": 10,
                        "time_frequency": 0.1
                    }
                }
            }
            ''')
        process = single_mesh_primal_output_process.Factory(
            settings, self.model)
        process.ExecuteBeforeSolutionLoop()
        for time in [0.09999999, 0.19999998]:
            self.model_part.CloneTimeStep(time)
            process.ExecuteFinalizeSolutionStep()
        self.assertEqual(self.HDF5NodalSolutionStepBossakIO.call_count, 3)
        self.assertEqual(
            self.HDF5NodalSolutionStepBossakIO.call_args[0][0]['prefix'].GetString(), '/ResultsData')
        self.assertEqual(
            self.HDF5NodalSolutionStepBossakIO.call_args[0][0]['list_of_variables'].size(), 0)
        self.assertEqual(
            self.HDF5NodalSolutionStepBossakIO.call_args[0][1], self.HDF5FileSerial.return_value)
        self.assertEqual(
            self.HDF5NodalSolutionStepBossakIO.return_value.SetAlphaBossak.call_count, 3)
        self.HDF5NodalSolutionStepBossakIO.return_value.SetAlphaBossak.assert_called_with(
            -0.2)
        self.assertEqual(
            self.HDF5NodalSolutionStepBossakIO.return_value.WriteNodalResults.call_count, 3)
        self.HDF5NodalSolutionStepBossakIO.return_value.WriteNodalResults.assert_called_with(
            self.model_part.Nodes)

    def test_InitializationFromHDF5Process(self):
        settings = KratosMultiphysics.Parameters('''
            {
                "Parameters": {
                    "model_part_name": "test_model_part",
                    "file_settings": {
                        "file_name": "kratos"
                    }
                }
            }
            ''')
        process = initialization_from_hdf5_process.Factory(
            settings, self.model)
        process.ExecuteInitialize()
        self.assertEqual(self.HDF5FileSerial.call_count, 1)
        self.assertEqual(
            self.HDF5FileSerial.call_args[0][0]["file_name"].GetString(), "kratos.h5")
        self.assertEqual(
            self.HDF5FileSerial.call_args[0][0]["file_access_mode"].GetString(), "read_only")
        self.assertEqual(self.HDF5NodalSolutionStepDataIO.call_count, 1)
        self.assertEqual(
            self.HDF5NodalSolutionStepDataIO.return_value.ReadNodalResults.call_count, 1)
        self.HDF5NodalSolutionStepDataIO.return_value.ReadNodalResults.assert_called_with(
            self.model_part.Nodes, self.model_part.GetCommunicator(), 0)
        self.assertEqual(
            self.HDF5NodalDataValueIO.return_value.ReadNodalResults.call_count, 1)
        self.HDF5NodalDataValueIO.return_value.ReadNodalResults.assert_called_with(
            self.model_part.Nodes, self.model_part.GetCommunicator())
        self.assertEqual(
            self.HDF5NodalFlagValueIO.return_value.ReadNodalFlags.call_count, 1)
        self.HDF5NodalFlagValueIO.return_value.ReadNodalFlags.assert_called_with(
            self.model_part.Nodes, self.model_part.GetCommunicator())
        self.assertEqual(
            self.HDF5ElementDataValueIO.return_value.ReadElementResults.call_count, 1)
        self.HDF5ElementDataValueIO.return_value.ReadElementResults.assert_called_with(
            self.model_part.Elements, self.model_part.GetCommunicator())
        self.assertEqual(
            self.HDF5ElementFlagValueIO.return_value.ReadElementFlags.call_count, 1)
        self.HDF5ElementFlagValueIO.return_value.ReadElementFlags.assert_called_with(
            self.model_part.Elements, self.model_part.GetCommunicator())
        self.assertEqual(
            self.HDF5ConditionDataValueIO.return_value.ReadConditionResults.call_count, 1)
        self.HDF5ConditionDataValueIO.return_value.ReadConditionResults.assert_called_with(
            self.model_part.Conditions, self.model_part.GetCommunicator())
        self.assertEqual(
            self.HDF5ConditionFlagValueIO.return_value.ReadConditionFlags.call_count, 1)
        self.HDF5ConditionFlagValueIO.return_value.ReadConditionFlags.assert_called_with(
            self.model_part.Conditions, self.model_part.GetCommunicator())

    def test_SingleMeshTemporalInputProcess(self):
        settings = KratosMultiphysics.Parameters('''
            {
                "Parameters": {
                    "model_part_name": "test_model_part"
                }
            }
            ''')
        process = single_mesh_temporal_input_process.Factory(
            settings, self.model)
        for time in [0.1, 0.2]:
            self.model_part.CloneTimeStep(time)
            process.ExecuteInitializeSolutionStep()
        self.assertEqual(self.HDF5FileSerial.call_count, 2)
        self.assertEqual(
            self.HDF5FileSerial.call_args[0][0]["file_name"].GetString(), "test_model_part-0.2000.h5")
        self.assertEqual(
            self.HDF5FileSerial.call_args[0][0]["file_access_mode"].GetString(), "read_only")
        self.assertEqual(self.HDF5NodalSolutionStepDataIO.call_count, 2)
        self.assertEqual(
            self.HDF5NodalSolutionStepDataIO.return_value.ReadNodalResults.call_count, 2)
        self.HDF5NodalSolutionStepDataIO.return_value.ReadNodalResults.assert_called_with(
            self.model_part.Nodes, self.model_part.GetCommunicator(), 0)
        self.assertEqual(
            self.HDF5NodalDataValueIO.return_value.ReadNodalResults.call_count, 2)
        self.HDF5NodalDataValueIO.return_value.ReadNodalResults.assert_called_with(
            self.model_part.Nodes, self.model_part.GetCommunicator())
        self.assertEqual(
            self.HDF5ElementDataValueIO.return_value.ReadElementResults.call_count, 2)
        self.HDF5ElementDataValueIO.return_value.ReadElementResults.assert_called_with(
            self.model_part.Elements, self.model_part.GetCommunicator())
        self.assertEqual(
            self.HDF5NodalFlagValueIO.return_value.ReadNodalFlags.call_count, 2)
        self.HDF5NodalFlagValueIO.return_value.ReadNodalFlags.assert_called_with(
            self.model_part.Nodes, self.model_part.GetCommunicator())
        self.assertEqual(
            self.HDF5ElementFlagValueIO.return_value.ReadElementFlags.call_count, 2)
        self.HDF5ElementFlagValueIO.return_value.ReadElementFlags.assert_called_with(
            self.model_part.Elements, self.model_part.GetCommunicator())
        self.assertEqual(
            self.HDF5ConditionDataValueIO.return_value.ReadConditionResults.call_count, 2)
        self.HDF5ConditionDataValueIO.return_value.ReadConditionResults.assert_called_with(
            self.model_part.Conditions, self.model_part.GetCommunicator())
        self.assertEqual(
            self.HDF5NodalFlagValueIO.return_value.ReadNodalFlags.call_count, 2)
        self.HDF5NodalFlagValueIO.return_value.ReadNodalFlags.assert_called_with(
            self.model_part.Nodes, self.model_part.GetCommunicator())
        self.assertEqual(
            self.HDF5ConditionFlagValueIO.return_value.ReadConditionFlags.call_count, 2)
        self.HDF5ConditionFlagValueIO.return_value.ReadConditionFlags.assert_called_with(
            self.model_part.Conditions, self.model_part.GetCommunicator())

    def test_SingleMeshXdmfOutputProcess(self):
        settings = KratosMultiphysics.Parameters('''
            {
                "Parameters": {
                    "model_part_name": "test_model_part",
                    "output_time_settings": {
                        "step_frequency": 10,
                        "time_frequency": 0.1
                    }
                }
            }
            ''')
        patcher1 = patch(
            'KratosMultiphysics.HDF5Application.xdmf_utils.WriteMultifileTemporalAnalysisToXdmf', autospec=True)
        patcher2 = patch(
            'KratosMultiphysics.kratos_utilities.DeleteFileIfExisting', autospec=True)
        patcher3 = patch('os.listdir', autospec=True)
        WriteMultifileTemporalAnalysisToXdmf = patcher1.start()
        DeleteFileIfExisting = patcher2.start()
        listdir = patcher3.start()
        listdir.return_value = [
            'test_model_part-0.0000.h5', 'test_model_part-0.1000.h5']
        process = single_mesh_xdmf_output_process.Factory(settings, self.model)
        process.ExecuteInitialize()
        process.ExecuteBeforeSolutionLoop()
        for time in [0.09999999, 0.19999998]:
            self.model_part.CloneTimeStep(time)
            process.ExecuteFinalizeSolutionStep()
        self.assertEqual(WriteMultifileTemporalAnalysisToXdmf.call_count, 2)
        WriteMultifileTemporalAnalysisToXdmf.assert_called_with(
            'test_model_part.h5', '/ModelData', '/ResultsData')
        DeleteFileIfExisting.assert_called_once_with(
            './test_model_part-0.1000.h5')
        patcher1.stop()
        patcher2.stop()
        patcher3.stop()


if __name__ == "__main__":
    KratosUnittest.main()
