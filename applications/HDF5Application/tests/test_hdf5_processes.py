# --- Kratos Imports ---
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities
import KratosMultiphysics.HDF5Application.single_mesh_temporal_output_process as single_mesh_temporal_output_process
import KratosMultiphysics.HDF5Application.multiple_mesh_temporal_output_process as multiple_mesh_temporal_output_process
import KratosMultiphysics.HDF5Application.single_mesh_primal_output_process as single_mesh_primal_output_process
import KratosMultiphysics.HDF5Application.initialization_from_hdf5_process as initialization_from_hdf5_process
import KratosMultiphysics.HDF5Application.single_mesh_temporal_input_process as single_mesh_temporal_input_process
import KratosMultiphysics.HDF5Application.single_mesh_xdmf_output_process as single_mesh_xdmf_output_process
import KratosMultiphysics.HDF5Application.import_model_part_from_hdf5_process as import_model_part_from_hdf5_process

# --- STD Imports ---
from unittest.mock import patch
import pathlib


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
        self.patcher11 = patch(
            'KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5ElementGaussPointOutput', autospec=True)
        self.patcher12 = patch(
            'KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5ConditionGaussPointOutput', autospec=True)
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
        self.HDF5ElementGaussPointOutput = self.patcher11.start()
        self.HDF5ConditionGaussPointOutput = self.patcher12.start()

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
                    "element_gauss_point_value_settings"      : {
                        "prefix": "/ResultsData/ElementGaussPointValues"
                    },
                    "condition_data_value_settings": {
                        "prefix": "/ResultsData/ConditionDataValues"
                    },
                    "condition_flag_value_settings": {
                        "prefix": "/ResultsData/ConditionFlagValues"
                    },
                    "condition_gauss_point_value_settings"      : {
                        "prefix": "/ResultsData/ConditionGaussPointValues"
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
            self.model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
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
        self.assertEqual(self.HDF5NodalSolutionStepDataIO.call_count, 1)
        self.assertEqual(self.HDF5NodalSolutionStepDataIO.call_args[0][0]['prefix'].GetString(), '/ResultsData/test_model_part/0.20')
        self.assertEqual(
            self.HDF5NodalSolutionStepDataIO.call_args[0][0]['list_of_variables'][0].GetString(), 'DISPLACEMENT')
        self.assertEqual(
            self.HDF5NodalSolutionStepDataIO.return_value.WriteNodalResults.call_count, 1)
        self.HDF5NodalSolutionStepDataIO.return_value.WriteNodalResults.assert_called_with(
            self.model_part, 0)
        self.assertEqual(self.HDF5NodalDataValueIO.call_count, 1)
        self.assertEqual(
            self.HDF5NodalDataValueIO.call_args[0][0]['prefix'].GetString(), '/ResultsData')
        self.assertEqual(
            self.HDF5NodalDataValueIO.call_args[0][0]['list_of_variables'].size(), 0)
        self.assertEqual(
            self.HDF5NodalDataValueIO.return_value.WriteNodalResults.call_count, 1)
        self.HDF5NodalDataValueIO.return_value.WriteNodalResults.assert_called_with(
            self.model_part.Nodes)

        self.assertEqual(self.HDF5NodalFlagValueIO.call_count, 1)
        self.assertEqual(
            self.HDF5NodalFlagValueIO.call_args[0][0]['prefix'].GetString(), '/ResultsData/NodalFlagValues')
        self.assertEqual(
            self.HDF5NodalFlagValueIO.call_args[0][0]['list_of_variables'].size(), 0)
        self.assertEqual(
            self.HDF5NodalFlagValueIO.return_value.WriteNodalFlags.call_count, 1)
        self.HDF5NodalFlagValueIO.return_value.WriteNodalFlags.assert_called_with(
            self.model_part.Nodes)

        self.assertEqual(self.HDF5ElementDataValueIO.call_count, 1)
        self.assertEqual(self.HDF5ElementDataValueIO.call_args[0][0]['prefix'].GetString(
        ), '/ResultsData/ElementDataValues')
        self.assertEqual(
            self.HDF5ElementDataValueIO.call_args[0][0]['list_of_variables'].size(), 0)
        self.assertEqual(
            self.HDF5ElementDataValueIO.return_value.WriteElementResults.call_count, 1)
        self.HDF5ElementDataValueIO.return_value.WriteElementResults.assert_called_with(
            self.model_part.Elements)

        self.assertEqual(self.HDF5ElementFlagValueIO.call_count, 1)
        self.assertEqual(self.HDF5ElementFlagValueIO.call_args[0][0]['prefix'].GetString(
        ), '/ResultsData/ElementFlagValues')
        self.assertEqual(
            self.HDF5ElementFlagValueIO.call_args[0][0]['list_of_variables'].size(), 0)
        self.assertEqual(
            self.HDF5ElementFlagValueIO.return_value.WriteElementFlags.call_count, 1)
        self.HDF5ElementFlagValueIO.return_value.WriteElementFlags.assert_called_with(
            self.model_part.Elements)

        self.assertEqual(self.HDF5ElementGaussPointOutput.call_count, 1)
        self.assertEqual(self.HDF5ElementGaussPointOutput.call_args[0][0]['prefix'].GetString(
        ), '/ResultsData/ElementGaussPointValues')
        self.assertEqual(
            self.HDF5ElementGaussPointOutput.call_args[0][0]['list_of_variables'].size(), 0)
        self.assertEqual(
            self.HDF5ElementGaussPointOutput.return_value.WriteElementGaussPointValues.call_count, 1)
        self.HDF5ElementGaussPointOutput.return_value.WriteElementGaussPointValues.assert_called_with(
            self.model_part.Elements, self.model_part.GetCommunicator().GetDataCommunicator(), self.model_part.ProcessInfo)

        self.assertEqual(self.HDF5ConditionDataValueIO.call_count, 1)
        self.assertEqual(self.HDF5ConditionDataValueIO.call_args[0][0]['prefix'].GetString(
        ), '/ResultsData/ConditionDataValues')
        self.assertEqual(
            self.HDF5ConditionDataValueIO.call_args[0][0]['list_of_variables'].size(), 0)
        self.assertEqual(
            self.HDF5ConditionDataValueIO.return_value.WriteConditionResults.call_count, 1)
        self.HDF5ConditionDataValueIO.return_value.WriteConditionResults.assert_called_with(
            self.model_part.Conditions)

        self.assertEqual(self.HDF5ConditionFlagValueIO.call_count, 1)
        self.assertEqual(self.HDF5ConditionFlagValueIO.call_args[0][0]['prefix'].GetString(
        ), '/ResultsData/ConditionFlagValues')
        self.assertEqual(
            self.HDF5ConditionFlagValueIO.call_args[0][0]['list_of_variables'].size(), 0)
        self.assertEqual(
            self.HDF5ConditionFlagValueIO.return_value.WriteConditionFlags.call_count, 1)
        self.HDF5ConditionFlagValueIO.return_value.WriteConditionFlags.assert_called_with(
            self.model_part.Conditions)

        self.assertEqual(self.HDF5ConditionGaussPointOutput.call_count, 1)
        self.assertEqual(self.HDF5ConditionGaussPointOutput.call_args[0][0]['prefix'].GetString(
        ), '/ResultsData/ConditionGaussPointValues')
        self.assertEqual(
            self.HDF5ConditionGaussPointOutput.call_args[0][0]['list_of_variables'].size(), 0)
        self.assertEqual(
            self.HDF5ConditionGaussPointOutput.return_value.WriteConditionGaussPointValues.call_count, 1)
        self.HDF5ConditionGaussPointOutput.return_value.WriteConditionGaussPointValues.assert_called_with(
            self.model_part.Conditions, self.model_part.GetCommunicator().GetDataCommunicator(), self.model_part.ProcessInfo)

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
        for time in [0.09999999, 0.19999998]:
            self.model_part.CloneTimeStep(time)
            self.model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
            process.ExecuteFinalizeSolutionStep()
        self.assertEqual(self.HDF5FileSerial.call_count, 2)
        self.assertEqual(
            self.HDF5FileSerial.call_args[0][0]['file_name'].GetString(), 'kratos-0.2000.h5')
        self.assertEqual(
            self.HDF5FileSerial.call_args[0][0]['file_access_mode'].GetString(), 'exclusive')
        self.assertEqual(
            self.HDF5FileSerial.call_args[0][0]['echo_level'].GetInt(), 0)
        self.assertEqual(self.HDF5ModelPartIO.call_count, 2)
        self.HDF5ModelPartIO.assert_called_with(
            self.HDF5FileSerial.return_value, '/ModelData')
        self.assertEqual(
            self.HDF5ModelPartIO.return_value.WriteModelPart.call_count, 2)
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
        self.assertEqual(self.HDF5NodalSolutionStepBossakIO.call_count, 2)
        self.assertEqual(
            self.HDF5NodalSolutionStepBossakIO.call_args[0][0]['prefix'].GetString(), '/ResultsData')
        self.assertEqual(
            self.HDF5NodalSolutionStepBossakIO.call_args[0][0]['list_of_variables'].size(), 0)
        self.assertEqual(
            self.HDF5NodalSolutionStepBossakIO.call_args[0][1], self.HDF5FileSerial.return_value)
        self.assertEqual(
            self.HDF5NodalSolutionStepBossakIO.return_value.SetAlphaBossak.call_count, 2)
        self.HDF5NodalSolutionStepBossakIO.return_value.SetAlphaBossak.assert_called_with(
            -0.2)
        self.assertEqual(
            self.HDF5NodalSolutionStepBossakIO.return_value.WriteNodalResults.call_count, 2)
        self.HDF5NodalSolutionStepBossakIO.return_value.WriteNodalResults.assert_called_with(
            self.model_part)

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
            self.model_part, 0)
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
            self.model_part, 0)
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
        with patch('KratosMultiphysics.HDF5Application.xdmf_utils.WriteMultifileTemporalAnalysisToXdmf', autospec=True) as WriteMultifileTemporalAnalysisToXdmf:
            with patch('KratosMultiphysics.kratos_utilities.DeleteFileIfExisting', autospec=True) as DeleteFileIfExisting:
                with patch('pathlib.Path.glob', autospec=True) as listdir:
                    self.HDF5FileSerial().GetFileName.return_value = settings["Parameters"]["model_part_name"].GetString() + ".h5"
                    listdir.return_value = [
                        pathlib.Path('test_model_part-0.0000.h5').absolute(),
                        pathlib.Path('test_model_part-0.1000.h5').absolute()]
                    process = single_mesh_xdmf_output_process.Factory(settings, self.model)
                    process.ExecuteInitialize()
                    process.ExecuteBeforeSolutionLoop()
                    for time in [0.09999999, 0.19999998]:
                        self.model_part.CloneTimeStep(time)
                        process.ExecuteFinalizeSolutionStep()
                    self.assertEqual(WriteMultifileTemporalAnalysisToXdmf.call_count, 2)
                    WriteMultifileTemporalAnalysisToXdmf.assert_called_with('test_model_part.h5',
                                                                            '/ModelData',
                                                                            '/ResultsData')
                    DeleteFileIfExisting.assert_called_once_with(
                        str(pathlib.Path('./test_model_part-0.1000.h5').absolute()))

    def test_ImportModelPartFromHDF5Process(self):
        settings = KratosMultiphysics.Parameters('''
            {
                "Parameters": {
                    "model_part_name": "test_model_part"
                }
            }
            ''')
        process = import_model_part_from_hdf5_process.Factory(
            settings, self.model)
        process.ExecuteInitialize()
        self.assertEqual(self.HDF5FileSerial.call_count, 1)
        self.assertEqual(
            self.HDF5FileSerial.call_args[0][0]["file_name"].GetString(), "test_model_part.h5")
        self.assertEqual(
            self.HDF5FileSerial.call_args[0][0]["file_access_mode"].GetString(), "read_only")
        self.assertEqual(self.HDF5NodalSolutionStepDataIO.call_count, 1)
        self.assertEqual(
            self.HDF5NodalSolutionStepDataIO.return_value.ReadNodalResults.call_count, 1)
        self.HDF5NodalSolutionStepDataIO.return_value.ReadNodalResults.assert_called_with(
            self.model_part, 0)
        self.assertEqual(
            self.HDF5NodalDataValueIO.return_value.ReadNodalResults.call_count, 1)
        self.HDF5NodalDataValueIO.return_value.ReadNodalResults.assert_called_with(
            self.model_part.Nodes, self.model_part.GetCommunicator())
        self.assertEqual(
            self.HDF5ElementDataValueIO.return_value.ReadElementResults.call_count, 1)
        self.HDF5ElementDataValueIO.return_value.ReadElementResults.assert_called_with(
            self.model_part.Elements, self.model_part.GetCommunicator())
        self.assertEqual(
            self.HDF5NodalFlagValueIO.return_value.ReadNodalFlags.call_count, 1)
        self.HDF5NodalFlagValueIO.return_value.ReadNodalFlags.assert_called_with(
            self.model_part.Nodes, self.model_part.GetCommunicator())
        self.assertEqual(
            self.HDF5ElementFlagValueIO.return_value.ReadElementFlags.call_count, 1)
        self.HDF5ElementFlagValueIO.return_value.ReadElementFlags.assert_called_with(
            self.model_part.Elements, self.model_part.GetCommunicator())
        self.assertEqual(
            self.HDF5ConditionDataValueIO.return_value.ReadConditionResults.call_count, 1)
        self.HDF5ConditionDataValueIO.return_value.ReadConditionResults.assert_called_with(
            self.model_part.Conditions, self.model_part.GetCommunicator())
        self.assertEqual(
            self.HDF5NodalFlagValueIO.return_value.ReadNodalFlags.call_count, 1)
        self.HDF5NodalFlagValueIO.return_value.ReadNodalFlags.assert_called_with(
            self.model_part.Nodes, self.model_part.GetCommunicator())
        self.assertEqual(
            self.HDF5ConditionFlagValueIO.return_value.ReadConditionFlags.call_count, 1)
        self.HDF5ConditionFlagValueIO.return_value.ReadConditionFlags.assert_called_with(
            self.model_part.Conditions, self.model_part.GetCommunicator())


    @KratosUnittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
    def test_OutputProcess(self):
        """Test whether HDF5 output processes conform to the OutputProcess concept."""
        # Define a scoped mdpa file context
        class ScopedMDPA:
            def __init__(self, model_part_name: str):
                self.model_part_name = model_part_name
                KratosMultiphysics.ModelPartIO(self.model_part_name, KratosMultiphysics.IO.WRITE).WriteModelPart(KratosMultiphysics.Model().CreateModelPart(model_part_name))

            def __enter__(self) -> None:
                pass

            def __exit__(self, *args) -> None:
                # Delete all files with the model part in their names
                for file_path in pathlib.Path(".").glob("*{}*".format(self.model_part_name)):
                    kratos_utilities.DeleteFileIfExisting(str(file_path))

        settings = KratosMultiphysics.Parameters('''
            {
                "problem_data" : {
                    "problem_name" : "test_OutputProcess",
                    "start_time" : 0.0,
                    "end_time" : 5.0,
                    "echo_level" : 0,
                    "parallel_type" : "OpenMP"
                },
                "solver_settings" : {
                    "model_part_name" : "test_OutputProcess",
                    "domain_size" : 2,
                    "solver_type" : "dynamic",
                    "time_integration_method" : "explicit",
                    "time_stepping" : {
                        "time_step" : 0.5
                    },
                    "model_import_settings" : {
                        "input_type" : "mdpa",
                        "input_filename" : "test_OutputProcess"
                    }
                },
                "processes" : {},
                "output_processes" : {
                    "hdf5_output" : [{
                    "python_module" : "single_mesh_temporal_output_process",
                    "kratos_module" : "KratosMultiphysics.HDF5Application",
                    "process_name" : "SingleMeshTemporalOutputProcess",
                    "Parameters": {
                        "model_part_name": "test_OutputProcess",
                        "file_settings": {
                            "file_access_mode": "truncate",
                            "echo_level": 1
                        },
                        "model_part_output_settings": {
                            "prefix": "/ModelData/<model_part_name>"
                        },
                        "output_time_settings": {
                            "step_frequency" : 2
                        }
                    }
                }]}}
            ''')

        with patch("KratosMultiphysics.HDF5Application.core.controllers.Controller.ExecuteOperation") as mocked_execute:
            self.assertEqual(mocked_execute.call_count, 0)
            with ScopedMDPA("test_OutputProcess"):
                from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
                model = KratosMultiphysics.Model()
                simulation = StructuralMechanicsAnalysis(model, settings)
                simulation.Run()
            self.assertEqual(mocked_execute.call_count, 1 + 5)


if __name__ == "__main__":
    KratosUnittest.main()
