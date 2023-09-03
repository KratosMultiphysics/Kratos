from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.HDF5Application as KratosHDF5
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils
from KratosMultiphysics.compare_two_files_check_process import CompareTwoFilesCheckProcess

from KratosMultiphysics.HDF5Application.core.dataset_generator import SingleMeshMultiFileSameDatasetsGenerator
from KratosMultiphysics.HDF5Application.core.dataset_generator import SingleFileDatasetsGenerator
from KratosMultiphysics.HDF5Application.core.dataset_generator import MultiFileDatasetsGenerator
from KratosMultiphysics.HDF5Application.xdmf_utils import WriteDataSetsToXdmf

from KratosMultiphysics.HDF5Application.core.processes import HDF5OutputProcess
from KratosMultiphysics.HDF5Application.single_mesh_temporal_output_process import SingleMeshTemporalOutputProcess
from KratosMultiphysics.HDF5Application.multiple_mesh_temporal_output_process import MultipleMeshTemporalOutputProcess

class TestDatasetGenerator(KratosUnittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")

        cls.model_part.CreateNewNode(1, 0, 0, 0)
        cls.model_part.CreateNewNode(2, 1, 0, 0)
        cls.model_part.CreateNewNode(3, 0, 1, 0)

        properties = cls.model_part.CreateNewProperties(1)
        cls.model_part.CreateNewElement("Element2D3N", 1, [1, 2, 3], properties)
        cls.model_part.CreateNewCondition("LineCondition2D2N", 1, [1, 2], properties)

        for entity in cls.model_part.Nodes:
            entity.SetValue(Kratos.PRESSURE, entity.Id)
            entity.SetValue(Kratos.VELOCITY, [entity.Id, entity.Id + 1, entity.Id + 2])

        for entity in cls.model_part.Conditions:
            entity.SetValue(Kratos.PRESSURE, entity.Id)
            entity.SetValue(Kratos.VELOCITY, [entity.Id, entity.Id + 1, entity.Id + 2])

        for entity in cls.model_part.Elements:
            entity.SetValue(Kratos.PRESSURE, entity.Id)
            entity.SetValue(Kratos.VELOCITY, [entity.Id, entity.Id + 1, entity.Id + 2])

    def setUp(self) -> None:
        self.model_part.ProcessInfo[Kratos.TIME] = 0
        self.model_part.ProcessInfo[Kratos.STEP] = 0

    def testSingleMeshMultiFileSameDatasetsGenerator(self):
        number_of_time_steps = 5

        for t in range(1, number_of_time_steps + 1, 1):
            self.addCleanup(kratos_utils.DeleteFileIfExisting, str((Path(__file__).parent / f"auxiliar_files/xdmf/test-{t:0.7f}.h5").absolute()))

        with KratosUnittest.WorkFolderScope("auxiliar_files/xdmf", __file__):
            parameters = Kratos.Parameters("""{
                "model_part_name": "test",
                "file_settings": {
                    "file_name"        : "<model_part_name>-<time>.h5",
                    "file_access_mode" : "exclusive",
                    "time_format"      : "0.7f"
                },
                "output_time_settings": {
                    "output_control_type": "time",
                    "output_interval"    : 1.0
                },
                "nodal_data_value_settings": {
                    "list_of_variables": ["PRESSURE", "VELOCITY"]
                },
                "element_data_value_settings": {
                    "list_of_variables": ["PRESSURE", "VELOCITY"]
                },
                "condition_data_value_settings": {
                    "list_of_variables": ["PRESSURE", "VELOCITY"]
                }
            }""")

            process = SingleMeshTemporalOutputProcess(self.model, parameters)
            process.ExecuteInitialize()
            process.Check()
            for t in range(1, number_of_time_steps + 1, 1):
                self.model_part.ProcessInfo[Kratos.TIME] = t
                if process.IsOutputStep():
                    process.PrintOutput()

            generator = SingleMeshMultiFileSameDatasetsGenerator("test-<time>.h5")
            WriteDataSetsToXdmf(generator, "single_mesh1.xdmf")
            CompareTwoFilesCheckProcess(Kratos.Parameters("""
            {
                "reference_file_name"   : "single_mesh1.orig.xdmf",
                "output_file_name"      : "single_mesh1.xdmf",
                "remove_output_file"    : true,
                "comparison_type"       : "deterministic"
            }""")).Execute()

        with self.assertRaises(RuntimeError):
            generator = SingleMeshMultiFileSameDatasetsGenerator("test-1.h5")

        with self.assertRaises(RuntimeError):
            generator = SingleMeshMultiFileSameDatasetsGenerator("test-1:/ModelData<step>.h5")

        with self.assertRaises(RuntimeError):
            generator = SingleMeshMultiFileSameDatasetsGenerator("test-<step>:/ModelData<step>.h5")

    def testSingleFileDatasetsGenerator1(self):
        number_of_time_steps = 5
        self.addCleanup(kratos_utils.DeleteFileIfExisting, str((Path(__file__).parent / f"auxiliar_files/xdmf/test.h5").absolute()))

        with KratosUnittest.WorkFolderScope("auxiliar_files/xdmf", __file__):
            parameters = Kratos.Parameters("""{
                "model_part_name": "test",
                "file_settings": {
                    "file_name"        : "<model_part_name>.h5",
                    "file_access_mode" : "read_write"
                },
                "output_time_settings": {
                    "output_control_type": "time",
                    "output_interval"    : 1.0
                },
                "nodal_data_value_settings": {
                    "prefix": "/Results-<time>/Nodes-<time>/",
                    "list_of_variables": ["PRESSURE", "VELOCITY"]
                },
                "element_data_value_settings": {
                    "prefix": "/Results-<time>/Elements-<time>/",
                    "list_of_variables": ["PRESSURE", "VELOCITY"]
                },
                "condition_data_value_settings": {
                    "prefix": "/Results-<time>/Conditions-<time>/",
                    "list_of_variables": ["PRESSURE", "VELOCITY"]
                }
            }""")

            process = SingleMeshTemporalOutputProcess(self.model, parameters)
            process.ExecuteInitialize()
            process.Check()
            for t in range(1, number_of_time_steps + 1, 1):
                self.model_part.ProcessInfo[Kratos.TIME] = t
                if process.IsOutputStep():
                    process.PrintOutput()

            generator = SingleFileDatasetsGenerator("test.h5:/Results-<time>")
            WriteDataSetsToXdmf(generator, "single_mesh2.xdmf")
            CompareTwoFilesCheckProcess(Kratos.Parameters("""
            {
                "reference_file_name"   : "single_mesh2.orig.xdmf",
                "output_file_name"      : "single_mesh2.xdmf",
                "remove_output_file"    : true,
                "comparison_type"       : "deterministic"
            }""")).Execute()

        with self.assertRaises(RuntimeError):
            generator = SingleFileDatasetsGenerator("test-<step>.h5")

        with self.assertRaises(RuntimeError):
            generator = SingleFileDatasetsGenerator("test-<step>:/ModelData.h5")

        with self.assertRaises(RuntimeError):
            generator = SingleFileDatasetsGenerator("test-<step>:/ModelData<step>.h5")

    def testSingleFileDatasetsGenerator2(self):
        number_of_time_steps = 5
        self.addCleanup(kratos_utils.DeleteFileIfExisting, str((Path(__file__).parent / f"auxiliar_files/xdmf/test.h5").absolute()))

        with KratosUnittest.WorkFolderScope("auxiliar_files/xdmf", __file__):
            parameters = Kratos.Parameters("""{
                "model_part_name": "test",
                "file_settings": {
                    "file_name"        : "<model_part_name>.h5",
                    "file_access_mode" : "read_write"
                },
                "output_time_settings": {
                    "output_control_type": "step",
                    "output_interval"    : 1.0
                },
                "model_part_output_settings": {
                    "prefix": "/ModelData-<step>/<step>"
                },
                "nodal_data_value_settings": {
                    "prefix": "/Results-<step>/Nodes-<step>/",
                    "list_of_variables": ["PRESSURE", "VELOCITY"]
                },
                "element_data_value_settings": {
                    "prefix": "/Results-<step>/Elements-<step>/",
                    "list_of_variables": ["PRESSURE", "VELOCITY"]
                },
                "condition_data_value_settings": {
                    "prefix": "/Results-<step>/Conditions-<step>/",
                    "list_of_variables": ["PRESSURE", "VELOCITY"]
                }
            }""")

            process = MultipleMeshTemporalOutputProcess(self.model, parameters)
            process.ExecuteInitialize()
            process.Check()
            for t in range(1, number_of_time_steps + 1, 1):
                self.model_part.ProcessInfo[Kratos.STEP] = t
                if process.IsOutputStep():
                    process.PrintOutput()

            generator = SingleFileDatasetsGenerator("test.h5:/Results-<step>")
            WriteDataSetsToXdmf(generator, "single_mesh3.xdmf")
            CompareTwoFilesCheckProcess(Kratos.Parameters("""
            {
                "reference_file_name"   : "single_mesh3.orig.xdmf",
                "output_file_name"      : "single_mesh3.xdmf",
                "remove_output_file"    : true,
                "comparison_type"       : "deterministic"
            }""")).Execute()

    def testMultiFileDatasetsGenerator1(self):
        number_of_time_steps = 5

        for t in range(1, number_of_time_steps + 1, 1):
            self.addCleanup(kratos_utils.DeleteFileIfExisting, str((Path(__file__).parent / f"auxiliar_files/xdmf/test-{t:0.7f}.h5").absolute()))

        with KratosUnittest.WorkFolderScope("auxiliar_files/xdmf", __file__):
            parameters = Kratos.Parameters("""{
                "model_part_name": "test",
                "file_settings": {
                    "file_name"        : "<model_part_name>-<time>.h5",
                    "file_access_mode" : "exclusive",
                    "time_format"      : "0.7f"
                },
                "output_time_settings": {
                    "output_control_type": "time",
                    "output_interval"    : 1.0
                },
                "nodal_data_value_settings": {
                    "list_of_variables": ["PRESSURE", "VELOCITY"]
                },
                "element_data_value_settings": {
                    "list_of_variables": ["PRESSURE", "VELOCITY"]
                },
                "condition_data_value_settings": {
                    "list_of_variables": ["PRESSURE", "VELOCITY"]
                }
            }""")

            process = SingleMeshTemporalOutputProcess(self.model, parameters)
            process.ExecuteInitialize()
            process.Check()
            for t in range(1, number_of_time_steps + 1, 1):
                self.model_part.ProcessInfo[Kratos.TIME] = t
                if process.IsOutputStep():
                    process.PrintOutput()

            generator = MultiFileDatasetsGenerator("test-<time>.h5")
            WriteDataSetsToXdmf(generator, "single_mesh4.xdmf")
            CompareTwoFilesCheckProcess(Kratos.Parameters("""
            {
                "reference_file_name"   : "single_mesh1.orig.xdmf",
                "output_file_name"      : "single_mesh4.xdmf",
                "remove_output_file"    : true,
                "comparison_type"       : "deterministic"
            }""")).Execute()

        with self.assertRaises(RuntimeError):
            generator = SingleMeshMultiFileSameDatasetsGenerator("test-1.h5")

        with self.assertRaises(RuntimeError):
            generator = SingleMeshMultiFileSameDatasetsGenerator("test-1:/ModelData<step>.h5")

        with self.assertRaises(RuntimeError):
            generator = SingleMeshMultiFileSameDatasetsGenerator("test-<step>:/ModelData<step>.h5")

    def testMultiFileDatasetsGenerator2(self):
        number_of_time_steps = 5

        for t in range(1, number_of_time_steps + 1, 1):
            self.addCleanup(kratos_utils.DeleteFileIfExisting, str((Path(__file__).parent / f"auxiliar_files/xdmf/test-{t}.h5").absolute()))

        with KratosUnittest.WorkFolderScope("auxiliar_files/xdmf", __file__):
            parameters = Kratos.Parameters("""{
                "model_part_name": "test",
                "file_settings": {
                    "file_name"        : "<model_part_name>-<step>.h5",
                    "file_access_mode" : "exclusive"
                },
                "output_time_settings": {
                    "output_control_type": "step",
                    "output_interval"    : 1.0
                },
                "nodal_data_value_settings": {
                    "list_of_variables": ["PRESSURE", "VELOCITY"]
                },
                "element_data_value_settings": {
                    "list_of_variables": ["PRESSURE", "VELOCITY"]
                },
                "condition_data_value_settings": {
                    "list_of_variables": ["PRESSURE", "VELOCITY"]
                }
            }""")

            process = MultipleMeshTemporalOutputProcess(self.model, parameters)
            process.ExecuteInitialize()
            process.Check()
            for t in range(1, number_of_time_steps + 1, 1):
                self.model_part.ProcessInfo[Kratos.STEP] = t
                if process.IsOutputStep():
                    process.PrintOutput()

            generator = MultiFileDatasetsGenerator("test-<step>.h5")
            WriteDataSetsToXdmf(generator, "single_mesh5.xdmf")
            CompareTwoFilesCheckProcess(Kratos.Parameters("""
            {
                "reference_file_name"   : "single_mesh4.orig.xdmf",
                "output_file_name"      : "single_mesh5.xdmf",
                "remove_output_file"    : true,
                "comparison_type"       : "deterministic"
            }""")).Execute()

    def testSingleMeshMultiFileSameDatasetsGenerator(self):
        number_of_time_steps = 5

        for t in range(1, number_of_time_steps + 1, 1):
            self.addCleanup(kratos_utils.DeleteFileIfExisting, str((Path(__file__).parent / f"auxiliar_files/xdmf/test-{t:0.7f}.h5").absolute()))

        nodal_cexp = Kratos.Expression.NodalExpression(self.model_part)
        cond_cexp = Kratos.Expression.ConditionExpression(self.model_part)
        elem_cexp = Kratos.Expression.ElementExpression(self.model_part)

        Kratos.Expression.VariableExpressionIO.Read(nodal_cexp, Kratos.VELOCITY, False)
        Kratos.Expression.VariableExpressionIO.Read(cond_cexp, Kratos.VELOCITY)
        Kratos.Expression.VariableExpressionIO.Read(elem_cexp, Kratos.VELOCITY)

        with KratosUnittest.WorkFolderScope("auxiliar_files/xdmf", __file__):
            process = HDF5OutputProcess(self.model_part)
            process_id = process.GetProcessId()

            parameters = Kratos.Parameters("""{
                "file_name"        : "<model_part_name>-<time>.h5",
                "file_access_mode" : "exclusive",
                "custom_attributes": {
                    "__hdf5_process_id": []
                }
            }""")
            parameters["custom_attributes"]["__hdf5_process_id"].Append(process_id[0])
            parameters["custom_attributes"]["__hdf5_process_id"].Append(process_id[1])
            model_part_written = False
            for t in range(1, number_of_time_steps + 1, 1):
                parameters["file_name"].SetString(f"test-{t:0.7f}.h5")
                h5_file = KratosHDF5.HDF5File(self.model_part.GetCommunicator().GetDataCommunicator(), parameters)
                if not model_part_written:
                    model_param = Kratos.Parameters("""{
                        "prefix": "/ModelData",
                        "custom_attributes": {
                            "__hdf5_process_id": []
                        }
                    }""")
                    model_param["custom_attributes"]["__hdf5_process_id"].Append(process_id[0])
                    model_param["custom_attributes"]["__hdf5_process_id"].Append(process_id[1])
                    KratosHDF5.HDF5ModelPartIO(model_param, h5_file).WriteModelPart(self.model_part)
                    model_part_written = True

                nodal_cexp += nodal_cexp
                cond_cexp += cond_cexp
                elem_cexp += elem_cexp

                attribs = Kratos.Parameters("""{
                    "__hdf5_process_id": []
                }""")
                attribs["__hdf5_process_id"].Append(process_id[0])
                attribs["__hdf5_process_id"].Append(process_id[1])
                KratosHDF5.ExpressionIO(Kratos.Parameters("""{"prefix":"/Results/"}"""), h5_file).Write("nodal_values", nodal_cexp, attribs)
                KratosHDF5.ExpressionIO(Kratos.Parameters("""{"prefix":"/Results/"}"""), h5_file).Write("cond_values", cond_cexp, attribs)
                KratosHDF5.ExpressionIO(Kratos.Parameters("""{"prefix":"/Results/"}"""), h5_file).Write("elem_values", elem_cexp, attribs)

            generator = SingleMeshMultiFileSameDatasetsGenerator("test-<time>.h5")
            WriteDataSetsToXdmf(generator, "single_mesh5.xdmf")
            CompareTwoFilesCheckProcess(Kratos.Parameters("""
            {
                "reference_file_name"   : "single_mesh5.orig.xdmf",
                "output_file_name"      : "single_mesh5.xdmf",
                "remove_output_file"    : true,
                "comparison_type"       : "deterministic"
            }""")).Execute()

if __name__ == "__main__":
    KratosUnittest.main()