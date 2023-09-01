from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils
from KratosMultiphysics.compare_two_files_check_process import CompareTwoFilesCheckProcess

from KratosMultiphysics.HDF5Application.core.xdmf_dataset_generator import SingleMeshMultiFileSameDatasetsGenerator
from KratosMultiphysics.HDF5Application.core.xdmf_dataset_generator import SingleFileDatasetsGenerator
from KratosMultiphysics.HDF5Application.core.xdmf_dataset_generator import MultiFileDatasetsGenerator
from KratosMultiphysics.HDF5Application.xdmf_utils import WriteDataSetsToXdmf

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

        for entity in cls.model_part.Elements:
            entity.SetValue(Kratos.PRESSURE, entity.Id)
            entity.SetValue(Kratos.VELOCITY, [entity.Id, entity.Id + 1, entity.Id + 2])

    def setUp(self) -> None:
        self.model_part.ProcessInfo[Kratos.STEP] = 0

    def testSingleMeshMultiFileSameDatasetsGenerator(self):
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

            process = SingleMeshTemporalOutputProcess(self.model, parameters)
            process.ExecuteInitialize()
            process.Check()
            for t in range(1, number_of_time_steps + 1, 1):
                self.model_part.ProcessInfo[Kratos.STEP] = t
                if process.IsOutputStep():
                    process.PrintOutput()

            generator = SingleMeshMultiFileSameDatasetsGenerator("test-<step>.h5")
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
                    "output_control_type": "step",
                    "output_interval"    : 1.0
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

            process = SingleMeshTemporalOutputProcess(self.model, parameters)
            process.ExecuteInitialize()
            process.Check()
            for t in range(1, number_of_time_steps + 1, 1):
                self.model_part.ProcessInfo[Kratos.STEP] = t
                if process.IsOutputStep():
                    process.PrintOutput()

            generator = SingleFileDatasetsGenerator("test.h5:/Results-<step>")
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

if __name__ == "__main__":
    KratosUnittest.main()