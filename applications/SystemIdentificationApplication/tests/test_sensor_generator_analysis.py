import KratosMultiphysics as Kratos
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as UnitTest
from KratosMultiphysics.SystemIdentificationApplication.sensor_generator_analysis import SensorGeneratorAnalysis
from KratosMultiphysics.compare_two_files_check_process import CompareTwoFilesCheckProcess

class TestSensorGeneratorAnalysis(UnitTest.TestCase):
    def test_MeshBased(self):
        parameters = Kratos.Parameters("""{
            "model_part_name"     : "test.ElemMat4",
            "mdpa_file_name"      : "auxiliary_files/structure",
            "output_file_name"    : "sensor_data.json",
            "csv_output_file_name": "sensor_data.csv",
            "sensor_groups"   : [
                {
                    "generation_method": "mesh_based",
                    "container_type"   : "elements",
                    "location"         : "center",
                    "sensor_settings"  : {
                        "name"    : "<ENTITY_ID>_dummy",
                        "PRESSURE": 1.0,
                        "STEP"    : 2,
                        "VELOCITY": [1, 1, 1]
                    }
                }
            ]
        }""")

        model = Kratos.Model()
        SensorGeneratorAnalysis(model, parameters).Run()

        check_params = Kratos.Parameters("""{
            "reference_file_name"   : "auxiliary_files/sensor_data_mesh_based_ref.csv",
            "output_file_name"      : "sensor_data.csv",
            "remove_output_file"    : true,
            "comparison_type"       : "deterministic",
            "tolerance"             : 1e-6,
            "relative_tolerance"    : 1e-9,
            "dimension"             : 3
        }""")
        CompareTwoFilesCheckProcess(check_params).Execute()

        check_params = Kratos.Parameters("""{
            "reference_file_name"   : "auxiliary_files/sensor_data_mesh_based_ref.json",
            "output_file_name"      : "sensor_data.json",
            "remove_output_file"    : true,
            "comparison_type"       : "deterministic",
            "tolerance"             : 1e-6,
            "relative_tolerance"    : 1e-9,
            "dimension"             : 3
        }""")
        CompareTwoFilesCheckProcess(check_params).Execute()

    def test_BoundingSurfaceBased(self):
        parameters = Kratos.Parameters("""{
            "model_part_name"     : "test.ElemMat2",
            "mdpa_file_name"      : "auxiliary_files/structure",
            "output_file_name"    : "sensor_data.json",
            "csv_output_file_name": "sensor_data.csv",
            "sensor_groups"   : [
                {
                    "generation_method"        : "bounding_surface",
                    "bounding_surface_corner_1": [15.0, 30.0, 0.0],
                    "bounding_surface_corner_2": [20.0, 25.0, 0.0],
                    "number_of_sensors"        : [3, 3, 1],
                    "sensor_settings"  : {
                        "name"    : "<ENTITY_ID>_dummy",
                        "PRESSURE": 1.0,
                        "STEP"    : 2,
                        "VELOCITY": [1, 1, 1]
                    }
                }
            ]
        }""")

        model = Kratos.Model()
        SensorGeneratorAnalysis(model, parameters).Run()

        check_params = Kratos.Parameters("""{
            "reference_file_name"   : "auxiliary_files/sensor_data_bb_based_ref.csv",
            "output_file_name"      : "sensor_data.csv",
            "remove_output_file"    : true,
            "comparison_type"       : "deterministic",
            "tolerance"             : 1e-6,
            "relative_tolerance"    : 1e-9,
            "dimension"             : 3
        }""")
        CompareTwoFilesCheckProcess(check_params).Execute()

        check_params = Kratos.Parameters("""{
            "reference_file_name"   : "auxiliary_files/sensor_data_bb_based_ref.json",
            "output_file_name"      : "sensor_data.json",
            "remove_output_file"    : true,
            "comparison_type"       : "deterministic",
            "tolerance"             : 1e-6,
            "relative_tolerance"    : 1e-9,
            "dimension"             : 3
        }""")
        CompareTwoFilesCheckProcess(check_params).Execute()


if __name__ == '__main__':
    UnitTest.main()