import KratosMultiphysics as Kratos
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as UnitTest
from KratosMultiphysics.testing.utilities import ReadModelPart
from KratosMultiphysics.SystemIdentificationApplication.processes.sensor_output_process import SensorOutputProcess
from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting

class TestSensorOutputProcess(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3
        cls.model_part.AddNodalSolutionStepVariable(Kratos.DISPLACEMENT)
        ReadModelPart("auxiliary_files/structure", cls.model_part)

        params = Kratos.Parameters("""{
            "model_part_name" : "test",
            "output_file_name": "auxiliary_files/measured_data.csv",
            "properties_list": [
                "type",
                "name",
                "location",
                "value"
            ],
            "list_of_sensors": [
                {
                    "type": "displacement_sensor",
                    "location": [
                        58.54103333333333,
                        15.060899999999998,
                        0.0
                    ],
                    "direction": [
                        1.0,
                        0.0,
                        0.0
                    ],
                    "name": "disp_x_1",
                    "weight": 1.0
                },
                {
                    "type": "displacement_sensor",
                    "name": "disp_x_2",
                    "location": [
                        59.34613333333333,
                        26.242666666666665,
                        0.0
                    ],
                    "direction": [
                        1.0,
                        0.0,
                        0.0
                    ],
                    "weight": 1.0
                }
            ]
        }""")
        cls.sensor_output_process = SensorOutputProcess(cls.model, params)

        for node in cls.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.DISPLACEMENT, Kratos.Array3([node.Id, node.Id + 1, node.Id + 2]))

    def test_PrintOutput(self):
        self.addCleanup(DeleteFileIfExisting, "auxiliary_files/measured_data.csv")
        self.sensor_output_process.PrintOutput()

        with open("auxiliary_files/measured_data.csv", "r") as file_input:
            lines = file_input.readlines()
        with open("auxiliary_files/measured_data_ref.csv", "r") as file_input:
            ref_lines = file_input.readlines()
        self.assertEqual(lines, ref_lines)

if __name__ == '__main__':
    UnitTest.main()