from math import sqrt
import KratosMultiphysics as Kratos
import KratosMultiphysics.DigitalTwinApplication as KratosDT
import KratosMultiphysics.kratos_utilities as kratos_utils
import KratosMultiphysics.KratosUnittest as UnitTest
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_io import OpenSensorFile
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import GetSensors
from KratosMultiphysics.DigitalTwinApplication.sensor_io.sensor_data_model_part_controller import SensorDataModelPartController
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView

class TestSensorDataModelPartController(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.mask_model_part = cls.model.CreateModelPart("mask")

        """
          (0,0)
            1------2------3------4
            |      |      |      |
            |   1  |  2   |  3   |
            |      |      |      |
            5------6------7------8
            |      |      |      |
            |   4  |  5   |  6   |
            |      |      |      |
            9-----10-----11-----12
                                (6,4)
        """

        cls.model = Kratos.Model()
        cls.mask_model_part = cls.model.CreateModelPart("mask")
        cls.mask_model_part.AddNodalSolutionStepVariable(Kratos.DISPLACEMENT)

        cls.mask_model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        cls.mask_model_part.CreateNewNode(2, 2.0, 0.0, 0.0)
        cls.mask_model_part.CreateNewNode(3, 4.0, 0.0, 0.0)
        cls.mask_model_part.CreateNewNode(4, 6.0, 0.0, 0.0)
        cls.mask_model_part.CreateNewNode(5, 0.0, 2.0, 0.0)
        cls.mask_model_part.CreateNewNode(6, 2.0, 2.0, 0.0)
        cls.mask_model_part.CreateNewNode(7, 4.0, 2.0, 0.0)
        cls.mask_model_part.CreateNewNode(8, 6.0, 2.0, 0.0)
        cls.mask_model_part.CreateNewNode(9, 0.0, 4.0, 0.0)
        cls.mask_model_part.CreateNewNode(10, 2.0, 4.0, 0.0)
        cls.mask_model_part.CreateNewNode(11, 4.0, 4.0, 0.0)
        cls.mask_model_part.CreateNewNode(12, 6.0, 4.0, 0.0)

        prop = cls.mask_model_part.CreateNewProperties(1)

        cls.mask_model_part.CreateNewElement("Element2D4N", 1, [1, 2, 6, 5], prop)
        cls.mask_model_part.CreateNewElement("Element2D4N", 2, [2, 3, 7, 6], prop)
        cls.mask_model_part.CreateNewElement("Element2D4N", 3, [3, 4, 8, 7], prop)
        cls.mask_model_part.CreateNewElement("Element2D4N", 4, [5, 6, 10, 9], prop)
        cls.mask_model_part.CreateNewElement("Element2D4N", 5, [6, 7, 11, 10], prop)
        cls.mask_model_part.CreateNewElement("Element2D4N", 6, [7, 8, 12, 11], prop)

        params = Kratos.Parameters("""{
            "target_model_part_name": "mask",
            "sensor_data_file_name" : "sensor_data.h5",
            "sensor_data_prefix"    : "/sensor_data",
            "sensor_model_part_name": "sensors",
            "list_of_sensors"       : [
                {
                    "type"         : "displacement_sensor",
                    "name"         : "disp_x_1",
                    "value"        : 0,
                    "location"     : [1, 1, 0.0],
                    "direction"    : [1.0, 0.0, 0.0],
                    "weight"       : 1.0,
                    "variable_data": {}
                },
                {
                    "type"         : "displacement_sensor",
                    "name"         : "disp_x_2",
                    "value"        : 0,
                    "location"     : [3, 1, 0.0],
                    "direction"    : [1.0, 0.0, 0.0],
                    "weight"       : 1.0,
                    "variable_data": {}
                },
                {
                    "type"         : "displacement_sensor",
                    "name"         : "disp_x_3",
                    "value"        : 0,
                    "location"     : [3, 3],
                    "direction"    : [1.0, 1.0, 0.0],
                    "weight"       : 1.0,
                    "variable_data": {}
                },
                {
                    "type"         : "displacement_sensor",
                    "name"         : "disp_x_4",
                    "value"        : 0,
                    "location"     : [5, 3, 0.0],
                    "direction"    : [1.0, 1.0, 0.0],
                    "weight"       : 1.0,
                    "variable_data": {}
                }
            ]
        }""")

        sensors = GetSensors(cls.mask_model_part, params["list_of_sensors"].values())
        with OpenSensorFile(cls.mask_model_part, "sensor_data.h5", "/sensor_data", "w") as sensor_io:
            for i, sensor in enumerate(sensors):
                mask = Kratos.Expression.NodalExpression(cls.mask_model_part)
                Kratos.Expression.LiteralExpressionIO.SetData(mask, i)
                sensor.AddContainerExpression("mask", mask)
                sensor_io.Write(sensor)

        cls.optimization_problem = OptimizationProblem()
        cls.sensor_data_model_part_controller = SensorDataModelPartController(cls.model, params, cls.optimization_problem)
        cls.sensor_data_model_part_controller.Initialize()

    def test_ImportModelPart(self):
        self.sensor_data_model_part_controller.ImportModelPart()
        list_of_sensors: 'list[KratosDT.Sensors.Sensor]' = ComponentDataView("sensor", self.optimization_problem).GetUnBufferedData().GetValue("list_of_sensors")
        self.assertEqual(len(list_of_sensors), 4)

        for i, sensor in enumerate(list_of_sensors):
            mask = sensor.GetNodalExpression("mask")
            self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(mask), sqrt((i ** 2) * self.mask_model_part.NumberOfNodes()))

    @classmethod
    def tearDownClass(cls) -> None:
        cls.sensor_data_model_part_controller.Finalize()
        kratos_utils.DeleteFileIfExisting("sensor_data.h5")

    def test_Test(self):
        pass

if __name__ == '__main__':
    UnitTest.main()
