import numpy as np
import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.KratosUnittest as UnitTest
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import CreateSensors
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import SetSensors
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import AddMaskStatusController
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import GetMaskStatusControllers
from KratosMultiphysics.SystemIdentificationApplication.responses.sensor_localization_response import SensorLocalizationResponse

class TestSensorLocalizationResponse(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
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
        cls.domain_model_part = cls.model.CreateModelPart("mask")
        cls.domain_model_part.AddNodalSolutionStepVariable(Kratos.DISPLACEMENT)

        cls.domain_model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        cls.domain_model_part.CreateNewNode(2, 2.0, 0.0, 0.0)
        cls.domain_model_part.CreateNewNode(3, 4.0, 0.0, 0.0)
        cls.domain_model_part.CreateNewNode(4, 6.0, 0.0, 0.0)
        cls.domain_model_part.CreateNewNode(5, 0.0, 2.0, 0.0)
        cls.domain_model_part.CreateNewNode(6, 2.0, 2.0, 0.0)
        cls.domain_model_part.CreateNewNode(7, 4.0, 2.0, 0.0)
        cls.domain_model_part.CreateNewNode(8, 6.0, 2.0, 0.0)
        cls.domain_model_part.CreateNewNode(9, 0.0, 4.0, 0.0)
        cls.domain_model_part.CreateNewNode(10, 2.0, 4.0, 0.0)
        cls.domain_model_part.CreateNewNode(11, 4.0, 4.0, 0.0)
        cls.domain_model_part.CreateNewNode(12, 6.0, 4.0, 0.0)

        prop = cls.domain_model_part.CreateNewProperties(1)

        cls.domain_model_part.CreateNewElement("Element2D4N", 1, [1, 2, 6, 5], prop)
        cls.domain_model_part.CreateNewElement("Element2D4N", 2, [2, 3, 7, 6], prop)
        cls.domain_model_part.CreateNewElement("Element2D4N", 3, [3, 4, 8, 7], prop)
        cls.domain_model_part.CreateNewElement("Element2D4N", 4, [5, 6, 10, 9], prop)
        cls.domain_model_part.CreateNewElement("Element2D4N", 5, [6, 7, 11, 10], prop)
        cls.domain_model_part.CreateNewElement("Element2D4N", 6, [7, 8, 12, 11], prop)

        for node in cls.domain_model_part.Nodes:
            node.SetSolutionStepValue(Kratos.DISPLACEMENT, [node.Id, node.Id + 1, node.Id + 2])

        parameters = [
            Kratos.Parameters("""{

                "type"         : "displacement_sensor",
                "name"         : "disp_x_1",
                "value"        : 0,
                "location"     : [1, 1, 0.0],
                "direction"    : [1.0, 0.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "displacement_sensor",
                "name"         : "disp_x_2",
                "value"        : 0,
                "location"     : [3, 1, 0.0],
                "direction"    : [1.0, 0.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "displacement_sensor",
                "name"         : "disp_x_3",
                "value"        : 0,
                "location"     : [3, 3, 0],
                "direction"    : [1.0, 1.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "displacement_sensor",
                "name"         : "disp_x_4",
                "value"        : 0,
                "location"     : [5, 3, 0.0],
                "direction"    : [1.0, 1.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }""")
        ]

        cls.optimization_problem = OptimizationProblem()

        cls.sensor_model_part = cls.model.CreateModelPart("sensors")
        cls.sensors = CreateSensors(cls.sensor_model_part, cls.domain_model_part, parameters)

        cls.sensor_group_data = ComponentDataView("sensors", cls.optimization_problem)

        SetSensors(cls.sensor_group_data, cls.sensors)

        for sensor in cls.sensors:
            sensor.GetNode().SetValue(KratosSI.SENSOR_STATUS, (sensor.GetNode().Id % 2) / 2)

        elem_exp = Kratos.Expression.ElementExpression(cls.domain_model_part)

        elem_np = np.array([1, 1, 0, 1, 1, 0], dtype=np.float64)
        Kratos.Expression.CArrayExpressionIO.Read(elem_exp, elem_np)
        cls.sensors[0].AddContainerExpression("mask_exp", elem_exp.Clone())

        elem_np = np.array([1, 1, 0, 1, 0, 1], dtype=np.float64)
        Kratos.Expression.CArrayExpressionIO.Read(elem_exp, elem_np)
        cls.sensors[1].AddContainerExpression("mask_exp", elem_exp.Clone())

        elem_np = np.array([0, 0, 0, 1, 0, 0], dtype=np.float64)
        Kratos.Expression.CArrayExpressionIO.Read(elem_exp, elem_np)
        cls.sensors[2].AddContainerExpression("mask_exp", elem_exp.Clone())

        elem_np = np.array([0, 0, 1, 0, 1, 1], dtype=np.float64)
        Kratos.Expression.CArrayExpressionIO.Read(elem_exp, elem_np)
        cls.sensors[3].AddContainerExpression("mask_exp", elem_exp.Clone())

        # add the mask status controller
        cls.sensor_mask_status = KratosSI.SensorMaskStatus(cls.sensor_model_part, [sensor.GetContainerExpression("mask_exp").Clone() for sensor in cls.sensors], 0)
        cls.sensor_mask_status_kd_tree = KratosSI.SensorMaskStatusKDTree(cls.sensor_mask_status, 100, 0)
        AddMaskStatusController(cls.sensor_group_data, "mask_exp", cls.sensor_mask_status)
        AddMaskStatusController(cls.sensor_group_data, "mask_exp", cls.sensor_mask_status_kd_tree)

        params = Kratos.Parameters("""{
            "sensor_group_name": "sensors",
            "sensor_mask_name" : "mask_exp",
            "beta"             : 500
        }""")

        cls.response = SensorLocalizationResponse("test", cls.model, params, cls.optimization_problem)
        cls.response.Initialize()

    def test_CalculateValue1(self):
        self.sensor_model_part.GetNode(1).SetValue(KratosSI.SENSOR_STATUS, 1)
        self.sensor_model_part.GetNode(2).SetValue(KratosSI.SENSOR_STATUS, 0)
        self.sensor_model_part.GetNode(3).SetValue(KratosSI.SENSOR_STATUS, 1)
        self.sensor_model_part.GetNode(4).SetValue(KratosSI.SENSOR_STATUS, 0)

        """
        status  1   0   1   0
        mask    m1  m2  m3  m4
        e1      1   1   0   0
        e2      1   1   0   0
        e3      0   0   0   1
        e4      1   1   1   0
        e5      1   0   0   1
        e6      0   1   0   1

        So the clustering will be based on m1 and m3

        c_1 = e1, e2, e5
        c_2 = e1, e2, e5
        c_3 = e3, e6
        c_4 = e4
        c_5 = e1, e2, e5
        c_6 = e3, e6
        """
        self.__UpdateStatusDependents()
        self.assertAlmostEqual(0.5, self.response.CalculateValue())

        self.sensor_model_part.GetNode(1).SetValue(KratosSI.SENSOR_STATUS, 1)
        self.sensor_model_part.GetNode(2).SetValue(KratosSI.SENSOR_STATUS, 1)
        self.sensor_model_part.GetNode(3).SetValue(KratosSI.SENSOR_STATUS, 1)
        self.sensor_model_part.GetNode(4).SetValue(KratosSI.SENSOR_STATUS, 1)

        """
        status  1   1   1   1
        mask    m1  m2  m3  m4
        e1      1   1   0   0
        e2      1   1   0   0
        e3      0   0   0   1
        e4      1   1   1   0
        e5      1   0   0   1
        e6      0   1   0   1

        So the clustering will be based on m1, m2, m3, m4

        c_1 = e1, e2
        c_2 = e1, e2
        c_3 = e3
        c_4 = e4
        c_5 = e5
        c_6 = e6
        """

        self.__UpdateStatusDependents()
        self.assertAlmostEqual(1 / 3, self.response.CalculateValue())

    def test_CalculateValue2(self):
        for node in self.sensor_model_part.Nodes:
            node.SetValue(KratosSI.SENSOR_STATUS, (node.Id % 3) / 2)
        self.__UpdateStatusDependents()
        self.assertAlmostEqual(self.response.CalculateValue(), 0.5)

    def test_CalculateValue3(self):
        params = Kratos.Parameters("""{
            "sensor_group_name": "sensors",
            "sensor_mask_name" : "mask_exp",
            "beta"             : 503
        }""")
        response = SensorLocalizationResponse("test1", self.model, params, self.optimization_problem)
        response.Initialize()

        self.sensor_model_part.GetNode(1).SetValue(KratosSI.SENSOR_STATUS, 1)
        self.sensor_model_part.GetNode(2).SetValue(KratosSI.SENSOR_STATUS, 0)
        self.sensor_model_part.GetNode(3).SetValue(KratosSI.SENSOR_STATUS, 1)
        self.sensor_model_part.GetNode(4).SetValue(KratosSI.SENSOR_STATUS, 0)

        """
        status  1   0   1   0
        mask    m1  m2  m3  m4
        e1      1   1   0   0
        e2      1   1   0   0
        e3      0   0   0   1
        e4      1   1   1   0
        e5      1   0   0   1
        e6      0   1   0   1

        So the clustering will be based on m1 and m3

        c_1 = e1, e2, e5
        c_2 = e1, e2, e5
        c_3 = e3, e6
        c_4 = e4
        c_5 = e1, e2, e5
        c_6 = e3, e6
        """
        self.__UpdateStatusDependents()

        # the exact max cluster size ratio is 0.5
        self.assertAlmostEqual(response.CalculateValue(), 0.5)

    def test_CalculateGradient1(self):
        for node in self.sensor_model_part.Nodes:
            node.SetValue(KratosSI.SENSOR_STATUS, (node.Id % 4) / 4)
        self.__UpdateStatusDependents()

        ref_value = self.response.CalculateValue()
        collective_exp = KratosOA.CollectiveExpression()
        collective_exp.Add(Kratos.Expression.NodalExpression(self.sensor_model_part))
        self.response.CalculateGradient({KratosSI.SENSOR_STATUS: collective_exp})
        analytical_gradient = collective_exp.GetContainerExpressions()[0].Evaluate()

        delta = 1e-8
        for i, node in enumerate(self.sensor_model_part.Nodes):
            node.SetValue(KratosSI.SENSOR_STATUS, node.GetValue(KratosSI.SENSOR_STATUS) + delta)
            self.__UpdateStatusDependents()
            fd_sensitivity = (self.response.CalculateValue() - ref_value) / delta
            node.SetValue(KratosSI.SENSOR_STATUS, node.GetValue(KratosSI.SENSOR_STATUS) - delta)
            self.assertAlmostEqual(fd_sensitivity, analytical_gradient[i], 5)

    def __UpdateStatusDependents(self) -> None:
        for mask_status_controller in GetMaskStatusControllers(self.sensor_group_data, "mask_exp"):
            mask_status_controller.Update()

if __name__ == '__main__':
    UnitTest.main()