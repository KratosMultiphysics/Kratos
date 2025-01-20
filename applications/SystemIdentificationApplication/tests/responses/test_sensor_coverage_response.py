import numpy as np
import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import CreateSensors
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import SetSensors
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import AddMaskStatusController
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import GetMaskStatusControllers
from KratosMultiphysics.SystemIdentificationApplication.responses.sensor_coverage_response import SensorCoverageResponse
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView

class TestSensorCoverageResponse(UnitTest.TestCase):
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

        exps_list: 'list[Kratos.Expression.ElementExpression]' = []
        for i, sensor in enumerate(cls.sensors):
            node = sensor.GetNode()
            node.SetValue(KratosSI.SENSOR_STATUS, (node.Id % 3) / 2)

            elem_np = np.zeros(cls.domain_model_part.NumberOfElements())
            for j in range(elem_np.shape[0]):
                if j % (i + 2) == 0:
                    elem_np[j] = 1
            elem_exp = Kratos.Expression.ElementExpression(cls.domain_model_part)
            Kratos.Expression.CArrayExpressionIO.Read(elem_exp, elem_np)
            sensor.AddContainerExpression("mask_exp", elem_exp)
            exps_list.append(elem_exp.Clone())

        # add the mask status controller
        AddMaskStatusController(cls.sensor_group_data, "mask_exp", KratosSI.SensorMaskStatus(cls.sensor_model_part, exps_list, 0))

        params = Kratos.Parameters("""{
            "sensor_group_name": "sensors",
            "sensor_mask_name" : "mask_exp",
            "echo_level"       : 0
        }""")

        cls.response = SensorCoverageResponse("test", cls.model, params, cls.optimization_problem)
        cls.response.Initialize()

    def test_CalculateValue(self):
        domain_size_exp = Kratos.Expression.ElementExpression(self.domain_model_part)
        Kratos.Expression.DomainSizeExpressionIO.Read(domain_size_exp)

        total_mask = Kratos.Expression.ElementExpression(self.domain_model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(total_mask, 0)
        for sensor in self.sensors:
            sensor_status = sensor.GetNode().GetValue(KratosSI.SENSOR_STATUS)
            mask = sensor.GetContainerExpression("mask_exp")
            total_mask += mask * sensor_status

        total_mask = KratosSI.ElementSmoothClamper(0, 1).ProjectForward(total_mask)
        for mask_status_controller in GetMaskStatusControllers(self.sensor_group_data, "mask_exp"):
            mask_status_controller.Update()
        response_value = self.response.CalculateValue()
        self.assertNotEqual(response_value, 0.0)
        self.assertAlmostEqual(response_value, Kratos.Expression.Utils.InnerProduct(domain_size_exp, total_mask) / Kratos.Expression.Utils.Sum(domain_size_exp))

    def test_CalculateGradient(self):
        for mask_status_controller in GetMaskStatusControllers(self.sensor_group_data, "mask_exp"):
            mask_status_controller.Update()

        ref_value = self.response.CalculateValue()
        collective_exp = KratosOA.CollectiveExpression()
        collective_exp.Add(Kratos.Expression.NodalExpression(self.sensor_model_part))
        self.response.CalculateGradient({KratosSI.SENSOR_STATUS: collective_exp})
        analytical_gradient = collective_exp.GetContainerExpressions()[0].Evaluate()

        delta = 1e-8
        for i, node in enumerate(self.sensor_model_part.Nodes):
            node.SetValue(KratosSI.SENSOR_STATUS, node.GetValue(KratosSI.SENSOR_STATUS) + delta)
            for mask_status_controller in GetMaskStatusControllers(self.sensor_group_data, "mask_exp"):
                    mask_status_controller.Update()
            fd_sensitivity = (self.response.CalculateValue() - ref_value) / delta
            node.SetValue(KratosSI.SENSOR_STATUS, node.GetValue(KratosSI.SENSOR_STATUS) - delta)
            self.assertAlmostEqual(fd_sensitivity, analytical_gradient[i])

if __name__ == '__main__':
    UnitTest.main()