import numpy as np
import KratosMultiphysics as Kratos
import KratosMultiphysics.DigitalTwinApplication as KratosDT
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.KratosUnittest as UnitTest
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import GetSensors
from KratosMultiphysics.DigitalTwinApplication.controls.sensor_status_control import SensorStatusControl
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView

class TestSensorStatusControl(UnitTest.TestCase):
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

        for node in cls.mask_model_part.Nodes:
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

        cls.sensors = GetSensors(cls.mask_model_part, parameters)
        cls.sensor_model_part = cls.model.CreateModelPart("sensors")
        for i, sensor in enumerate(cls.sensors):
            loc = sensor.GetLocation()
            node: Kratos.Node = cls.sensor_model_part.CreateNewNode(i + 1, loc[0], loc[1], loc[2])
            node.SetValue(KratosDT.SENSOR_STATUS, (node.Id % 3) / 2)

            elem_np = np.zeros(cls.mask_model_part.NumberOfElements())
            for j in range(elem_np.shape[0]):
                if j % (i + 2) == 0:
                    elem_np[j] = 1
            elem_exp = Kratos.Expression.ElementExpression(cls.mask_model_part)
            Kratos.Expression.CArrayExpressionIO.Read(elem_exp, elem_np)
            sensor.AddElementExpression("mask_exp", elem_exp)

        params = Kratos.Parameters("""{
            "controlled_model_part_names": ["sensors"],
            "mask_expression_name"       : "mask_exp",
            "beta_settings": {
                "initial_value": 5,
                "max_value"    : 30,
                "adaptive"     : true,
                "increase_fac" : 1.05,
                "update_period": 50
            }
        }""")

        cls.optimization_problem = OptimizationProblem()
        ComponentDataView("sensor", cls.optimization_problem).GetUnBufferedData().SetValue("list_of_sensors", cls.sensors)
        cls.control = SensorStatusControl("test", cls.model, params, cls.optimization_problem)
        cls.optimization_problem.AddComponent(cls.control)
        cls.control.Initialize()
        cls.control.Check()

        cls.initial_control_field = cls.control.GetControlField()

    def setUp(self) -> None:
        self.control.Update(self.initial_control_field)

    def test_Update(self):
        current_beta = self.control.beta
        control_field = self.control.GetControlField()
        ref_sensor_status = Kratos.Expression.NodalExpression(self.sensor_model_part)
        Kratos.Expression.VariableExpressionIO.Read(ref_sensor_status, KratosDT.SENSOR_STATUS, False)
        self.control.Update(control_field)

        self.assertAlmostEqual(current_beta, self.control.beta)
        sensor_status = Kratos.Expression.NodalExpression(self.sensor_model_part)
        Kratos.Expression.VariableExpressionIO.Read(sensor_status, KratosDT.SENSOR_STATUS, False)
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(sensor_status - ref_sensor_status), 0.0)

        control_field *= 2.0
        for i in range(50):
            self.optimization_problem.AdvanceStep()
        self.control.Update(control_field)
        self.assertAlmostEqual(current_beta * 1.05, self.control.beta)
        Kratos.Expression.VariableExpressionIO.Read(sensor_status, KratosDT.SENSOR_STATUS, False)
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(sensor_status), 1.0)

    def test_MapGradient(self):
        sensitivity_field = Kratos.Expression.NodalExpression(self.sensor_model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(sensitivity_field, 2.0)
        mapped_gradient = self.control.MapGradient({KratosDT.SENSOR_STATUS: sensitivity_field})
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(mapped_gradient), 12.37796844095402)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.control.Finalize()

if __name__ == '__main__':
    UnitTest.main()