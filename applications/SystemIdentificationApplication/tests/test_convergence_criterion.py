import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
import KratosMultiphysics.KratosUnittest as UnitTest
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import CreateSensors
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import SetSensors
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import GetSensors
import KratosMultiphysics.SystemIdentificationApplication.convergence_criteria.sensor_error as sensor_error
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView

class TestSensorErrorConvCriterion(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("Test")
        cls.sensor_model_part = cls.model.CreateModelPart("SensorModelPart")
        cls.model_part.AddNodalSolutionStepVariable(Kratos.DISPLACEMENT)

        cls.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
        cls.model_part.CreateNewNode(4, 0.0, 1.0, 0.0)

        prop = cls.model_part.CreateNewProperties(1)

        cls.model_part.CreateNewElement("Element2D3N", 1, [1, 2, 4], prop).SetValue(Kratos.PRESSURE, 10.0)
        cls.model_part.CreateNewElement("Element2D3N", 2, [2, 3, 4], prop).SetValue(Kratos.PRESSURE, 100.0)

        for node in cls.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.DISPLACEMENT, [node.Id, node.Id + 1, node.Id + 2])

        parameters = [
            Kratos.Parameters("""{

                "type"         : "displacement_sensor",
                "name"         : "disp_x_1",
                "value"        : 0,
                "location"     : [0.3333333333333, 0.3333333333333, 0.0],
                "direction"    : [1.0, 0.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "displacement_sensor",
                "name"         : "disp_x_2",
                "value"        : 0,
                "location"     : [0.6666666666667, 0.6666666666667, 0.0],
                "direction"    : [1.0, 0.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "displacement_sensor",
                "name"         : "disp_x_3",
                "value"        : 0,
                "location"     : [0.3333333333333, 0.3333333333333, 0.0],
                "direction"    : [1.0, 1.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "displacement_sensor",
                "name"         : "disp_x_4",
                "value"        : 0,
                "location"     : [0.6666666666667, 0.6666666666667, 0.0],
                "direction"    : [1.0, 1.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }""")
        ]

        cls.sensors = CreateSensors(cls.sensor_model_part, cls.model_part, parameters)

        for sensor in cls.sensors:
            sensor.GetNode().SetValue(KratosSI.SENSOR_ERROR, sensor.GetNode().Id * 1e-7 + 2e-7)

        cls.optimization_problem = OptimizationProblem()

    def test_SensorMaxErrorCriterionFalse(self):
        param = Kratos.Parameters("""{
                    "sensor_group_name": "SensorModelPart",
                    "tolerance"     : 1e-9,
                    "output_to_file" : false
                }""")
        SetSensors(ComponentDataView("SensorModelPart", self.optimization_problem), self.sensors)
        convergence_criterium = sensor_error.MaxSensorErrorCriterion(param, self.optimization_problem)
        convergence_criterium.Initialize()
        self.assertFalse(convergence_criterium.IsConverged())
        convergence_criterium.Finalize()

    def test_SensorMaxErrorCriterionTrue(self):
        param = Kratos.Parameters("""{
                    "sensor_group_name": "SensorModelPart",
                    "tolerance"     : 1e-4,
                    "output_to_file" : false
                }""")
        convergence_criterium = sensor_error.MaxSensorErrorCriterion(param, self.optimization_problem)
        convergence_criterium.Initialize()
        self.assertTrue(convergence_criterium.IsConverged())
        convergence_criterium.Finalize()


if __name__ == '__main__':
    UnitTest.main()