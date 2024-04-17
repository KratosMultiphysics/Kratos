from math import sqrt
import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
import KratosMultiphysics.StructuralMechanicsApplication as KratosStruct
import KratosMultiphysics.KratosUnittest as UnitTest

from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import GetSensors
class TestDisplacementSensor(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("Test")
        cls.model_part.AddNodalSolutionStepVariable(Kratos.DISPLACEMENT)

        cls.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
        cls.model_part.CreateNewNode(4, 0.0, 1.0, 0.0)

        prop = cls.model_part.CreateNewProperties(1)

        cls.model_part.CreateNewElement("Element2D3N", 1, [1, 2, 4], prop)
        cls.model_part.CreateNewElement("Element2D3N", 2, [2, 3, 4], prop)

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
                "name"         : "disp_x_1",
                "value"        : 0,
                "location"     : [0.3333333333333, 0.3333333333333, 0.0],
                "direction"    : [1.0, 1.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "displacement_sensor",
                "name"         : "disp_x_2",
                "value"        : 0,
                "location"     : [0.6666666666667, 0.6666666666667, 0.0],
                "direction"    : [1.0, 1.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }""")
        ]

        cls.sensors = GetSensors(cls.model_part, parameters)
        cls.ref_values = [7/3, 3, (7/3 + 10/3)/sqrt(2), (3 + 4)/sqrt(2)]

    def test_CalculateValue(self):
        values = [sensor.CalculateValue(self.model_part) for sensor in self.sensors]
        self.assertVectorAlmostEqual(values, self.ref_values, 7)

    def test_CalculateGradient(self):
        residual_matrix = Kratos.Matrix(18, 18)
        response_sensitivities = Kratos.Vector()
        for i, sensor in enumerate(self.sensors):
            ref_value = self.ref_values[i]
            delta = 1e-5

            element: Kratos.Element = self.model_part.GetElement(sensor.GetValue(KratosSI.SENSOR_ELEMENT_ID))
            sensor.CalculateGradient(element, residual_matrix, response_sensitivities, self.model_part.ProcessInfo)
            for i, node in enumerate(element.GetGeometry()):
                node.SetSolutionStepValue(Kratos.DISPLACEMENT_X, node.GetSolutionStepValue(Kratos.DISPLACEMENT_X) + delta)
                perturbed_value = sensor.CalculateValue(self.model_part)
                sensitivity = (perturbed_value - ref_value) / delta
                self.assertAlmostEqual(sensitivity, response_sensitivities[i * 6])
                node.SetSolutionStepValue(Kratos.DISPLACEMENT_X, node.GetSolutionStepValue(Kratos.DISPLACEMENT_X) - delta)

                node.SetSolutionStepValue(Kratos.DISPLACEMENT_Y, node.GetSolutionStepValue(Kratos.DISPLACEMENT_Y) + delta)
                perturbed_value = sensor.CalculateValue(self.model_part)
                sensitivity = (perturbed_value - ref_value) / delta
                self.assertAlmostEqual(sensitivity, response_sensitivities[i * 6 + 1])
                node.SetSolutionStepValue(Kratos.DISPLACEMENT_Y, node.GetSolutionStepValue(Kratos.DISPLACEMENT_Y) - delta)


class TestStrainSensor(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("Test")
        cls.model_part.ProcessInfo[KratosStruct.PERTURBATION_SIZE] = 1e-10
        cls.model_part.AddNodalSolutionStepVariable(Kratos.DISPLACEMENT)
        cls.model_part.AddNodalSolutionStepVariable(KratosStruct.ADJOINT_DISPLACEMENT)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.ROTATION)
        cls.model_part.AddNodalSolutionStepVariable(KratosStruct.ADJOINT_ROTATION)

        cls.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
        cls.model_part.CreateNewNode(4, 0.0, 1.0, 0.0)

        prop = cls.model_part.CreateNewProperties(1)

        cls.model_part.CreateNewElement("ShellThinElement3D3N", 1, [1, 2, 4], prop)
        cls.model_part.CreateNewElement("ShellThinElement3D3N", 2, [2, 3, 4], prop)

        # set the constitutive laws
        material_settings = Kratos.Parameters("""{"Parameters": {"materials_filename": "auxiliary_files/shell_material_properties.json"}}""")
        Kratos.ReadMaterialsUtility(material_settings, cls.model)

        for element in cls.model_part.Elements:
            element.Initialize(cls.model_part.ProcessInfo)

        Kratos.VariableUtils().AddDof(Kratos.DISPLACEMENT_X, cls.model_part)
        Kratos.VariableUtils().AddDof(Kratos.DISPLACEMENT_Y, cls.model_part)
        Kratos.VariableUtils().AddDof(Kratos.DISPLACEMENT_Z, cls.model_part)
        Kratos.VariableUtils().AddDof(Kratos.ROTATION_X, cls.model_part)
        Kratos.VariableUtils().AddDof(Kratos.ROTATION_Y, cls.model_part)
        Kratos.VariableUtils().AddDof(Kratos.ROTATION_Z, cls.model_part)
        Kratos.VariableUtils().AddDof(KratosStruct.ADJOINT_DISPLACEMENT_X, cls.model_part)
        Kratos.VariableUtils().AddDof(KratosStruct.ADJOINT_DISPLACEMENT_Y, cls.model_part)
        Kratos.VariableUtils().AddDof(KratosStruct.ADJOINT_DISPLACEMENT_Z, cls.model_part)
        Kratos.VariableUtils().AddDof(KratosStruct.ADJOINT_ROTATION_X, cls.model_part)
        Kratos.VariableUtils().AddDof(KratosStruct.ADJOINT_ROTATION_Y, cls.model_part)
        Kratos.VariableUtils().AddDof(KratosStruct.ADJOINT_ROTATION_Z, cls.model_part)

        for node in cls.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.DISPLACEMENT, [node.Id, node.Id + 1, node.Id + 2])
            node.SetSolutionStepValue(Kratos.ROTATION, [node.Id * 2, node.Id * 2 + 1, node.Id * 2 + 2])

        parameters = [
            Kratos.Parameters("""{
                "type"           : "strain_sensor",
                "name"           : "strain_x_1",
                "value"          : 0,
                "location"       : [0.3333333333333, 0.3333333333333, 0.0],
                "strain_type"    : "strain_xx",
                "strain_variable": "SHELL_STRAIN",
                "weight"         : 1.0,
                "variable_data"  : {}
            }"""),
            Kratos.Parameters("""{

                "type"           : "strain_sensor",
                "name"           : "strain_x_2",
                "value"          : 0,
                "location"       : [0.6666666666667, 0.6666666666667, 0.0],
                "strain_type"    : "strain_xx",
                "strain_variable": "SHELL_STRAIN",
                "weight"         : 1.0,
                "variable_data"  : {}
            }"""),
            Kratos.Parameters("""{

                "type"           : "strain_sensor",
                "name"           : "strain_y_1",
                "value"          : 0,
                "location"       : [0.3333333333333, 0.3333333333333, 0.0],
                "strain_type"    : "strain_yy",
                "strain_variable": "SHELL_STRAIN",
                "weight"         : 1.0,
                "variable_data"  : {}
            }"""),
            Kratos.Parameters("""{

                "type"           : "strain_sensor",
                "name"           : "strain_y_2",
                "value"          : 0,
                "location"       : [0.6666666666667, 0.6666666666667, 0.0],
                "strain_type"    : "strain_yy",
                "strain_variable": "SHELL_STRAIN",
                "weight"         : 1.0,
                "variable_data"  : {}
            }""")
        ]

        cls.sensors = GetSensors(cls.model_part, parameters)
        cls.ref_values = [0.5, -1.5, 4.5, 0.5]

    def test_CalculateValue(self):
        values = [sensor.CalculateValue(self.model_part) for sensor in self.sensors]
        self.assertVectorAlmostEqual(values, self.ref_values, 3)

    def test_CalculateGradient(self):
        residual_matrix = Kratos.Matrix(18, 18)
        response_sensitivities = Kratos.Vector()
        for i, sensor in enumerate(self.sensors):
            ref_value = self.ref_values[i]
            delta = 1e-5

            element: Kratos.Element = self.model_part.GetElement(sensor.GetValue(KratosSI.SENSOR_ELEMENT_ID))
            sensor.CalculateGradient(element, residual_matrix, response_sensitivities, self.model_part.ProcessInfo)
            for i, node in enumerate(element.GetGeometry()):
                node.SetSolutionStepValue(Kratos.DISPLACEMENT_X, node.GetSolutionStepValue(Kratos.DISPLACEMENT_X) + delta)
                perturbed_value = sensor.CalculateValue(self.model_part)
                sensitivity = (perturbed_value - ref_value) / delta
                self.assertAlmostEqual(sensitivity, response_sensitivities[i * 6], 4)
                node.SetSolutionStepValue(Kratos.DISPLACEMENT_X, node.GetSolutionStepValue(Kratos.DISPLACEMENT_X) - delta)

                node.SetSolutionStepValue(Kratos.DISPLACEMENT_Y, node.GetSolutionStepValue(Kratos.DISPLACEMENT_Y) + delta)
                perturbed_value = sensor.CalculateValue(self.model_part)
                sensitivity = (perturbed_value - ref_value) / delta
                self.assertAlmostEqual(sensitivity, response_sensitivities[i * 6 + 1], 4)
                node.SetSolutionStepValue(Kratos.DISPLACEMENT_Y, node.GetSolutionStepValue(Kratos.DISPLACEMENT_Y) - delta)


if __name__ == '__main__':
    UnitTest.main()