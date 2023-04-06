
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.OptimizationApplication.responses.mass_response_function import MassResponseFunction

class TestStandardizedResponses(kratos_unittest.TestCase):
    @classmethod
    def GetParameters(cls):
        return Kratos.Parameters("""{
            "evaluated_model_part_names": ["test"]
        }""")

    @classmethod
    def CreateElements(cls):
        cls.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
        cls.model_part.CreateNewNode(4, 0.0, 1.0, 0.0)

        properties = cls.model_part.CreateNewProperties(1)
        properties[Kratos.DENSITY] = 2.0
        properties[Kratos.THICKNESS] = 3.0
        cls.model_part.CreateNewElement("Element2D3N", 1, [1, 2, 3], properties)

        properties = cls.model_part.CreateNewProperties(2)
        properties[Kratos.DENSITY] = 4.0
        properties[Kratos.THICKNESS] = 6.0
        cls.model_part.CreateNewElement("Element2D3N", 2, [4, 1, 3], properties)

    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")

        cls.response_function = MassResponseFunction(cls.model, cls.GetParameters(), None)
        cls.CreateElements()

        cls.response_function.Initialize()
        cls.response_function.Check()
        cls.ref_value = cls.response_function.CalculateValue()

    def _CalculateSensitivity(self, sensitivity_variable):
        self.response_function.CalculateSensitivity({sensitivity_variable: [self.model_part]})

    def _CheckSensitivity(self, response_function, entities, sensitivity_method, update_method, delta, precision):
        ref_value = response_function.CalculateValue()
        for entity in entities:
            v = sensitivity_method(entity)
            update_method(entity, delta)
            value = response_function.CalculateValue()
            sensitivity = (value - ref_value)/delta
            update_method(entity, -delta)
            self.assertAlmostEqual(v, sensitivity, precision)

    def _UpdateProperties(self, variable, entity, delta):
        entity.Properties[variable] += delta

    def _UpdateNodalPositions(self, direction, entity, delta):
        if direction == 0:
            entity.X += delta
        if direction == 1:
            entity.Y += delta
        if direction == 2:
            entity.Z += delta
