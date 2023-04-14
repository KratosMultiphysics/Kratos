
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.OptimizationApplication.responses.mass_response_function import MassResponseFunction

class TestMassResponseFunctionBase(kratos_unittest.TestCase):
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

class TestMassResponseFunctionBeams(TestMassResponseFunctionBase):
    @classmethod
    def GetParameters(cls):
        return Kratos.Parameters("""{
            "evaluated_model_part_names": ["test"]
        }""")

    @classmethod
    def CreateElements(cls):
        cls.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(3, 4.0, 4.0, 0.0)

        for i in range(2):
            node_ids = [(i % 3) + 1, ((i + 1) % 3) + 1]
            properties = cls.model_part.CreateNewProperties(i)
            properties[Kratos.DENSITY] = 2.0 * (i + 1)
            properties[KratosOA.CROSS_AREA] = 3.0 * (i + 1)
            cls.model_part.CreateNewElement("Element2D2N", i, node_ids, properties)

    def test_CalculateValue(self):
        self.assertAlmostEqual(self.ref_value, 126, 12)

    def test_CalculateShapeSensitivity(self):
        self._CalculateSensitivity(Kratos.SHAPE_SENSITIVITY)
        # calculate nodal shape sensitivities
        self._CheckSensitivity(
            self.response_function,
            self.model_part.Nodes,
            lambda x: x.GetValue(Kratos.SHAPE_SENSITIVITY_X),
            lambda x, y: self._UpdateNodalPositions(0, x, y),
            1e-6,
            4)

        self._CheckSensitivity(
            self.response_function,
            self.model_part.Nodes,
            lambda x: x.GetValue(Kratos.SHAPE_SENSITIVITY_Y),
            lambda x, y: self._UpdateNodalPositions(1, x, y),
            1e-6,
            4)

    def test_CalculateDensitySensitivity(self):
        self._CalculateSensitivity(KratosOA.DENSITY_SENSITIVITY)

        # calculate element density sensitivity
        self._CheckSensitivity(
            self.response_function,
            self.model_part.Elements,
            lambda x: x.Properties[KratosOA.DENSITY_SENSITIVITY],
            lambda x, y: self._UpdateProperties(Kratos.DENSITY, x, y),
            1e-6,
            6)

    def test_CalculateCrossAreaSensitivity(self):
        self._CalculateSensitivity(KratosOA.CROSS_AREA_SENSITIVITY)

        # calculate element cross area sensitivity
        self._CheckSensitivity(
            self.response_function,
            self.model_part.Elements,
            lambda x: x.Properties[KratosOA.CROSS_AREA_SENSITIVITY],
            lambda x, y: self._UpdateProperties(KratosOA.CROSS_AREA, x, y),
            1e-6,
            6)

class TestMassResponseFunctionShells(TestMassResponseFunctionBase):
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

    def test_CalculateValue(self):
        self.assertAlmostEqual(self.ref_value, 15, 12)

    def test_CalculateShapeSensitivity(self):
        self._CalculateSensitivity(Kratos.SHAPE_SENSITIVITY)

        # calculate nodal shape sensitivities
        self._CheckSensitivity(
            self.response_function,
            self.model_part.Nodes,
            lambda x: x.GetValue(Kratos.SHAPE_SENSITIVITY_X),
            lambda x, y: self._UpdateNodalPositions(0, x, y),
            1e-6,
            4)

        self._CheckSensitivity(
            self.response_function,
            self.model_part.Nodes,
            lambda x: x.GetValue(Kratos.SHAPE_SENSITIVITY_Y),
            lambda x, y: self._UpdateNodalPositions(1, x, y),
            1e-6,
            4)

    def test_CalculateDensitySensitivity(self):
        self._CalculateSensitivity(KratosOA.DENSITY_SENSITIVITY)

        # calculate element density sensitivity
        self._CheckSensitivity(
            self.response_function,
            self.model_part.Elements,
            lambda x: x.Properties[KratosOA.DENSITY_SENSITIVITY],
            lambda x, y: self._UpdateProperties(Kratos.DENSITY, x, y),
            1e-6,
            6)

    def test_CalculateThicknessSensitivity(self):
        self._CalculateSensitivity(KratosOA.THICKNESS_SENSITIVITY)

        # calculate element cross area sensitivity
        self._CheckSensitivity(
            self.response_function,
            self.model_part.Elements,
            lambda x: x.Properties[KratosOA.THICKNESS_SENSITIVITY],
            lambda x, y: self._UpdateProperties(Kratos.THICKNESS, x, y),
            1e-7,
            6)

class TestMassResponseFunctionSolids(TestMassResponseFunctionBase):
    @classmethod
    def GetParameters(cls):
        return Kratos.Parameters("""{
            "evaluated_model_part_names": ["test"]
        }""")

    @classmethod
    def CreateElements(cls):
        cls.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(2, 2.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(3, 2.0, 2.0, 0.0)
        cls.model_part.CreateNewNode(4, 2.0, 2.0, 1.0)
        cls.model_part.CreateNewNode(5, 2.0, 2.0, -1.0)

        properties = cls.model_part.CreateNewProperties(1)
        properties[Kratos.DENSITY] = 2.0
        cls.model_part.CreateNewElement("Element3D4N", 1, [1, 2, 3, 4], properties)

        properties = cls.model_part.CreateNewProperties(2)
        properties[Kratos.DENSITY] = 4.0
        cls.model_part.CreateNewElement("Element3D4N", 2, [5, 1, 2, 3], properties)

    def test_CalculateValue(self):
        v = 0.0
        for element in self.model_part.Elements:
            v += element.GetGeometry().DomainSize() * element.Properties[Kratos.DENSITY]
        self.assertAlmostEqual(self.ref_value, v, 12)

    def test_CalculateShapeSensitivity(self):
        self._CalculateSensitivity(Kratos.SHAPE_SENSITIVITY)
        # calculate nodal shape sensitivities
        self._CheckSensitivity(
            self.response_function,
            self.model_part.Nodes,
            lambda x: x.GetValue(Kratos.SHAPE_SENSITIVITY_X),
            lambda x, y: self._UpdateNodalPositions(0, x, y),
            1e-6,
            4)

        self._CheckSensitivity(
            self.response_function,
            self.model_part.Nodes,
            lambda x: x.GetValue(Kratos.SHAPE_SENSITIVITY_Y),
            lambda x, y: self._UpdateNodalPositions(1, x, y),
            1e-6,
            4)

        self._CheckSensitivity(
            self.response_function,
            self.model_part.Nodes,
            lambda x: x.GetValue(Kratos.SHAPE_SENSITIVITY_Z),
            lambda x, y: self._UpdateNodalPositions(2, x, y),
            1e-6,
            4)

    def test_CalculateDensitySensitivity(self):
        self._CalculateSensitivity(KratosOA.DENSITY_SENSITIVITY)

        # calculate element density sensitivity
        self._CheckSensitivity(
            self.response_function,
            self.model_part.Elements,
            lambda x: x.Properties[KratosOA.DENSITY_SENSITIVITY],
            lambda x, y: self._UpdateProperties(Kratos.DENSITY, x, y),
            1e-6,
            6)


if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()