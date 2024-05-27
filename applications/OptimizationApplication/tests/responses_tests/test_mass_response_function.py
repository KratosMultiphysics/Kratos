from abc import ABC, abstractmethod
import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.OptimizationApplication.responses.mass_response_function import MassResponseFunction
import KratosMultiphysics.StructuralMechanicsApplication

class TestMassResponseFunctionBase(kratos_unittest.TestCase, ABC):
    @classmethod
    @abstractmethod
    def GetParameters(cls) -> Kratos.Parameters:
        pass

    @classmethod
    @abstractmethod
    def CreateElements(cls) -> None:
        pass

    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")

        cls.response_function = MassResponseFunction("mass", cls.model, cls.GetParameters())
        cls.CreateElements()

        cls.response_function.Initialize()
        cls.response_function.Check()
        cls.ref_value = cls.response_function.CalculateValue()

    def _CheckSensitivity(self, response_function: MassResponseFunction, entities, sensitivity_method, update_method, container_expression_data, delta, precision):
        ref_value = response_function.CalculateValue()
        list_of_sensitivities = []
        for entity in entities:
            v = sensitivity_method(entity)
            update_method(entity, delta)
            value = response_function.CalculateValue()
            sensitivity = (value - ref_value)/delta
            update_method(entity, -delta)
            self.assertAlmostEqual(v, sensitivity, precision)
            list_of_sensitivities.append(v)

        list_of_sensitivities = numpy.array(list_of_sensitivities)
        self.assertVectorAlmostEqual(list_of_sensitivities, container_expression_data, precision)

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
    def GetParameters(cls) -> Kratos.Parameters:
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
        sensitivity = KratosOA.CollectiveExpression([Kratos.Expression.NodalExpression(self.model_part)])
        self.response_function.CalculateGradient({KratosOA.SHAPE: sensitivity})

        # calculate nodal shape sensitivities
        self._CheckSensitivity(
            self.response_function,
            self.model_part.Nodes,
            lambda x: x.GetValue(Kratos.SHAPE_SENSITIVITY_X),
            lambda x, y: self._UpdateNodalPositions(0, x, y),
            sensitivity.GetContainerExpressions()[0].Evaluate()[:, 0],
            1e-6,
            4)

        self._CheckSensitivity(
            self.response_function,
            self.model_part.Nodes,
            lambda x: x.GetValue(Kratos.SHAPE_SENSITIVITY_Y),
            lambda x, y: self._UpdateNodalPositions(1, x, y),
            sensitivity.GetContainerExpressions()[0].Evaluate()[:, 1],
            1e-6,
            4)

    def test_CalculateDensitySensitivity(self):
        sensitivity = KratosOA.CollectiveExpression([Kratos.Expression.ElementExpression(self.model_part)])
        self.response_function.CalculateGradient({Kratos.DENSITY: sensitivity})

        # calculate element density sensitivity
        self._CheckSensitivity(
            self.response_function,
            self.model_part.Elements,
            lambda x: x.Properties[KratosOA.DENSITY_SENSITIVITY],
            lambda x, y: self._UpdateProperties(Kratos.DENSITY, x, y),
            sensitivity.GetContainerExpressions()[0].Evaluate(),
            1e-6,
            6)

    def test_CalculateCrossAreaSensitivity(self):
        sensitivity = KratosOA.CollectiveExpression([Kratos.Expression.ElementExpression(self.model_part)])
        self.response_function.CalculateGradient({KratosOA.CROSS_AREA: sensitivity})

        # calculate element cross area sensitivity
        self._CheckSensitivity(
            self.response_function,
            self.model_part.Elements,
            lambda x: x.Properties[KratosOA.CROSS_AREA_SENSITIVITY],
            lambda x, y: self._UpdateProperties(KratosOA.CROSS_AREA, x, y),
            sensitivity.GetContainerExpressions()[0].Evaluate(),
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
        sensitivity = KratosOA.CollectiveExpression([Kratos.Expression.NodalExpression(self.model_part)])
        self.response_function.CalculateGradient({KratosOA.SHAPE: sensitivity})

        # calculate nodal shape sensitivities
        self._CheckSensitivity(
            self.response_function,
            self.model_part.Nodes,
            lambda x: x.GetValue(Kratos.SHAPE_SENSITIVITY_X),
            lambda x, y: self._UpdateNodalPositions(0, x, y),
            sensitivity.GetContainerExpressions()[0].Evaluate()[:, 0],
            1e-6,
            4)

        self._CheckSensitivity(
            self.response_function,
            self.model_part.Nodes,
            lambda x: x.GetValue(Kratos.SHAPE_SENSITIVITY_Y),
            lambda x, y: self._UpdateNodalPositions(1, x, y),
            sensitivity.GetContainerExpressions()[0].Evaluate()[:, 1],
            1e-6,
            4)

    def test_CalculateDensitySensitivity(self):
        sensitivity = KratosOA.CollectiveExpression([Kratos.Expression.ElementExpression(self.model_part)])
        self.response_function.CalculateGradient({Kratos.DENSITY: sensitivity})

        # calculate element density sensitivity
        self._CheckSensitivity(
            self.response_function,
            self.model_part.Elements,
            lambda x: x.Properties[KratosOA.DENSITY_SENSITIVITY],
            lambda x, y: self._UpdateProperties(Kratos.DENSITY, x, y),
            sensitivity.GetContainerExpressions()[0].Evaluate(),
            1e-6,
            6)

    def test_CalculateThicknessSensitivity(self):
        sensitivity = KratosOA.CollectiveExpression([Kratos.Expression.ElementExpression(self.model_part)])
        self.response_function.CalculateGradient({Kratos.THICKNESS: sensitivity})

        # calculate element cross area sensitivity
        self._CheckSensitivity(
            self.response_function,
            self.model_part.Elements,
            lambda x: x.Properties[KratosOA.THICKNESS_SENSITIVITY],
            lambda x, y: self._UpdateProperties(Kratos.THICKNESS, x, y),
            sensitivity.GetContainerExpressions()[0].Evaluate(),
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
        sensitivity = KratosOA.CollectiveExpression([Kratos.Expression.NodalExpression(self.model_part)])
        self.response_function.CalculateGradient({KratosOA.SHAPE: sensitivity})

        # calculate nodal shape sensitivities
        self._CheckSensitivity(
            self.response_function,
            self.model_part.Nodes,
            lambda x: x.GetValue(Kratos.SHAPE_SENSITIVITY_X),
            lambda x, y: self._UpdateNodalPositions(0, x, y),
            sensitivity.GetContainerExpressions()[0].Evaluate()[:, 0],
            1e-6,
            4)

        self._CheckSensitivity(
            self.response_function,
            self.model_part.Nodes,
            lambda x: x.GetValue(Kratos.SHAPE_SENSITIVITY_Y),
            lambda x, y: self._UpdateNodalPositions(1, x, y),
            sensitivity.GetContainerExpressions()[0].Evaluate()[:, 1],
            1e-6,
            4)

        self._CheckSensitivity(
            self.response_function,
            self.model_part.Nodes,
            lambda x: x.GetValue(Kratos.SHAPE_SENSITIVITY_Z),
            lambda x, y: self._UpdateNodalPositions(2, x, y),
            sensitivity.GetContainerExpressions()[0].Evaluate()[:, 2],
            1e-6,
            4)

    def test_CalculateDensitySensitivity(self):
        sensitivity = KratosOA.CollectiveExpression([Kratos.Expression.ElementExpression(self.model_part)])
        self.response_function.CalculateGradient({Kratos.DENSITY: sensitivity})

        # calculate element density sensitivity
        self._CheckSensitivity(
            self.response_function,
            self.model_part.Elements,
            lambda x: x.Properties[KratosOA.DENSITY_SENSITIVITY],
            lambda x, y: self._UpdateProperties(Kratos.DENSITY, x, y),
            sensitivity.GetContainerExpressions()[0].Evaluate(),
            1e-6,
            6)

class TestMassResponseFunctionQuads(TestMassResponseFunctionBase):
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
        cls.model_part.CreateNewNode(4, 0.0, 2.0, 0.0)

        properties = cls.model_part.CreateNewProperties(1)
        properties[Kratos.DENSITY] = 2.0
        properties[Kratos.THICKNESS] = 0.01
        cls.model_part.CreateNewElement("ShellThickElement3D4N", 1, [1, 2, 3, 4], properties)

    def test_CalculateValue(self):
        v = 0.0
        for element in self.model_part.Elements:
            v += element.GetGeometry().DomainSize() * element.Properties[Kratos.DENSITY] * element.Properties[Kratos.THICKNESS]
        self.assertAlmostEqual(self.ref_value, v, 12)

    def test_CalculateShapeSensitivity(self):
        sensitivity = KratosOA.CollectiveExpression([Kratos.Expression.NodalExpression(self.model_part)])
        self.response_function.CalculateGradient({KratosOA.SHAPE: sensitivity})

        # calculate nodal shape sensitivities
        self._CheckSensitivity(
            self.response_function,
            self.model_part.Nodes,
            lambda x: x.GetValue(Kratos.SHAPE_SENSITIVITY_X),
            lambda x, y: self._UpdateNodalPositions(0, x, y),
            sensitivity.GetContainerExpressions()[0].Evaluate()[:, 0],
            1e-6,
            4)

        self._CheckSensitivity(
            self.response_function,
            self.model_part.Nodes,
            lambda x: x.GetValue(Kratos.SHAPE_SENSITIVITY_Y),
            lambda x, y: self._UpdateNodalPositions(1, x, y),
            sensitivity.GetContainerExpressions()[0].Evaluate()[:, 1],
            1e-6,
            4)

        self._CheckSensitivity(
            self.response_function,
            self.model_part.Nodes,
            lambda x: x.GetValue(Kratos.SHAPE_SENSITIVITY_Z),
            lambda x, y: self._UpdateNodalPositions(2, x, y),
            sensitivity.GetContainerExpressions()[0].Evaluate()[:, 2],
            1e-6,
            4)

    def test_CalculateDensitySensitivity(self):
        sensitivity = KratosOA.CollectiveExpression([Kratos.Expression.ElementExpression(self.model_part)])
        self.response_function.CalculateGradient({Kratos.DENSITY: sensitivity})

        # calculate element density sensitivity
        self._CheckSensitivity(
            self.response_function,
            self.model_part.Elements,
            lambda x: x.Properties[KratosOA.DENSITY_SENSITIVITY],
            lambda x, y: self._UpdateProperties(Kratos.DENSITY, x, y),
            sensitivity.GetContainerExpressions()[0].Evaluate(),
            1e-6,
            6)


if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()