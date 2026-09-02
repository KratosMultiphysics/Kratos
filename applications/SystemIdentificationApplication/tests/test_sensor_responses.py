# --- Kratos Imports ---
import KratosMultiphysics
import KratosMultiphysics.SystemIdentificationApplication
from KratosMultiphysics.KratosUnittest import TestCase, main, skipIf

has_opt_app: bool
try:
    import KratosMultiphysics.OptimizationApplication
    has_opt_app = True
except ImportError as _:
    has_opt_app = False

has_structural_app: bool
try:
    import KratosMultiphysics.StructuralMechanicsApplication
    has_structural_app = True
except ImportError as _:
    has_structural_app = False

# --- STD Imports ---
import math


class TestSensorResponses(TestCase):

    @skipIf(not has_opt_app or not has_structural_app, "Requires the Optimization and StructuralMechanics applications")
    def test_DisplacementSensorResponse(self) -> None:
        """
             element 1
            4 ------ 5 ------ 6
            |        |        |
            |   s1   |   s2   |
            |        |        |
            1 ------ 2 ------ 3
                      element 2

            Displacement sensors s1 and s2 at the center
            of element 1 and element 2 respectively.
        """
        # Construct model part.
        model: KratosMultiphysics.Model = KratosMultiphysics.Model()
        model_part: KratosMultiphysics.ModelPart = model.CreateModelPart("model_part")

        # Initialize model part.
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)

        # Construct nodes.
        nodes: list[KratosMultiphysics.Node] = [
            model_part.CreateNewNode(1, 1, 1, 0),
            model_part.CreateNewNode(2, 2, 1, 0),
            model_part.CreateNewNode(3, 3, 1, 0),
            model_part.CreateNewNode(4, 1, 2, 0),
            model_part.CreateNewNode(5, 2, 2, 0),
            model_part.CreateNewNode(6, 3, 2, 0)]

        # Fill nodes' initial state.
        for node in nodes:
            node.SetSolutionStepValue(
                KratosMultiphysics.DISPLACEMENT,
                10 * node.GetInitialPosition())

        # Construct elements.
        properties: KratosMultiphysics.Properties = KratosMultiphysics.Properties(1)
        properties[KratosMultiphysics.YOUNG_MODULUS] = 2.1e9
        properties[KratosMultiphysics.POISSON_RATIO] = 0.3
        properties[KratosMultiphysics.CONSTITUTIVE_LAW] = KratosMultiphysics.StructuralMechanicsApplication.LinearElasticPlaneStress2DLaw()
        elements: list[KratosMultiphysics.Element] = [
            model_part.CreateNewElement("SmallDisplacementElement2D4N", 1, [1, 2, 5, 4], properties),
            model_part.CreateNewElement("SmallDisplacementElement2D4N", 2, [2, 3, 6, 5], properties)]

        # Construct sensors.
        prototype: KratosMultiphysics.SystemIdentificationApplication.DisplacementSensorResponse = KratosMultiphysics.SystemIdentificationApplication.DisplacementSensorResponse()
        sensors: list[KratosMultiphysics.SystemIdentificationApplication.DisplacementSensorResponse] = [
            prototype.Create(
                model_part,
                model_part,
                7,
                KratosMultiphysics.Parameters(R"""{
                    "name" : "A",
                    "position" : [1.5, 1.5, 0.0],
                    "direction" : [1, 2, 0],
                    "design_variables" : ["YOUNG_MODULUS", "THICKNESS"]
                }""")),
            prototype.Create(
                model_part,
                model_part,
                8,
                KratosMultiphysics.Parameters(R"""{
                    "name" : "B",
                    "position" : [2.5, 1.5, 0.0],
                    "direction" : [-2, 1, 0],
                    "design_variables" : ["SHAPE_X", "SHAPE_Y", "SHAPE_Z"]
                }"""))]

        # Let the sensors precompute whatever they need to.
        for sensor in sensors:
            sensor.ComputeCache(model_part)

        # Check whether the sensors are correctly located.
        self.assertEqual(
            sensors[0].GetContainingElement().Id,
            1)
        self.assertEqual(
            sensors[1].GetContainingElement().Id,
            2)

        # Check computed values.
        def ProjectDisplacements(
            direction: list[float],
            geometry: KratosMultiphysics.Geometry) -> float:
                interpolated: KratosMultiphysics.Array3 = KratosMultiphysics.Array3([0, 0, 0])
                for node in geometry:
                    interpolated += node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
                interpolated /= 4
                output: float = 0.0
                norm: float = 0.0
                for direction_component, displacement_component in zip(direction, interpolated):
                    output += direction_component * displacement_component
                    norm += direction_component * direction_component
                return output / math.sqrt(norm)

        self.assertAlmostEqual(
            sensors[0].ComputeInnerValue(
                elements[0],
                model_part.ProcessInfo,
                0),
            ProjectDisplacements(
                [1.0, 2.0 , 0.0],
                elements[0].GetGeometry()))

        self.assertAlmostEqual(
            sensors[0].ComputeInnerValue(
                elements[1],
                model_part.ProcessInfo,
                0),
            0.0)

        self.assertAlmostEqual(
            sensors[1].ComputeInnerValue(
                elements[0],
                model_part.ProcessInfo,
                0),
            0.0)

        self.assertAlmostEqual(
            sensors[1].ComputeInnerValue(
                elements[1],
                model_part.ProcessInfo,
                0),
            ProjectDisplacements(
                [-2.0, 1.0 , 0.0],
                elements[1].GetGeometry()))

        # Collect sensor 0 state variables.
        state_variables: list[KratosMultiphysics.IAdjoint.DynamicVariable] = []
        for component in (KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.DISPLACEMENT_Z):
            for i_node in range(4):
                state_variables.append(KratosMultiphysics.IAdjoint.DynamicVariable(
                    component,
                    i_node))

        # Check sensor 0 state derivatives.
        state_variables = sensors[0].GetStateVariables(
            elements[0],
            model_part.ProcessInfo)
        state_variables = state_variables[::-1] # <== check a different order than what the sensor provides
        for variable in state_variables:
            self.assertIn(variable, state_variables)

        state_derivatives = sensors[0].ComputeInnerDerivative(
            elements[0],
            state_variables,
            model_part.ProcessInfo,
            0)
        for variable, derivative in zip(state_variables, state_derivatives):
            reference: float = 0.25 * (1.0 if variable == KratosMultiphysics.DISPLACEMENT_X else 2.0) / math.sqrt(5)
            self.assertAlmostEqual(
                derivative,
                reference)

        # Check sensor 0 design variables.
        design_variables = sensors[0].GetDesignVariables(
            elements[0],
            model_part.ProcessInfo)
        self.assertEqual(
            len(design_variables),
            1)
        self.assertEqual(
            design_variables[0],
            KratosMultiphysics.YOUNG_MODULUS)

        # Check sensor 0 design derivatives.
        design_derivatives = sensors[0].ComputeInnerDerivative(
            elements[0],
            design_variables,
            model_part.ProcessInfo,
            0)
        self.assertEqual(
            len(design_derivatives),
            1)
        self.assertAlmostEqual(
            design_derivatives[0],
            0.0)

        # Collect sensor 0 state variables.
        state_variables: list[KratosMultiphysics.IAdjoint.DynamicVariable] = []
        for component in (KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.DISPLACEMENT_Z):
            for i_node in range(4):
                state_variables.append(KratosMultiphysics.IAdjoint.DynamicVariable(
                    component,
                    i_node))

        # Check sensor 0 state derivatives.
        state_variables = sensors[1].GetStateVariables(
            elements[1],
            model_part.ProcessInfo)
        state_variables = state_variables[::-1] # <== check a different order than what the sensor provides
        for variable in state_variables:
            self.assertIn(variable, state_variables)

        state_derivatives = sensors[1].ComputeInnerDerivative(
            elements[1],
            state_variables,
            model_part.ProcessInfo,
            0)
        for variable, derivative in zip(state_variables, state_derivatives):
            reference: float = 0.25 * (-2.0 if variable == KratosMultiphysics.DISPLACEMENT_X else 1.0) / math.sqrt(5)
            self.assertAlmostEqual(
                derivative,
                reference)

        # Check sensor 0 design variables.
        design_variables = sensors[1].GetDesignVariables(
            elements[1],
            model_part.ProcessInfo)
        self.assertEqual(
            len(design_variables),
            4 * 3)
        for i_node in range(4):
            for variable_type in (KratosMultiphysics.SHAPE_X, KratosMultiphysics.SHAPE_Y, KratosMultiphysics.SHAPE_Z):
                self.assertIn(
                    KratosMultiphysics.IAdjoint.DynamicVariable(variable_type, i_node),
                    design_variables)

        # Check sensor 0 design derivatives.
        design_derivatives = sensors[1].ComputeInnerDerivative(
            elements[1],
            design_variables,
            model_part.ProcessInfo,
            0)
        self.assertEqual(
            len(design_derivatives),
            4 * 3)
        self.assertVectorAlmostEqual(
            design_derivatives,
            [0.0 for _ in range(4 * 3)])

        # Clear sensors' caches.
        for sensor in sensors:
            sensor.ClearCache()


    @skipIf(not has_opt_app or not has_structural_app, "Requires the Optimization and StructuralMechanics applications")
    def test_SensorAggregateResponse(self) -> None:
        """
                element 1
            4 ------ 5 ------ 6
            |        |        |
            |   s1   |   s2   |
            |        |        |
            1 ------ 2 ------ 3
                        element 2

            Displacement sensors s1 and s2 at the center
            of element 1 and element 2 respectively.
            Response aggregates both sensors.
        """
        # Construct model part.
        model: KratosMultiphysics.Model = KratosMultiphysics.Model()
        model_part: KratosMultiphysics.ModelPart = model.CreateModelPart("model_part")

        # Initialize model part.
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)

        # Construct nodes.
        nodes: list[KratosMultiphysics.Node] = [
            model_part.CreateNewNode(1, 1, 1, 0),
            model_part.CreateNewNode(2, 2, 1, 0),
            model_part.CreateNewNode(3, 3, 1, 0),
            model_part.CreateNewNode(4, 1, 2, 0),
            model_part.CreateNewNode(5, 2, 2, 0),
            model_part.CreateNewNode(6, 3, 2, 0)]

        # Fill nodes' initial state.
        for node in nodes:
            node.SetSolutionStepValue(
                KratosMultiphysics.DISPLACEMENT,
                10 * node.GetInitialPosition())

        # Construct elements.
        properties: KratosMultiphysics.Properties = KratosMultiphysics.Properties(1)
        properties[KratosMultiphysics.YOUNG_MODULUS] = 2.1e9
        properties[KratosMultiphysics.POISSON_RATIO] = 0.3
        properties[KratosMultiphysics.CONSTITUTIVE_LAW] = KratosMultiphysics.StructuralMechanicsApplication.LinearElasticPlaneStress2DLaw()
        elements: list[KratosMultiphysics.Element] = [
            model_part.CreateNewElement("SmallDisplacementElement2D4N", 1, [1, 2, 5, 4], properties),
            model_part.CreateNewElement("SmallDisplacementElement2D4N", 2, [2, 3, 6, 5], properties)]

        # Construct sensors.
        prototype: KratosMultiphysics.SystemIdentificationApplication.DisplacementSensorResponse = KratosMultiphysics.SystemIdentificationApplication.DisplacementSensorResponse()
        sensors: list[KratosMultiphysics.SystemIdentificationApplication.DisplacementSensorResponse] = [
            prototype.Create(
                model_part,
                model_part,
                7,
                KratosMultiphysics.Parameters(R"""{
                    "name" : "A",
                    "position" : [1.5, 1.5, 0.0],
                    "direction" : [1, 2, 0],
                    "design_variables" : ["YOUNG_MODULUS", "THICKNESS"]
                }""")),
            prototype.Create(
                model_part,
                model_part,
                8,
                KratosMultiphysics.Parameters(R"""{
                    "name" : "B",
                    "position" : [2.5, 1.5, 0.0],
                    "direction" : [-2, 1, 0],
                    "design_variables" : ["SHAPE_X", "SHAPE_Y", "SHAPE_Z"]
                }"""))]

        # Set sensor measured values.
        for i_sensor, sensor in enumerate(sensors, 1):
            sensor.GetNode().SetValue(
                KratosMultiphysics.SystemIdentificationApplication.SENSOR_MEASURED_VALUE,
                float(i_sensor))
            sensor.GetNode().SetValue(
                KratosMultiphysics.SystemIdentificationApplication.SENSOR_NORMALIZATION_FACTOR,
                100 * float(i_sensor))

        # Construct response.
        exponent: int = 4
        response = KratosMultiphysics.SystemIdentificationApplication.SensorAggregateResponse(
            exponent)
        response.AddSensors(sensors)
        response.SetDesignVariableTypes([
            KratosMultiphysics.YOUNG_MODULUS,
            KratosMultiphysics.THICKNESS,
            KratosMultiphysics.SHAPE_X,
            KratosMultiphysics.SHAPE_Y,
            KratosMultiphysics.SHAPE_Z])

        # Let the response precompute whatever it needs.
        response.ComputeCache(model_part)

        # Check computed values.
        def ProjectDisplacements(
            direction: list[float],
            geometry: KratosMultiphysics.Geometry) -> float:
                interpolated: KratosMultiphysics.Array3 = KratosMultiphysics.Array3([0, 0, 0])
                for node in geometry:
                    interpolated += node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
                interpolated /= 4
                output: float = 0.0
                norm: float = 0.0
                for direction_component, displacement_component in zip(direction, interpolated):
                    output += direction_component * displacement_component
                    norm += direction_component * direction_component
                return output / math.sqrt(norm)

        self.assertAlmostEqual(
            response.ComputeInnerValue(
                elements[0],
                model_part.ProcessInfo,
                0),
            ((ProjectDisplacements(
                [1.0, 2.0, 0.0],
                elements[0].GetGeometry()) - 1) * 100) ** exponent,
            places = 0)

        self.assertAlmostEqual(
            response.ComputeInnerValue(
                elements[1],
                model_part.ProcessInfo,
                0),
            ((ProjectDisplacements(
                [-2.0, 1.0, 0.0],
                elements[1].GetGeometry()) - 2) * 200) ** exponent,
            places = 0)

        # Check state derivatives for element 0.
        reference_state_variables: list[KratosMultiphysics.IAdjoint.DynamicVariable] = []
        for component in (KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.DISPLACEMENT_Z):
            for i_node in range(4):
                reference_state_variables.append(KratosMultiphysics.IAdjoint.DynamicVariable(
                    component,
                    i_node))

        state_variables: list[KratosMultiphysics.IAdjoint.DynamicVariable] = response.GetStateVariables(
            elements[1],
            model_part.ProcessInfo)
        state_variables = state_variables[::-1] # <== check a different order than what the sensor provides
        for variable in state_variables:
            self.assertIn(variable, reference_state_variables)

        state_derivatives = response.ComputeInnerDerivative(
            elements[0],
            state_variables,
            model_part.ProcessInfo,
            0)
        for variable, derivative in zip(state_variables, state_derivatives):
            sensor_value: float = ProjectDisplacements(
                [1.0, 2.0, 0.0],
                elements[0].GetGeometry())
            sensor_derivative: float = 0.25 * (1.0 if variable == KratosMultiphysics.DISPLACEMENT_X else 2.0) / math.sqrt(5)
            reference: float = exponent * ((sensor_value - 1) ** (exponent - 1)) * (100 ** exponent) * sensor_derivative
            self.assertAlmostEqual(
                derivative,
                reference)

        # Check state derivatives for element 1.
        state_variables = response.GetStateVariables(
            elements[1],
            model_part.ProcessInfo)
        state_variables = state_variables[::-1] # <== check a different order than what the sensor provides
        for variable in state_variables:
            self.assertIn(variable, reference_state_variables)

        state_derivatives = response.ComputeInnerDerivative(
            elements[1],
            state_variables,
            model_part.ProcessInfo,
            0)
        for variable, derivative in zip(state_variables, state_derivatives):
            sensor_value: float = ProjectDisplacements(
                [-2.0, 1.0, 0.0],
                elements[1].GetGeometry())
            sensor_derivative: float = 0.25 * (-2.0 if variable == KratosMultiphysics.DISPLACEMENT_X else 1.0) / math.sqrt(5)
            reference: float = exponent * ((sensor_value - 2) ** (exponent - 1)) * (200 ** exponent) * sensor_derivative
            self.assertAlmostEqual(
                derivative,
                reference,
                places = 2)

        # Check design variables.
        for element in elements:
            design_variables = response.GetDesignVariables(
                element,
                model_part.ProcessInfo)
            self.assertEqual(
                len(design_variables),
                1 + 4 * 3)
            self.assertIn(
                KratosMultiphysics.YOUNG_MODULUS,
                design_variables)
            for i_node in range(4):
                for variable_type in (KratosMultiphysics.SHAPE_X, KratosMultiphysics.SHAPE_Y, KratosMultiphysics.SHAPE_Z):
                    self.assertIn(
                        KratosMultiphysics.IAdjoint.DynamicVariable(variable_type, i_node),
                        design_variables)

        # Check design derivatives.
        for element in elements:
            design_derivatives = response.ComputeInnerDerivative(
                element,
                design_variables,
                model_part.ProcessInfo,
                0)
            self.assertEqual(
                len(design_derivatives),
                1 + 4 * 3)
            self.assertVectorAlmostEqual(
                design_derivatives,
                [0.0 for _ in range(1 + 4 * 3)])

        # Clear sensors' caches.
        for sensor in sensors:
            sensor.ClearCache()



if __name__ == "__main__":
    main()
