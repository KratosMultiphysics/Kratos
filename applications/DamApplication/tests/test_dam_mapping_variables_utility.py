"""Regression tests for the DamApplication nodal result mapping.

The tests exercise the public Python interface (``MappingVariables2DUtilities``
and ``MappingVariables3DUtilities``) and pin the numerical behaviour that must
be preserved when the custom cell search is replaced by the Kratos Core point
locators.
"""

import os
import tempfile

import KratosMultiphysics as KM
import KratosMultiphysics.DamApplication as KratosDam
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.PoromechanicsApplication import NODAL_CAUCHY_STRESS_TENSOR

INITIAL_NODAL_CAUCHY_STRESS_TENSOR = (
    KratosDam.INITIAL_NODAL_CAUCHY_STRESS_TENSOR
)

_RELATIVE_TOLERANCE = 1.0e-12


def _AddMappingVariables(model_part):
    model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(KM.TEMPERATURE)
    model_part.AddNodalSolutionStepVariable(KratosDam.NODAL_REFERENCE_TEMPERATURE)
    model_part.AddNodalSolutionStepVariable(NODAL_CAUCHY_STRESS_TENSOR)
    model_part.AddNodalSolutionStepVariable(
        INITIAL_NODAL_CAUCHY_STRESS_TENSOR
    )


def _Temperature(node_id):
    return 10.0 * node_id + 2.5


def _ReferenceTemperature(node_id):
    return -3.0 * node_id + 7.25


def _DisplacementX(node_id):
    return 1.5 * node_id - 2.0


def _DisplacementY(node_id):
    return -0.75 * node_id + 4.0


def _DisplacementZ(node_id):
    return 0.5 * node_id + 1.0


def _StressComponent(node_id, i, j):
    return 0.5 * node_id + 10.0 * i + j


def _CreateOldModelPart(model, dimension, coordinates, connectivity, element_name):
    model_part = model.CreateModelPart("old", 2)
    _AddMappingVariables(model_part)

    for node_id, coordinate in enumerate(coordinates, start=1):
        model_part.CreateNewNode(node_id, *coordinate)

    properties = model_part.CreateNewProperties(0)
    for element_id, element_connectivity in enumerate(connectivity, start=1):
        model_part.CreateNewElement(
            element_name, element_id, element_connectivity, properties
        )

    for node in model_part.Nodes:
        node_id = node.Id
        stress = KM.Matrix(dimension, dimension)
        for i in range(dimension):
            for j in range(dimension):
                stress[i, j] = _StressComponent(node_id, i, j)

        node.SetSolutionStepValue(KM.TEMPERATURE, _Temperature(node_id))
        node.SetSolutionStepValue(
            KratosDam.NODAL_REFERENCE_TEMPERATURE,
            _ReferenceTemperature(node_id),
        )
        node.SetSolutionStepValue(KM.DISPLACEMENT_X, _DisplacementX(node_id))
        node.SetSolutionStepValue(KM.DISPLACEMENT_Y, _DisplacementY(node_id))
        node.SetSolutionStepValue(KM.DISPLACEMENT_Z, _DisplacementZ(node_id))
        initial_stress = KM.Matrix(dimension, dimension)
        for i in range(dimension):
            for j in range(dimension):
                initial_stress[i, j] = -stress[i, j]

        node.SetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR, stress)
        node.SetSolutionStepValue(
            INITIAL_NODAL_CAUCHY_STRESS_TENSOR, initial_stress
        )

    return model_part


def _CreateNewModelPart(model, coordinates, initial_coordinates=None):
    model_part = model.CreateModelPart("new", 2)
    _AddMappingVariables(model_part)

    for node_id, coordinate in enumerate(coordinates, start=1):
        node = model_part.CreateNewNode(node_id, *coordinate)
        if initial_coordinates is not None:
            node.SetInitialPosition(*initial_coordinates[node_id - 1])

    return model_part


def _CreateUtility(dimension):
    if dimension == 2:
        return KratosDam.MappingVariables2DUtilities()
    return KratosDam.MappingVariables3DUtilities()


def _Interpolate(nodal_values, weights):
    return sum(value * weight for value, weight in zip(nodal_values, weights))


def _ExpectedClamped(nodal_values, weights):
    return min(max(_Interpolate(nodal_values, weights), min(nodal_values)), max(nodal_values))


class TestDamMappingVariablesUtility(KratosUnittest.TestCase):
    def _MapAll(self, dimension, old_model_part, new_model_part, **flags):
        utility = _CreateUtility(dimension)
        utility.MappingMechanicalModelParts(
            old_model_part,
            new_model_part,
            flags.get("add_displacement", True),
            flags.get("add_stress", True),
        )
        utility.MappingThermalModelParts(
            old_model_part,
            new_model_part,
            flags.get("add_temperature", True),
            flags.get("add_reference_temperature", True),
        )

    def _AssertAllVariables(
        self,
        node,
        node_ids,
        dimension,
        weights,
        clamp=False,
    ):
        def expected(values):
            if clamp:
                return _ExpectedClamped(values, weights)
            return _Interpolate(values, weights)

        self.assertAlmostEqual(
            node.GetSolutionStepValue(KM.TEMPERATURE),
            expected([_Temperature(node_id) for node_id in node_ids]),
            delta=_RELATIVE_TOLERANCE,
        )
        self.assertAlmostEqual(
            node.GetSolutionStepValue(KratosDam.NODAL_REFERENCE_TEMPERATURE),
            expected(
                [_ReferenceTemperature(node_id) for node_id in node_ids]
            ),
            delta=_RELATIVE_TOLERANCE,
        )
        self.assertAlmostEqual(
            node.GetSolutionStepValue(KM.DISPLACEMENT_X),
            expected([_DisplacementX(node_id) for node_id in node_ids]),
            delta=_RELATIVE_TOLERANCE,
        )
        self.assertAlmostEqual(
            node.GetSolutionStepValue(KM.DISPLACEMENT_Y),
            expected([_DisplacementY(node_id) for node_id in node_ids]),
            delta=_RELATIVE_TOLERANCE,
        )
        if dimension == 3:
            self.assertAlmostEqual(
                node.GetSolutionStepValue(KM.DISPLACEMENT_Z),
                expected([_DisplacementZ(node_id) for node_id in node_ids]),
                delta=_RELATIVE_TOLERANCE,
            )

        stress = node.GetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR)
        initial_stress = node.GetSolutionStepValue(
            INITIAL_NODAL_CAUCHY_STRESS_TENSOR
        )
        for i in range(dimension):
            for j in range(dimension):
                expected_stress = expected(
                    [
                        _StressComponent(node_id, i, j)
                        for node_id in node_ids
                    ]
                )
                self.assertAlmostEqual(
                    stress[i, j], expected_stress, delta=_RELATIVE_TOLERANCE
                )
                self.assertAlmostEqual(
                    initial_stress[i, j],
                    expected_stress,
                    delta=_RELATIVE_TOLERANCE,
                )

    def _AssertPresetValues(self, node, dimension, offset=0.0):
        self.assertAlmostEqual(
            node.GetSolutionStepValue(KM.TEMPERATURE),
            101.0 + offset,
            delta=_RELATIVE_TOLERANCE,
        )
        self.assertAlmostEqual(
            node.GetSolutionStepValue(KratosDam.NODAL_REFERENCE_TEMPERATURE),
            202.0 + offset,
            delta=_RELATIVE_TOLERANCE,
        )
        self.assertAlmostEqual(
            node.GetSolutionStepValue(KM.DISPLACEMENT_X),
            303.0 + offset,
            delta=_RELATIVE_TOLERANCE,
        )
        self.assertAlmostEqual(
            node.GetSolutionStepValue(KM.DISPLACEMENT_Y),
            404.0 + offset,
            delta=_RELATIVE_TOLERANCE,
        )
        if dimension == 3:
            self.assertAlmostEqual(
                node.GetSolutionStepValue(KM.DISPLACEMENT_Z),
                505.0 + offset,
                delta=_RELATIVE_TOLERANCE,
            )
        stress = node.GetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR)
        initial_stress = node.GetSolutionStepValue(
            INITIAL_NODAL_CAUCHY_STRESS_TENSOR
        )
        for i in range(dimension):
            for j in range(dimension):
                self.assertAlmostEqual(
                    stress[i, j], 606.0 + offset, delta=_RELATIVE_TOLERANCE
                )
                self.assertAlmostEqual(
                    initial_stress[i, j],
                    707.0 + offset,
                    delta=_RELATIVE_TOLERANCE,
                )

    def _PresetNodeValues(self, node, dimension, offset=0.0):
        node.SetSolutionStepValue(KM.TEMPERATURE, 101.0 + offset)
        node.SetSolutionStepValue(
            KratosDam.NODAL_REFERENCE_TEMPERATURE, 202.0 + offset
        )
        node.SetSolutionStepValue(KM.DISPLACEMENT_X, 303.0 + offset)
        node.SetSolutionStepValue(KM.DISPLACEMENT_Y, 404.0 + offset)
        node.SetSolutionStepValue(KM.DISPLACEMENT_Z, 505.0 + offset)
        stress = KM.Matrix(dimension, dimension)
        initial_stress = KM.Matrix(dimension, dimension)
        for i in range(dimension):
            for j in range(dimension):
                stress[i, j] = 606.0 + offset
                initial_stress[i, j] = 707.0 + offset
        node.SetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR, stress)
        node.SetSolutionStepValue(
            INITIAL_NODAL_CAUCHY_STRESS_TENSOR, initial_stress
        )

    def test_2d_triangle_interpolation(self):
        model = KM.Model()
        old = _CreateOldModelPart(
            model,
            2,
            [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)],
            [[1, 2, 3]],
            "Element2D3N",
        )
        new = _CreateNewModelPart(
            model, [(0.2, 0.3, 0.0), (1.0 / 3.0, 1.0 / 3.0, 0.0)]
        )

        self._MapAll(2, old, new)

        weights = [1.0 - 0.2 - 0.3, 0.2, 0.3]
        self._AssertAllVariables(new.GetNode(1), [1, 2, 3], 2, weights)
        self._AssertAllVariables(
            new.GetNode(2), [1, 2, 3], 2, [1.0 / 3.0] * 3
        )

    def test_2d_quadrilateral_interpolation(self):
        model = KM.Model()
        old = _CreateOldModelPart(
            model,
            2,
            [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0)],
            [[1, 2, 3, 4]],
            "Element2D4N",
        )
        new = _CreateNewModelPart(model, [(0.25, 0.75, 0.0)])

        self._MapAll(2, old, new)

        x, y = 0.25, 0.75
        weights = [
            (1.0 - x) * (1.0 - y),
            x * (1.0 - y),
            x * y,
            (1.0 - x) * y,
        ]
        self._AssertAllVariables(new.GetNode(1), [1, 2, 3, 4], 2, weights)

    def test_2d_shared_triangle_boundary(self):
        model = KM.Model()
        old = _CreateOldModelPart(
            model,
            2,
            [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0)],
            [[1, 2, 3], [1, 3, 4]],
            "Element2D3N",
        )
        new = _CreateNewModelPart(model, [(0.5, 0.5, 0.0)])

        self._MapAll(2, old, new)

        # On the shared diagonal only nodes 1 and 3 contribute. Both adjacent
        # elements must therefore produce the same value.
        self._AssertAllVariables(
            new.GetNode(1), [1, 3], 2, [0.5, 0.5], clamp=False
        )

    def test_2d_shared_quadrilateral_edge(self):
        model = KM.Model()
        old = _CreateOldModelPart(
            model,
            2,
            [
                (0.0, 0.0, 0.0),
                (1.0, 0.0, 0.0),
                (2.0, 0.0, 0.0),
                (2.0, 1.0, 0.0),
                (1.0, 1.0, 0.0),
                (0.0, 1.0, 0.0),
            ],
            [[1, 2, 5, 6], [2, 3, 4, 5]],
            "Element2D4N",
        )
        new = _CreateNewModelPart(model, [(1.0, 0.5, 0.0)])

        self._MapAll(2, old, new)

        self._AssertAllVariables(
            new.GetNode(1), [2, 5], 2, [0.5, 0.5], clamp=False
        )

    def test_2d_outside_points_are_clamped_or_skipped(self):
        model = KM.Model()
        old = _CreateOldModelPart(
            model,
            2,
            [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)],
            [[1, 2, 3]],
            "Element2D3N",
        )
        new = _CreateNewModelPart(model, [(-0.1, -0.1, 0.0), (5.0, 5.0, 0.0)])
        self._PresetNodeValues(new.GetNode(2), 2)

        self._MapAll(2, old, new)

        # Just outside the mesh: the adaptive search tolerance locates the
        # element and the shape-function extrapolation is clamped to the
        # element nodal range.
        self._AssertAllVariables(
            new.GetNode(1), [1, 2, 3], 2, [1.2, -0.1, -0.1], clamp=True
        )

        # Far outside the mesh: the node is reported as not found. Thermal
        # values are left untouched, while the stress tensors are zeroed for
        # the mechanical mapping (historical Dam behaviour).
        missing_node = new.GetNode(2)
        self.assertAlmostEqual(
            missing_node.GetSolutionStepValue(KM.TEMPERATURE),
            101.0,
            delta=_RELATIVE_TOLERANCE,
        )
        self.assertAlmostEqual(
            missing_node.GetSolutionStepValue(
                KratosDam.NODAL_REFERENCE_TEMPERATURE
            ),
            202.0,
            delta=_RELATIVE_TOLERANCE,
        )
        self.assertAlmostEqual(
            missing_node.GetSolutionStepValue(KM.DISPLACEMENT_X),
            303.0,
            delta=_RELATIVE_TOLERANCE,
        )
        self.assertAlmostEqual(
            missing_node.GetSolutionStepValue(KM.DISPLACEMENT_Y),
            404.0,
            delta=_RELATIVE_TOLERANCE,
        )
        for tensor in (
            missing_node.GetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR),
            missing_node.GetSolutionStepValue(
                INITIAL_NODAL_CAUCHY_STRESS_TENSOR
            ),
        ):
            self.assertAlmostEqual(tensor[0, 0], 0.0, delta=_RELATIVE_TOLERANCE)
            self.assertAlmostEqual(tensor[0, 1], 0.0, delta=_RELATIVE_TOLERANCE)
            self.assertAlmostEqual(tensor[1, 0], 0.0, delta=_RELATIVE_TOLERANCE)
            self.assertAlmostEqual(tensor[1, 1], 0.0, delta=_RELATIVE_TOLERANCE)

    def test_3d_tetrahedron_interpolation(self):
        model = KM.Model()
        old = _CreateOldModelPart(
            model,
            3,
            [
                (0.0, 0.0, 0.0),
                (1.0, 0.0, 0.0),
                (0.0, 1.0, 0.0),
                (0.0, 0.0, 1.0),
            ],
            [[1, 2, 3, 4]],
            "Element3D4N",
        )
        new = _CreateNewModelPart(model, [(0.25, 0.25, 0.25)])

        self._MapAll(3, old, new)

        self._AssertAllVariables(new.GetNode(1), [1, 2, 3, 4], 3, [0.25] * 4)

    def test_3d_hexahedron_interpolation(self):
        model = KM.Model()
        old = _CreateOldModelPart(
            model,
            3,
            [
                (0.0, 0.0, 0.0),
                (1.0, 0.0, 0.0),
                (1.0, 1.0, 0.0),
                (0.0, 1.0, 0.0),
                (0.0, 0.0, 1.0),
                (1.0, 0.0, 1.0),
                (1.0, 1.0, 1.0),
                (0.0, 1.0, 1.0),
            ],
            [[1, 2, 3, 4, 5, 6, 7, 8]],
            "Element3D8N",
        )
        new = _CreateNewModelPart(model, [(0.25, 0.5, 0.75)])

        self._MapAll(3, old, new)

        x, y, z = 0.25, 0.5, 0.75
        weights = [
            (1.0 - x) * (1.0 - y) * (1.0 - z),
            x * (1.0 - y) * (1.0 - z),
            x * y * (1.0 - z),
            (1.0 - x) * y * (1.0 - z),
            (1.0 - x) * (1.0 - y) * z,
            x * (1.0 - y) * z,
            x * y * z,
            (1.0 - x) * y * z,
        ]
        self._AssertAllVariables(
            new.GetNode(1), [1, 2, 3, 4, 5, 6, 7, 8], 3, weights
        )

    def test_3d_shared_tetrahedron_face(self):
        model = KM.Model()
        old = _CreateOldModelPart(
            model,
            3,
            [
                (0.0, 0.0, 0.0),
                (1.0, 0.0, 0.0),
                (0.0, 1.0, 0.0),
                (0.0, 0.0, 1.0),
                (1.0, 1.0, 1.0),
            ],
            [[1, 2, 3, 4], [2, 3, 4, 5]],
            "Element3D4N",
        )
        new = _CreateNewModelPart(model, [(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0)])

        self._MapAll(3, old, new)

        self._AssertAllVariables(
            new.GetNode(1), [2, 3, 4], 3, [1.0 / 3.0] * 3, clamp=False
        )

    def test_3d_shared_hexahedron_face(self):
        model = KM.Model()
        old = _CreateOldModelPart(
            model,
            3,
            [
                (0.0, 0.0, 0.0),
                (1.0, 0.0, 0.0),
                (1.0, 1.0, 0.0),
                (0.0, 1.0, 0.0),
                (0.0, 0.0, 1.0),
                (1.0, 0.0, 1.0),
                (1.0, 1.0, 1.0),
                (0.0, 1.0, 1.0),
                (2.0, 0.0, 0.0),
                (2.0, 1.0, 0.0),
                (2.0, 0.0, 1.0),
                (2.0, 1.0, 1.0),
            ],
            [
                [1, 2, 3, 4, 5, 6, 7, 8],
                [2, 9, 10, 3, 6, 11, 12, 7],
            ],
            "Element3D8N",
        )
        new = _CreateNewModelPart(model, [(1.0, 0.5, 0.5)])

        self._MapAll(3, old, new)

        self._AssertAllVariables(
            new.GetNode(1), [2, 3, 6, 7], 3, [0.25] * 4, clamp=False
        )

    def test_3d_outside_points_are_clamped_or_skipped(self):
        model = KM.Model()
        old = _CreateOldModelPart(
            model,
            3,
            [
                (0.0, 0.0, 0.0),
                (1.0, 0.0, 0.0),
                (0.0, 1.0, 0.0),
                (0.0, 0.0, 1.0),
            ],
            [[1, 2, 3, 4]],
            "Element3D4N",
        )
        new = _CreateNewModelPart(
            model, [(-0.1, -0.1, -0.1), (5.0, 5.0, 5.0)]
        )
        self._PresetNodeValues(new.GetNode(2), 3)

        self._MapAll(3, old, new)

        self._AssertAllVariables(
            new.GetNode(1),
            [1, 2, 3, 4],
            3,
            [1.3, -0.1, -0.1, -0.1],
            clamp=True,
        )

        missing_node = new.GetNode(2)
        self.assertAlmostEqual(
            missing_node.GetSolutionStepValue(KM.TEMPERATURE),
            101.0,
            delta=_RELATIVE_TOLERANCE,
        )
        self.assertAlmostEqual(
            missing_node.GetSolutionStepValue(
                KratosDam.NODAL_REFERENCE_TEMPERATURE
            ),
            202.0,
            delta=_RELATIVE_TOLERANCE,
        )
        self.assertAlmostEqual(
            missing_node.GetSolutionStepValue(KM.DISPLACEMENT_X),
            303.0,
            delta=_RELATIVE_TOLERANCE,
        )
        self.assertAlmostEqual(
            missing_node.GetSolutionStepValue(KM.DISPLACEMENT_Y),
            404.0,
            delta=_RELATIVE_TOLERANCE,
        )
        self.assertAlmostEqual(
            missing_node.GetSolutionStepValue(KM.DISPLACEMENT_Z),
            505.0,
            delta=_RELATIVE_TOLERANCE,
        )
        for tensor in (
            missing_node.GetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR),
            missing_node.GetSolutionStepValue(
                INITIAL_NODAL_CAUCHY_STRESS_TENSOR
            ),
        ):
            for i in range(3):
                for j in range(3):
                    self.assertAlmostEqual(
                        tensor[i, j], 0.0, delta=_RELATIVE_TOLERANCE
                    )

    def test_variable_selection(self):
        for dimension, coordinates, connectivity, element_name in (
            (
                2,
                [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)],
                [[1, 2, 3]],
                "Element2D3N",
            ),
            (
                3,
                [
                    (0.0, 0.0, 0.0),
                    (1.0, 0.0, 0.0),
                    (0.0, 1.0, 0.0),
                    (0.0, 0.0, 1.0),
                ],
                [[1, 2, 3, 4]],
                "Element3D4N",
            ),
        ):
            with self.subTest(dimension=dimension):
                model = KM.Model()
                old = _CreateOldModelPart(
                    model, dimension, coordinates, connectivity, element_name
                )
                probe_coordinates = list((0.1,) * dimension)
                probe_coordinates += [0.0] * (3 - dimension)
                new = _CreateNewModelPart(
                    model, [tuple(probe_coordinates)]
                )
                node = new.GetNode(1)
                self._PresetNodeValues(node, dimension)

                utility = _CreateUtility(dimension)
                utility.MappingMechanicalModelParts(old, new, False, False)
                utility.MappingThermalModelParts(old, new, False, True)

                self.assertAlmostEqual(
                    node.GetSolutionStepValue(KM.TEMPERATURE),
                    101.0,
                    delta=_RELATIVE_TOLERANCE,
                )
                self.assertAlmostEqual(
                    node.GetSolutionStepValue(KM.DISPLACEMENT_X),
                    303.0,
                    delta=_RELATIVE_TOLERANCE,
                )
                self.assertAlmostEqual(
                    node.GetSolutionStepValue(KM.DISPLACEMENT_Y),
                    404.0,
                    delta=_RELATIVE_TOLERANCE,
                )
                self.assertAlmostEqual(
                    node.GetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR)[0, 0],
                    606.0,
                    delta=_RELATIVE_TOLERANCE,
                )
                # NODAL_REFERENCE_TEMPERATURE was requested and is mapped.
                reference_temperature = node.GetSolutionStepValue(
                    KratosDam.NODAL_REFERENCE_TEMPERATURE
                )
                self.assertNotAlmostEqual(
                    reference_temperature, 202.0, delta=1.0e-6
                )

    def test_current_coordinates_are_used_for_locating_nodes(self):
        model = KM.Model()
        old = _CreateOldModelPart(
            model,
            2,
            [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)],
            [[1, 2, 3]],
            "Element2D3N",
        )
        # The node is inside the mesh in the current configuration, even
        # though its initial position lies far outside it.
        new = _CreateNewModelPart(
            model,
            [(0.2, 0.3, 0.0)],
            initial_coordinates=[(100.0, 100.0, 0.0)],
        )

        self._MapAll(2, old, new)

        self._AssertAllVariables(
            new.GetNode(1), [1, 2, 3], 2, [0.5, 0.2, 0.3]
        )

    def test_add_previous_results_end_to_end(self):
        from KratosMultiphysics.DamApplication.dam_analysis import DamAnalysis

        case = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "joint_elastic_cohesive_2d_normal",
            "joint_elastic_cohesive_2d_normal",
        )
        with open(case + "_parameters.json", "r") as parameter_file:
            parameters = KM.Parameters(parameter_file.read())

        parameters["problem_data"]["problem_name"].SetString(
            "mapping_previous_results"
        )
        parameters["problem_data"]["end_time"].SetDouble(0.0)
        parameters["solver_settings"]["model_import_settings"][
            "input_filename"
        ].SetString(case)

        transfer = parameters["transfer_results_process"]
        transfer["add_previous_results"].SetBool(True)
        transfer["type_of_results"].SetString("Mechanical")
        transfer["add_displacement"].SetBool(True)
        transfer["add_stress"].SetBool(False)
        transfer["file_name_mechanical"].SetString(case + ".mdpa")

        with tempfile.TemporaryDirectory() as working_directory:
            with KratosUnittest.WorkFolderScope(working_directory, __file__):
                model = KM.Model()
                DamAnalysis(model, parameters).Run()

        model_part = model["MainModelPart"]
        self.assertEqual(model_part.NumberOfNodes(), 24)

    def test_python_mapping_wrapper(self):
        model = KM.Model()
        old = _CreateOldModelPart(
            model,
            2,
            [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)],
            [[1, 2, 3]],
            "Element2D3N",
        )
        new = _CreateNewModelPart(model, [(0.2, 0.3, 0.0)])

        from KratosMultiphysics.DamApplication import (
            dam_mapping_variables_utility,
        )

        wrapper = dam_mapping_variables_utility.MappingVariablesUtility(2)
        returned_model_part = (
            wrapper.AddPreviousModelPartThermoMechanical(
                new, old, old, True, True, True, True
            )
        )
        self.assertIs(returned_model_part, new)
        self._AssertAllVariables(
            new.GetNode(1), [1, 2, 3], 2, [0.5, 0.2, 0.3]
        )


if __name__ == "__main__":
    KratosUnittest.main()
