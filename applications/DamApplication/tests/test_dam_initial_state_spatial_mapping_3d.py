"""Spatial mapping of Dam initial state between nonmatching 3D meshes."""
import json
import math
import sys
import tempfile
from pathlib import Path

import KratosMultiphysics as KM
import KratosMultiphysics.DamApplication as Dam
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.ConvectionDiffusionApplication as CDA

Dam.__path__.insert(0, str(Path(__file__).resolve().parents[1] / "python_scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from KratosMultiphysics.DamApplication import dam_initial_state_file
from KratosMultiphysics.DamApplication import dam_initial_state_mapping
from KratosMultiphysics.DamApplication.export_dam_final_state_process import (
    ExportDamFinalStateProcess)
from KratosMultiphysics.DamApplication.import_dam_initial_state_process import (
    ImportDamInitialStateProcess)
from test_dam_initial_state_spatial_mapping import (
    _assign_analytic_stresses, _build_part, _build_two_material_part,
    _displacements, _import, _integrations, _make_strategy, _max_error,
    _shift_mesh, _stresses, _write_state)


# ---------------------------------------------------------------------------
# Mesh generators. Node ordering follows the Kratos geometry definitions.
# ---------------------------------------------------------------------------

def _mesh(nodes, elements, element, geometry_type):
    return {
        "nodes": nodes,
        "elements": elements,
        "element": element,
        "geometry_type": geometry_type,
        "dimension": 3,
    }


def _tet4_mesh():
    return _mesh(
        {1: [0.0, 0.0, 0.0], 2: [1.0, 0.0, 0.0],
         3: [0.0, 1.0, 0.0], 4: [0.0, 0.0, 1.0]},
        [(1, [1, 2, 3, 4])],
        "SmallDisplacementElement3D4N", 15)


def _small_tet4_mesh():
    return _mesh(
        {1: [0.1, 0.1, 0.1], 2: [0.5, 0.1, 0.1],
         3: [0.1, 0.5, 0.1], 4: [0.1, 0.1, 0.5]},
        [(1, [1, 2, 3, 4])],
        "SmallDisplacementElement3D4N", 15)


def _tet10_mesh(curved=False):
    nodes = {
        1: [0.0, 0.0, 0.0], 2: [1.0, 0.0, 0.0],
        3: [0.0, 1.0, 0.0], 4: [0.0, 0.0, 1.0],
        5: [0.5, 0.0, 0.0], 6: [0.5, 0.5, 0.0],
        7: [0.0, 0.5, 0.0], 8: [0.0, 0.0, 0.5],
        9: [0.5, 0.0, 0.5], 10: [0.0, 0.5, 0.5],
    }
    if curved:
        nodes[5] = [0.5, -0.06, -0.04]
        nodes[6] = [0.56, 0.54, -0.05]
        nodes[7] = [-0.05, 0.52, -0.03]
        nodes[8] = [-0.04, -0.05, 0.53]
        nodes[9] = [0.55, -0.04, 0.48]
        nodes[10] = [-0.03, 0.54, 0.51]
    return _mesh(nodes, [(1, list(range(1, 11)))],
                 "SmallDisplacementElement3D10N", 14)


def _prism6_mesh():
    return _mesh(
        {1: [0.0, 0.0, 0.0], 2: [1.0, 0.0, 0.0],
         3: [0.0, 1.0, 0.0], 4: [0.0, 0.0, 1.0],
         5: [1.0, 0.0, 1.0], 6: [0.0, 1.0, 1.0]},
        [(1, [1, 2, 3, 4, 5, 6])],
        "SmallDisplacementElement3D6N", 5)


def _prism15_mesh(curved=False):
    nodes = {
        1: [0.0, 0.0, 0.0], 2: [1.0, 0.0, 0.0],
        3: [0.0, 1.0, 0.0], 4: [0.0, 0.0, 1.0],
        5: [1.0, 0.0, 1.0], 6: [0.0, 1.0, 1.0],
        7: [0.5, 0.0, 0.0], 8: [0.5, 0.5, 0.0],
        9: [0.0, 0.5, 0.0], 10: [0.0, 0.0, 0.5],
        11: [1.0, 0.0, 0.5], 12: [0.0, 1.0, 0.5],
        13: [0.5, 0.0, 1.0], 14: [0.5, 0.5, 1.0],
        15: [0.0, 0.5, 1.0],
    }
    if curved:
        nodes[10] = [-0.04, -0.03, 0.5]
        nodes[13] = [0.5, -0.05, 1.04]
        nodes[14] = [0.54, 0.53, 1.03]
    return _mesh(nodes, [(1, list(range(1, 16)))],
                 "SmallDisplacementElement3D15N", 4)


def _hex20_mesh(curved=False):
    nodes = {
        1: [0.0, 0.0, 0.0], 2: [1.0, 0.0, 0.0],
        3: [1.0, 1.0, 0.0], 4: [0.0, 1.0, 0.0],
        5: [0.0, 0.0, 1.0], 6: [1.0, 0.0, 1.0],
        7: [1.0, 1.0, 1.0], 8: [0.0, 1.0, 1.0],
        9: [0.5, 0.0, 0.0], 10: [1.0, 0.5, 0.0],
        11: [0.5, 1.0, 0.0], 12: [0.0, 0.5, 0.0],
        13: [0.0, 0.0, 0.5], 14: [1.0, 0.0, 0.5],
        15: [1.0, 1.0, 0.5], 16: [0.0, 1.0, 0.5],
        17: [0.5, 0.0, 1.0], 18: [1.0, 0.5, 1.0],
        19: [0.5, 1.0, 1.0], 20: [0.0, 0.5, 1.0],
    }
    if curved:
        for node_id in range(9, 21):
            x, y, z = nodes[node_id]
            nodes[node_id] = [
                x + 0.04 * math.sin(3.1 * y + 2.0 * z),
                y + 0.04 * math.sin(2.3 * x + 1.7 * z),
                z + 0.04 * math.sin(1.9 * x + 2.7 * y)]
    return _mesh(nodes, [(1, list(range(1, 21)))],
                 "SmallDisplacementElement3D20N", 1)


def _hex27_mesh(curved=False):
    mesh = _hex20_mesh(curved)
    nodes = mesh["nodes"]
    nodes[21] = [0.5, 0.5, 0.0]
    nodes[22] = [0.5, 0.0, 0.5]
    nodes[23] = [1.0, 0.5, 0.5]
    nodes[24] = [0.5, 1.0, 0.5]
    nodes[25] = [0.0, 0.5, 0.5]
    nodes[26] = [0.5, 0.5, 1.0]
    nodes[27] = [0.5, 0.5, 0.5]
    return _mesh(nodes, [(1, list(range(1, 28)))],
                 "SmallDisplacementElement3D27N", 2)


def _hex8_mesh(nx, ny, nz, lx=1.0, ly=2.0, lz=1.0, curved=False):
    nodes = {}
    index = {}
    node_id = 1
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                x = lx * i / nx
                y = ly * j / ny
                z = lz * k / nz
                if curved and 0 < i < nx and 0 < j < ny and 0 < k < nz:
                    x += 0.03 * math.sin(3.0 * y + 2.0 * z)
                    y += 0.03 * math.sin(2.0 * x + 3.0 * z)
                    z += 0.03 * math.sin(3.0 * x + y)
                nodes[node_id] = [x, y, z]
                index[(i, j, k)] = node_id
                node_id += 1
    elements = []
    element_id = 1
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                connectivity = [
                    index[(i, j, k)], index[(i + 1, j, k)],
                    index[(i + 1, j + 1, k)], index[(i, j + 1, k)],
                    index[(i, j, k + 1)], index[(i + 1, j, k + 1)],
                    index[(i + 1, j + 1, k + 1)],
                    index[(i, j + 1, k + 1)]]
                elements.append((element_id, connectivity))
                element_id += 1
    return _mesh(nodes, elements, "SmallDisplacementElement3D8N", 3)


def _tet4_block_mesh(nx, ny, nz, lx=1.0, ly=2.0, lz=1.0):
    mesh = _hex8_mesh(nx, ny, nz, lx, ly, lz)
    index = {}
    node_id = 1
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                index[(i, j, k)] = node_id
                node_id += 1
    elements = []
    element_id = 1
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                b0 = index[(i, j, k)]
                b1 = index[(i + 1, j, k)]
                b2 = index[(i + 1, j + 1, k)]
                b3 = index[(i, j + 1, k)]
                t0 = index[(i, j, k + 1)]
                t1 = index[(i + 1, j, k + 1)]
                t2 = index[(i + 1, j + 1, k + 1)]
                t3 = index[(i, j + 1, k + 1)]
                for connectivity in (
                        [b0, b1, b2, t2], [b0, b2, b3, t2],
                        [b0, b3, t3, t2], [b0, t3, t0, t2],
                        [b0, t0, t1, t2], [b0, t1, b1, t2]):
                    elements.append((element_id, connectivity))
                    element_id += 1
    return _mesh(mesh["nodes"], elements, "SmallDisplacementElement3D4N", 15)


# ---------------------------------------------------------------------------
# Analytic fields.
# ---------------------------------------------------------------------------

def _constant_3d(value):
    return lambda x, y, z: list(value)


def _linear_3d(x, y, z):
    return [100.0 + 10.0 * x + 20.0 * y + 5.0 * z,
            -50.0 + 5.0 * x - 7.0 * y + 2.0 * z,
            30.0 + 7.0 * x - 3.0 * y + 4.0 * z,
            11.0 + x + 2.0 * y - 3.0 * z,
            -7.0 + 2.0 * x - y + 4.0 * z,
            5.0 - 3.0 * x + 2.0 * y + z]


def _quadratic_3d(x, y, z):
    return [
        100.0 + 10.0 * x + 20.0 * y + 5.0 * z + 2.0 * x * x +
        3.0 * y * y + 1.5 * z * z + x * y - 2.0 * x * z + y * z,
        -50.0 + 5.0 * x - 7.0 * y + 2.0 * z + x * x - 2.0 * y * y +
        0.5 * z * z + 2.0 * x * y - x * z,
        30.0 + 7.0 * x - 3.0 * y + 4.0 * z + 1.5 * x * x + y * y -
        z * z - x * y + 2.0 * y * z,
        11.0 + x + 2.0 * y - 3.0 * z + 0.5 * x * x - y * y + z * z + x * z,
        -7.0 + 2.0 * x - y + 4.0 * z + x * y + y * z,
        5.0 - 3.0 * x + 2.0 * y + z - x * z + 2.0 * y * z]


def _temperature_3d(x, y, z):
    return 10.0 + 2.0 * x - 3.0 * y + 4.0 * z


def _reference_temperature_3d(x, y, z):
    return 5.0 + x + y - 2.0 * z


class TestDamInitialStateSpatialMapping3D(KratosUnittest.TestCase):

    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.directory.cleanup()

    def _state_file(self, name):
        return Path(self.directory.name) / name

    def _roundtrip(self, source_mesh, destination_mesh, function,
                   source_law=None, destination_law=None,
                   source_order=None, destination_order=None,
                   tolerance=1.e-8):
        if source_law is None:
            source_law = SMA.LinearElastic3DLaw()
        if destination_law is None:
            destination_law = SMA.LinearElastic3DLaw()
        source_model = KM.Model()
        source = _build_part(
            source_model, "source", source_mesh, source_law,
            integration_order=source_order)
        _assign_analytic_stresses(source, function)
        path = self._state_file("state.json")
        _write_state(path, source, stresses=True)

        destination_model = KM.Model()
        destination = _build_part(
            destination_model, "destination", destination_mesh,
            destination_law, integration_order=destination_order)
        _import(destination_model, "destination", path, stresses=True)
        actual = _stresses(destination)
        expected = {
            element.Id: [list(function(*coordinates))
                         for coordinates in _integrations(destination)[element.Id]]
            for element in destination.Elements}
        self.assertLessEqual(_max_error(actual, expected), tolerance)
        for node in destination.Nodes:
            for value in node.GetSolutionStepValue(KM.DISPLACEMENT):
                self.assertEqual(value, 0.0)
        return destination

    # -- geometry and reconstruction ---------------------------------------

    def test_constant_stress_t4(self):
        self._roundtrip(
            _tet4_mesh(), _small_tet4_mesh(),
            _constant_3d([100.0, -50.0, 30.0, 11.0, -7.0, 5.0]))

    def test_t4_uses_constant_reconstruction_only(self):
        terms, _ = dam_initial_state_mapping._build_reconstruction_operator(
            [[0.25, 0.25, 0.25, 0.0]], 3)
        self.assertEqual(len(terms), 1)
        terms, _ = dam_initial_state_mapping._build_reconstruction_operator(
            [[0.138, 0.138, 0.138, 0.0], [0.585, 0.138, 0.138, 0.0],
             [0.138, 0.585, 0.138, 0.0], [0.138, 0.138, 0.585, 0.0]], 3)
        self.assertEqual(len(terms), 4)

    def test_linear_stress_t10(self):
        self._roundtrip(
            _tet10_mesh(), _small_tet4_mesh(),
            _constant_3d([100.0, -50.0, 30.0, 11.0, -7.0, 5.0]))
        self._roundtrip(_tet10_mesh(), _small_tet4_mesh(), _linear_3d)

    def test_constant_and_linear_stress_h8(self):
        self._roundtrip(
            _hex8_mesh(1, 2, 1), _hex8_mesh(2, 4, 2),
            _constant_3d([100.0, -50.0, 30.0, 11.0, -7.0, 5.0]))
        self._roundtrip(
            _hex8_mesh(1, 2, 1), _hex8_mesh(2, 4, 2), _linear_3d)

    def test_higher_order_stress_h20_and_h27(self):
        self._roundtrip(
            _hex20_mesh(), _hex8_mesh(2, 2, 2, ly=1.0), _quadratic_3d,
            source_order=3, tolerance=1.e-8)
        self._roundtrip(
            _hex27_mesh(), _hex8_mesh(2, 2, 2, ly=1.0), _quadratic_3d,
            source_order=3, tolerance=1.e-8)

    def test_prism_p6_and_p15(self):
        self._roundtrip(
            _prism6_mesh(), _small_tet4_mesh(), _linear_3d,
            source_order=2)
        self._roundtrip(
            _prism15_mesh(), _small_tet4_mesh(), _linear_3d,
            source_order=2)

    def test_mixed_source_destination_geometry(self):
        self._roundtrip(
            _hex20_mesh(), _hex8_mesh(2, 2, 2, ly=1.0), _quadratic_3d,
            source_order=3)
        self._roundtrip(
            _prism15_mesh(), _small_tet4_mesh(), _linear_3d,
            source_order=2)
        self._roundtrip(
            _hex8_mesh(1, 2, 1), _tet4_block_mesh(1, 2, 1), _linear_3d)

    def test_different_integration_rules(self):
        # One-point source supports constant fields only.
        self._roundtrip(
            _hex8_mesh(1, 2, 1), _hex8_mesh(2, 4, 2),
            _constant_3d([100.0, -50.0, 30.0, 11.0, -7.0, 5.0]),
            source_order=1, destination_order=2)
        # Over-integrated source still reproduces a linear field.
        self._roundtrip(
            _hex8_mesh(1, 2, 1), _hex8_mesh(2, 4, 2), _linear_3d,
            source_order=5, destination_order=1, tolerance=1.e-8)
        # GI_GAUSS_4 is also available for the 3D solid geometries.
        self._roundtrip(
            _hex8_mesh(1, 2, 1), _hex8_mesh(2, 4, 2), _linear_3d,
            source_order=4, destination_order=1, tolerance=1.e-8)

    def test_coarse_to_fine_h8(self):
        self._roundtrip(
            _hex8_mesh(1, 2, 1), _hex8_mesh(2, 4, 2), _linear_3d)

    def test_fine_to_coarse_h8(self):
        self._roundtrip(
            _hex8_mesh(2, 4, 2), _hex8_mesh(1, 2, 1), _linear_3d)

    def test_distorted_h8_linear(self):
        self._roundtrip(
            _hex8_mesh(2, 4, 2, curved=True), _hex8_mesh(3, 6, 3),
            _linear_3d, tolerance=1.e-8)

    def test_curved_higher_order_elements(self):
        self._roundtrip(
            _tet10_mesh(curved=True), _small_tet4_mesh(), _linear_3d,
            source_order=5, tolerance=1.e-8)
        self._roundtrip(
            _hex20_mesh(curved=True), _hex8_mesh(2, 2, 2, ly=1.0), _linear_3d,
            source_order=3, tolerance=1.e-8)

    def test_boundary_point_acceptance(self):
        source_model = KM.Model()
        source = _build_part(
            source_model, "source", _hex8_mesh(1, 2, 1),
            SMA.LinearElastic3DLaw(), temperature=20.0)
        path = self._state_file("boundary.json")
        _write_state(path, source, stresses=False, temperature=True)
        destination_mesh = _hex8_mesh(2, 4, 2)
        shifted = {
            node_id: [x + 2.e-10, y + 2.e-10, z + 2.e-10]
            for node_id, (x, y, z) in destination_mesh["nodes"].items()}
        destination_model = KM.Model()
        destination = _build_part(
            destination_model, "destination",
            {**destination_mesh, "nodes": shifted},
            SMA.LinearElastic3DLaw())
        _import(destination_model, "destination", path, stresses=False,
                temperature=True)
        for node in destination.Nodes:
            self.assertAlmostEqual(
                node.GetSolutionStepValue(KM.TEMPERATURE), 20.0, places=10)

    def test_outside_domain_rejection(self):
        source_model = KM.Model()
        source = _build_part(
            source_model, "source", _hex8_mesh(1, 2, 1),
            SMA.LinearElastic3DLaw(), temperature=20.0)
        path = self._state_file("boundary.json")
        _write_state(path, source, stresses=False, temperature=True)
        destination_mesh = _hex8_mesh(2, 4, 2)
        shifted = {
            node_id: [x + 1.e-6, y, z]
            for node_id, (x, y, z) in destination_mesh["nodes"].items()}
        destination_model = KM.Model()
        _build_part(
            destination_model, "destination",
            {**destination_mesh, "nodes": shifted},
            SMA.LinearElastic3DLaw())
        with self.assertRaisesRegex(RuntimeError, "outside the source domain"):
            _import(destination_model, "destination", path, stresses=False,
                    temperature=True)

    # -- thermal fields -----------------------------------------------------

    def _temperature_case(self, source_mesh, destination_mesh,
                          temperature_function, reference_function=None,
                          import_reference=True):
        source_model = KM.Model()
        source = _build_part(
            source_model, "source", source_mesh, SMA.LinearElastic3DLaw())
        for node in source.Nodes:
            node.SetSolutionStepValue(
                KM.TEMPERATURE, temperature_function(node.X0, node.Y0, node.Z0))
            if reference_function is not None:
                node.SetSolutionStepValue(
                    Dam.NODAL_REFERENCE_TEMPERATURE,
                    reference_function(node.X0, node.Y0, node.Z0))
        path = self._state_file("temperature.json")
        _write_state(
            path, source, stresses=False, temperature=True,
            reference=import_reference)

        destination_model = KM.Model()
        destination = _build_part(
            destination_model, "destination", destination_mesh,
            SMA.LinearElastic3DLaw(), temperature=0.0,
            reference_temperature=0.0)
        _import(
            destination_model, "destination", path, stresses=False,
            temperature=True, reference=import_reference)
        for node in destination.Nodes:
            self.assertAlmostEqual(
                node.GetSolutionStepValue(KM.TEMPERATURE),
                temperature_function(node.X0, node.Y0, node.Z0), places=10)
            if reference_function is not None:
                self.assertAlmostEqual(
                    node.GetSolutionStepValue(
                        Dam.NODAL_REFERENCE_TEMPERATURE),
                    reference_function(node.X0, node.Y0, node.Z0), places=10)
        return destination

    def test_uniform_temperature_3d(self):
        self._temperature_case(
            _hex8_mesh(1, 2, 1), _hex8_mesh(2, 4, 2),
            lambda x, y, z: 25.0, import_reference=False)

    def test_linear_temperature_3d(self):
        self._temperature_case(
            _hex8_mesh(1, 2, 1), _hex8_mesh(2, 4, 2), _temperature_3d,
            import_reference=False)

    def test_independent_reference_temperature_3d(self):
        self._temperature_case(
            _hex8_mesh(1, 2, 1), _hex8_mesh(2, 4, 2), _temperature_3d,
            reference_function=_reference_temperature_3d)

    def test_thermal_initial_strain_compensation_3d(self):
        value = [100.0, -50.0, 30.0, 11.0, -7.0, 5.0]
        source_model = KM.Model()
        source = _build_part(
            source_model, "source", _hex8_mesh(1, 2, 1),
            Dam.ThermalLinearElastic3DLaw(), thermal_expansion=1.e-5,
            temperature=30.0, reference_temperature=20.0)
        _assign_analytic_stresses(
            source, _constant_3d(value), compensate=True)
        path = self._state_file("thermal.json")
        _write_state(
            path, source, stresses=True, temperature=True, reference=True)

        destination_model = KM.Model()
        destination = _build_part(
            destination_model, "destination", _hex8_mesh(2, 4, 2),
            Dam.ThermalLinearElastic3DLaw(), thermal_expansion=1.e-5,
            temperature=30.0, reference_temperature=20.0)
        _import(
            destination_model, "destination", path, stresses=True,
            temperature=True, reference=True)
        for element in destination.Elements:
            for point in element.CalculateOnIntegrationPoints(
                    KM.CAUCHY_STRESS_VECTOR, destination.ProcessInfo):
                for actual, expected in zip(point, value):
                    self.assertAlmostEqual(actual, expected, places=8)
        for node in destination.Nodes:
            for displacement in node.GetSolutionStepValue(KM.DISPLACEMENT):
                self.assertEqual(displacement, 0.0)
            self.assertAlmostEqual(
                node.GetSolutionStepValue(KM.TEMPERATURE), 30.0, places=10)
            self.assertAlmostEqual(
                node.GetSolutionStepValue(
                    Dam.NODAL_REFERENCE_TEMPERATURE), 20.0, places=10)

    def test_conflicting_shared_node_temperature_rejected(self):
        source_a_model = KM.Model()
        source_a = _build_part(
            source_a_model, "source_a", _hex8_mesh(1, 2, 1),
            SMA.LinearElastic3DLaw(), temperature=20.0)
        path_a = self._state_file("a.json")
        _write_state(path_a, source_a, stresses=False, temperature=True)

        source_b_model = KM.Model()
        source_b = _build_part(
            source_b_model, "source_b", _hex8_mesh(1, 2, 1),
            SMA.LinearElastic3DLaw(), temperature=30.0)
        path_b = self._state_file("b.json")
        _write_state(path_b, source_b, stresses=False, temperature=True)

        destination_model = KM.Model()
        destination = _build_part(
            destination_model, "destination", _hex8_mesh(1, 2, 1),
            SMA.LinearElastic3DLaw())
        processes = []
        for path in (path_a, path_b):
            processes.append(ImportDamInitialStateProcess(
                destination_model, KM.Parameters(json.dumps({
                    "model_part_name": "destination",
                    "input_file_name": str(path),
                    "import_stresses": False,
                    "import_temperature": True,
                    "transfer_method": "spatial_mapping"}))))
        for process in processes:
            process.PrevalidateInitialState()
        with self.assertRaisesRegex(RuntimeError, "conflict"):
            dam_initial_state_file.ValidateCombinedImports(processes)

    # -- materials and constitutive behaviour ------------------------------

    def test_matching_mesh_3d(self):
        value = [100.0, -50.0, 30.0, 11.0, -7.0, 5.0]
        mesh = _hex8_mesh(2, 4, 2)
        source_model = KM.Model()
        source = _build_part(
            source_model, "source", mesh, SMA.LinearElastic3DLaw())
        _assign_analytic_stresses(source, _constant_3d(value))
        path = self._state_file("matching.json")
        _write_state(path, source, stresses=True)

        destination_model = KM.Model()
        destination = _build_part(
            destination_model, "destination", mesh, SMA.LinearElastic3DLaw())
        process = ImportDamInitialStateProcess(
            destination_model, KM.Parameters(json.dumps({
                "model_part_name": "destination",
                "input_file_name": str(path),
                "import_stresses": True,
                "transfer_method": "matching_mesh"})))
        process.PrevalidateInitialState()
        process.ExecuteBeforeSolutionLoop()
        for element in destination.Elements:
            for point in element.CalculateOnIntegrationPoints(
                    KM.CAUCHY_STRESS_VECTOR, destination.ProcessInfo):
                for actual, expected in zip(point, value):
                    self.assertAlmostEqual(actual, expected, places=10)
        for node in destination.Nodes:
            for displacement in node.GetSolutionStepValue(KM.DISPLACEMENT):
                self.assertEqual(displacement, 0.0)

    def test_six_component_stress_ordering(self):
        value = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
        destination = self._roundtrip(
            _hex8_mesh(1, 2, 1), _hex8_mesh(2, 4, 2),
            _constant_3d(value))
        for element in destination.Elements:
            for point in element.CalculateOnIntegrationPoints(
                    KM.CAUCHY_STRESS_VECTOR, destination.ProcessInfo):
                self.assertEqual(len(point), 6)
                for actual, expected in zip(point, value):
                    self.assertAlmostEqual(actual, expected, places=10)

    def test_exact_constitutive_initial_stress_3d(self):
        value = [100.0, -50.0, 30.0, 11.0, -7.0, 5.0]
        self._roundtrip(
            _hex8_mesh(1, 2, 1), _hex8_mesh(2, 4, 2),
            _constant_3d(value), tolerance=1.e-10)

    def test_zero_initial_displacement_3d(self):
        destination = self._roundtrip(
            _hex8_mesh(1, 2, 1), _hex8_mesh(2, 4, 2), _linear_3d)
        for node in destination.Nodes:
            for displacement in node.GetSolutionStepValue(KM.DISPLACEMENT):
                self.assertEqual(displacement, 0.0)

    def test_material_incompatibility_rejected_3d(self):
        source_model = KM.Model()
        source = _build_part(
            source_model, "source", _hex8_mesh(1, 2, 1),
            SMA.LinearElastic3DLaw(), young_modulus=1.e7,
            poisson_ratio=0.25)
        _assign_analytic_stresses(
            source, _constant_3d([100.0, -50.0, 30.0, 11.0, -7.0, 5.0]))
        path = self._state_file("incompatible.json")
        _write_state(path, source, stresses=True)

        destination_model = KM.Model()
        _build_part(
            destination_model, "destination", _hex8_mesh(2, 4, 2),
            SMA.LinearElastic3DLaw(), young_modulus=2.e7,
            poisson_ratio=0.25)
        with self.assertRaisesRegex(RuntimeError, "not present in the source"):
            _import(destination_model, "destination", path, stresses=True)

    def test_two_material_3d(self):
        concrete = [100.0, -50.0, 30.0, 11.0, -7.0, 5.0]
        soil = [-20.0, 40.0, 10.0, 3.0, 2.0, -1.0]
        concrete_mesh = _hex8_mesh(1, 2, 1, lx=0.5)
        soil_mesh = _shift_mesh(_hex8_mesh(1, 2, 1, lx=0.5), 0.5)

        concrete_model = KM.Model()
        concrete_part = _build_part(
            concrete_model, "Concrete", concrete_mesh,
            SMA.LinearElastic3DLaw(), young_modulus=1.e7, poisson_ratio=0.2)
        _assign_analytic_stresses(
            concrete_part, _constant_3d(concrete))
        concrete_path = self._state_file("concrete.json")
        _write_state(concrete_path, concrete_part, stresses=True)

        soil_model = KM.Model()
        soil_part = _build_part(
            soil_model, "Soil", soil_mesh, SMA.LinearElastic3DLaw(),
            young_modulus=5.e6, poisson_ratio=0.25)
        _assign_analytic_stresses(soil_part, _constant_3d(soil))
        soil_path = self._state_file("soil.json")
        _write_state(soil_path, soil_part, stresses=True)

        destination_mesh = _hex8_mesh(2, 4, 2)
        destination_model = KM.Model()
        _, destination_parts = _build_two_material_part(
            destination_model, "destination", destination_mesh,
            [("Parts_Concrete", lambda x, y, z: x < 0.5,
              {"law": SMA.LinearElastic3DLaw(),
               "young_modulus": 1.e7, "poisson_ratio": 0.2}),
             ("Parts_Soil", lambda x, y, z: x >= 0.5,
              {"law": SMA.LinearElastic3DLaw(),
               "young_modulus": 5.e6, "poisson_ratio": 0.25})])
        processes = []
        for subpart_name, path in (
                ("Parts_Concrete", concrete_path),
                ("Parts_Soil", soil_path)):
            processes.append(ImportDamInitialStateProcess(
                destination_model, KM.Parameters(json.dumps({
                    "model_part_name": "destination." + subpart_name,
                    "input_file_name": str(path),
                    "import_stresses": True,
                    "transfer_method": "spatial_mapping"}))))
        for process in processes:
            process.PrevalidateInitialState()
        dam_initial_state_file.ValidateCombinedImports(processes)
        for process in processes:
            process.ExecuteBeforeSolutionLoop()

        process_info = destination_model["destination"].ProcessInfo
        for name, expected in (("Parts_Concrete", concrete),
                               ("Parts_Soil", soil)):
            for element in destination_parts[name].Elements:
                for point in element.CalculateOnIntegrationPoints(
                        KM.CAUCHY_STRESS_VECTOR, process_info):
                    for actual, value in zip(point, expected):
                        self.assertAlmostEqual(actual, value, places=8)

    # -- equilibrium diagnostic --------------------------------------------

    def test_mechanical_condition_contributions_3d(self):
        model = KM.Model()
        part = _build_part(
            model, "part", _hex8_mesh(1, 2, 1),
            SMA.LinearElastic3DLaw())
        loaded_node = max(part.Nodes, key=lambda node: node.Y)
        condition = part.CreateNewCondition(
            "PointLoadCondition3D1N", 1, [loaded_node.Id],
            part.GetProperties(1))
        condition.SetValue(SMA.POINT_LOAD, [3.0, 7.0, 11.0])
        diagnostic = dam_initial_state_mapping.EvaluateMechanicalResidual(part)
        self.assertEqual(
            len(diagnostic["residual_by_node"][loaded_node.Id]), 3)
        self.assertAlmostEqual(
            diagnostic["residual_by_node"][loaded_node.Id][0], 3.0,
            places=10)
        self.assertAlmostEqual(
            diagnostic["residual_by_node"][loaded_node.Id][1], 7.0,
            places=10)
        self.assertAlmostEqual(
            diagnostic["residual_by_node"][loaded_node.Id][2], 11.0,
            places=10)
        self.assertAlmostEqual(
            diagnostic["max_free_residual"], 11.0, places=10)

    def test_surface_load_condition_3d(self):
        model = KM.Model()
        part = _build_part(
            model, "part", _hex8_mesh(1, 1, 1, ly=1.0),
            SMA.LinearElastic3DLaw())
        def node_at(x, y, z):
            return next(
                node.Id for node in part.Nodes
                if abs(node.X0 - x) < 1.e-14 and abs(node.Y0 - y) < 1.e-14
                and abs(node.Z0 - z) < 1.e-14)

        top_face = [
            node_at(0.0, 0.0, 1.0), node_at(1.0, 0.0, 1.0),
            node_at(1.0, 1.0, 1.0), node_at(0.0, 1.0, 1.0)]
        self.assertEqual(len(set(top_face)), 4)
        condition = part.CreateNewCondition(
            "SurfaceLoadCondition3D4N", 1, top_face, part.GetProperties(1))
        condition.SetValue(SMA.SURFACE_LOAD, [0.0, 0.0, -5.0])
        diagnostic = dam_initial_state_mapping.EvaluateMechanicalResidual(part)
        total = [0.0, 0.0, 0.0]
        for node_id in top_face:
            for component in range(3):
                total[component] += diagnostic["residual_by_node"][node_id][
                    component]
        self.assertAlmostEqual(total[0], 0.0, places=10)
        self.assertAlmostEqual(total[1], 0.0, places=10)
        self.assertAlmostEqual(total[2], -5.0, places=10)

    def test_thermal_conditions_excluded_3d(self):
        model = KM.Model()
        part = _build_part(
            model, "part", _hex8_mesh(1, 1, 1),
            SMA.LinearElastic3DLaw(),
            extra_variables=(KM.FACE_HEAT_FLUX,))
        for node in part.Nodes:
            node.AddDof(KM.TEMPERATURE)
        settings = KM.ConvectionDiffusionSettings()
        settings.SetUnknownVariable(KM.TEMPERATURE)
        settings.SetSurfaceSourceVariable(KM.FACE_HEAT_FLUX)
        part.ProcessInfo.SetValue(KM.CONVECTION_DIFFUSION_SETTINGS, settings)
        properties = part.GetProperties(1)
        properties[CDA.AMBIENT_TEMPERATURE] = 100.0
        properties[CDA.CONVECTION_COEFFICIENT] = 10.0
        properties[KM.EMISSIVITY] = 0.0
        top_nodes = sorted(
            node.Id for node in part.Nodes if abs(node.Z - 1.0) < 1.e-14)
        thermal = part.CreateNewCondition(
            "ThermalFace3D3N", 1, top_nodes[:3], properties)
        thermal_rhs = KM.Vector()
        thermal.CalculateRightHandSide(thermal_rhs, part.ProcessInfo)
        self.assertGreater(max(abs(value) for value in thermal_rhs), 1.0)

        condition = part.CreateNewCondition(
            "PointLoadCondition3D1N", 2, [top_nodes[0]],
            part.GetProperties(1))
        condition.SetValue(SMA.POINT_LOAD, [0.0, -5.0, 0.0])
        diagnostic = dam_initial_state_mapping.EvaluateMechanicalResidual(part)
        self.assertAlmostEqual(
            diagnostic["max_free_residual"], 5.0, places=10)
        self.assertAlmostEqual(
            diagnostic["residual_by_node"][top_nodes[0]][1], -5.0,
            places=10)

    def test_free_and_fixed_dofs_3d(self):
        model = KM.Model()
        part = _build_part(
            model, "part", _hex8_mesh(1, 1, 1),
            SMA.LinearElastic3DLaw())
        node = max(part.Nodes, key=lambda item: item.Id)
        node.Fix(KM.DISPLACEMENT_X)
        node.Fix(KM.DISPLACEMENT_Z)
        condition = part.CreateNewCondition(
            "PointLoadCondition3D1N", 1, [node.Id], part.GetProperties(1))
        condition.SetValue(SMA.POINT_LOAD, [3.0, 7.0, 11.0])
        diagnostic = dam_initial_state_mapping.EvaluateMechanicalResidual(part)
        self.assertEqual(len(diagnostic["reaction_sum"]), 3)
        self.assertEqual(
            len(diagnostic["free_residual"]),
            3 * part.NumberOfNodes() - 2)
        self.assertAlmostEqual(
            diagnostic["max_free_residual"], 7.0, places=10)
        self.assertAlmostEqual(
            diagnostic["residual_by_node"][node.Id][0], 3.0, places=10)
        self.assertAlmostEqual(
            diagnostic["residual_by_node"][node.Id][1], 7.0, places=10)
        self.assertAlmostEqual(
            diagnostic["residual_by_node"][node.Id][2], 11.0, places=10)

    def test_mechanical_dof_row_filter_3d_mixed_order(self):
        model = KM.Model()
        part = _build_part(
            model, "part", _hex8_mesh(1, 1, 1, ly=1.0),
            SMA.LinearElastic3DLaw())
        node = next(iter(part.Nodes))
        node.AddDof(KM.TEMPERATURE)

        class _Entity:
            def __init__(self, dofs):
                self.dofs = dofs

            def GetDofList(self, process_info):
                return self.dofs

        entity = _Entity([
            node.GetDof(KM.DISPLACEMENT_Z),
            node.GetDof(KM.TEMPERATURE),
            node.GetDof(KM.DISPLACEMENT_X),
            node.GetDof(KM.DISPLACEMENT_Y)])
        rows, number_of_dofs = \
            dam_initial_state_mapping._MechanicalDofRows(
                entity, part.ProcessInfo, 3)
        self.assertEqual(number_of_dofs, 4)
        self.assertEqual(
            rows, [(0, node.Id, 2), (2, node.Id, 0), (3, node.Id, 1)])
        rows, _ = dam_initial_state_mapping._MechanicalDofRows(
            entity, part.ProcessInfo, 2)
        self.assertEqual(rows, [(2, node.Id, 0), (3, node.Id, 1)])

    def test_residual_is_read_only_3d(self):
        model = KM.Model()
        part = _build_part(
            model, "part", _hex8_mesh(2, 4, 2),
            SMA.LinearElastic3DLaw())
        _assign_analytic_stresses(
            part, _constant_3d([100.0, -50.0, 30.0, 11.0, -7.0, 5.0]))
        before_stress = _stresses(part)
        before_displacement = _displacements(part)
        before_temperature = [
            node.GetSolutionStepValue(KM.TEMPERATURE) for node in part.Nodes]
        dam_initial_state_mapping.EvaluateMechanicalResidual(part)
        self.assertEqual(_stresses(part), before_stress)
        self.assertEqual(_displacements(part), before_displacement)
        self.assertEqual(
            [node.GetSolutionStepValue(KM.TEMPERATURE) for node in part.Nodes],
            before_temperature)

    def test_shared_interface_node_assembly_3d(self):
        model = KM.Model()
        part, subparts = _build_two_material_part(
            model, "part", _hex8_mesh(2, 2, 2),
            [("Parts_A", lambda x, y, z: x < 0.5,
              {"law": SMA.LinearElastic3DLaw(),
               "young_modulus": 1.e7, "poisson_ratio": 0.2}),
             ("Parts_B", lambda x, y, z: x >= 0.5,
              {"law": SMA.LinearElastic3DLaw(),
               "young_modulus": 5.e6, "poisson_ratio": 0.25})])
        interface_nodes = sorted(
            node.Id for node in part.Nodes if abs(node.X0 - 0.5) < 1.e-14)
        self.assertTrue(interface_nodes)
        for node_id in interface_nodes:
            self.assertIn(node_id, [
                node.Id for node in subparts["Parts_A"].Nodes])
            self.assertIn(node_id, [
                node.Id for node in subparts["Parts_B"].Nodes])
        diagnostic = dam_initial_state_mapping.EvaluateMechanicalResidual(part)
        self.assertEqual(
            set(diagnostic["residual_by_node"]),
            {node.Id for node in part.Nodes})

    # -- end-to-end 3D validation ------------------------------------------

    def _gravity_case(self, source_mesh, destination_mesh):
        body = [0.0, -9.81, 0.0]
        density = 2000.0
        weight = density * 9.81 * 1.0 * 2.0 * 1.0

        source_model = KM.Model()
        source = _build_part(
            source_model, "source", source_mesh, SMA.LinearElastic3DLaw(),
            body=body, fixed_bottom=True, density=density)
        _make_strategy(source).Solve()
        path = self._state_file("gravity.json")
        ExportDamFinalStateProcess(
            source_model, KM.Parameters(json.dumps({
                "model_part_name": "source",
                "export_stresses": True,
                "export_temperature": True,
                "export_reference_temperature": True,
                "output_file_name": str(path),
                "echo_level": 0}))).ExecuteFinalize()

        reference_model = KM.Model()
        reference = _build_part(
            reference_model, "reference", destination_mesh,
            SMA.LinearElastic3DLaw(), body=body, fixed_bottom=True,
            density=density)
        _make_strategy(reference).Solve()
        reference_stress = _stresses(reference)
        reference_displacement = math.sqrt(sum(
            value * value for value in _displacements(reference)))

        destination_model = KM.Model()
        destination = _build_part(
            destination_model, "destination", destination_mesh,
            SMA.LinearElastic3DLaw(), body=body, fixed_bottom=True,
            density=density)
        destination_strategy = _make_strategy(destination)
        _import(destination_model, "destination", path, stresses=True,
                temperature=True, reference=True)

        for node in destination.Nodes:
            for value in node.GetSolutionStepValue(KM.DISPLACEMENT):
                self.assertEqual(value, 0.0)

        initial = dam_initial_state_mapping.EvaluateMechanicalResidual(
            destination, force_scale=weight)
        destination_strategy.Solve()
        correction = math.sqrt(sum(
            value * value for value in _displacements(destination)))
        final = dam_initial_state_mapping.EvaluateMechanicalResidual(
            destination, force_scale=weight)
        stress_error = _max_error(
            _stresses(destination), reference_stress)
        scale = max(
            abs(value) for points in reference_stress.values()
            for point in points for value in point)
        return {
            "initial_residual": initial["free_residual_norm"],
            "normalized_initial_residual":
                initial["normalized_free_residual_norm"],
            "correction": correction,
            "reference_displacement": reference_displacement,
            "final_residual": final["free_residual_norm"],
            "stress_relative": stress_error / max(scale, 1.0),
            "reaction_y": final["reaction_sum"][1],
            "weight": weight,
        }

    def test_gravity_equilibrium_3d(self):
        metrics = self._gravity_case(
            _hex8_mesh(2, 4, 2), _hex8_mesh(3, 6, 3))
        self.assertLess(metrics["normalized_initial_residual"], 0.1)
        self.assertLess(
            metrics["correction"], 0.25 * metrics["reference_displacement"])
        self.assertLess(metrics["final_residual"], 1.e-6 * metrics["weight"])
        self.assertLess(metrics["stress_relative"], 0.5)
        self.assertAlmostEqual(
            metrics["reaction_y"], metrics["weight"], places=6)

    def test_gravity_coarsening_3d(self):
        metrics = self._gravity_case(
            _hex8_mesh(3, 6, 3), _hex8_mesh(2, 4, 2))
        self.assertLess(metrics["normalized_initial_residual"], 0.1)
        self.assertLess(
            metrics["correction"], 0.25 * metrics["reference_displacement"])
        self.assertLess(metrics["final_residual"], 1.e-6 * metrics["weight"])
        self.assertLess(metrics["stress_relative"], 0.5)
        self.assertAlmostEqual(
            metrics["reaction_y"], metrics["weight"], places=6)

    def test_thermomechanical_stage_c_3d(self):
        body = [0.0, -9.81, 0.0]
        density = 2000.0
        weight = density * 9.81 * 2.0
        source_model = KM.Model()
        source = _build_part(
            source_model, "source", _hex8_mesh(2, 4, 2),
            Dam.ThermalLinearElastic3DLaw(), body=body,
            fixed_bottom=True, density=density, thermal_expansion=1.e-5,
            temperature=30.0, reference_temperature=20.0)
        _make_strategy(source).Solve()
        path = self._state_file("thermal_gravity.json")
        ExportDamFinalStateProcess(
            source_model, KM.Parameters(json.dumps({
                "model_part_name": "source",
                "export_stresses": True,
                "export_temperature": True,
                "export_reference_temperature": True,
                "output_file_name": str(path),
                "echo_level": 0}))).ExecuteFinalize()

        reference_model = KM.Model()
        reference = _build_part(
            reference_model, "reference", _hex8_mesh(3, 6, 3),
            Dam.ThermalLinearElastic3DLaw(), body=body,
            fixed_bottom=True, density=density, thermal_expansion=1.e-5,
            temperature=30.0, reference_temperature=20.0)
        _make_strategy(reference).Solve()
        reference_stress = _stresses(reference)
        reference_displacement = math.sqrt(sum(
            value * value for value in _displacements(reference)))

        destination_model = KM.Model()
        destination = _build_part(
            destination_model, "destination", _hex8_mesh(3, 6, 3),
            Dam.ThermalLinearElastic3DLaw(), body=body,
            fixed_bottom=True, density=density, thermal_expansion=1.e-5,
            temperature=30.0, reference_temperature=20.0)
        destination_strategy = _make_strategy(destination)
        _import(destination_model, "destination", path, stresses=True,
                temperature=True, reference=True)

        for node in destination.Nodes:
            self.assertAlmostEqual(
                node.GetSolutionStepValue(KM.TEMPERATURE), 30.0, places=10)
            self.assertAlmostEqual(
                node.GetSolutionStepValue(
                    Dam.NODAL_REFERENCE_TEMPERATURE), 20.0, places=10)
            for value in node.GetSolutionStepValue(KM.DISPLACEMENT):
                self.assertEqual(value, 0.0)

        initial = dam_initial_state_mapping.EvaluateMechanicalResidual(
            destination)
        destination_strategy.Solve()
        correction = math.sqrt(sum(
            value * value for value in _displacements(destination)))
        final = dam_initial_state_mapping.EvaluateMechanicalResidual(
            destination)
        scale = max(
            abs(value) for points in reference_stress.values()
            for point in points for value in point)
        stress_error = _max_error(
            _stresses(destination), reference_stress)

        self.assertLess(initial["free_residual_norm"], 0.1 * weight)
        self.assertLess(correction, 0.25 * reference_displacement)
        self.assertLess(final["free_residual_norm"], 1.e-6 * weight)
        self.assertLess(stress_error / max(scale, 1.0), 0.5)
        self.assertAlmostEqual(final["reaction_sum"][1], weight, places=6)

    def test_sequential_abc_thermomechanical_3d(self):
        """Genuine A -> B -> C: B imports A, is heated, and C imports B."""
        law = Dam.ThermalLinearElastic3DLaw()
        body = [0.0, -9.81, 0.0]
        weight = 2000.0 * 9.81 * 2.0

        def build(model, name, mesh, temperature, reference_temperature):
            return _build_part(
                model, name, mesh, law, body=body, fixed_bottom=True,
                density=2000.0, young_modulus=1.e7, poisson_ratio=0.25,
                thermal_expansion=1.e-5, temperature=temperature,
                reference_temperature=reference_temperature)

        def export(part, path):
            ExportDamFinalStateProcess(
                part.GetModel(), KM.Parameters(json.dumps({
                    "model_part_name": part.Name,
                    "export_stresses": True,
                    "export_temperature": True,
                    "export_reference_temperature": True,
                    "output_file_name": str(path),
                    "echo_level": 0}))).ExecuteFinalize()

        def make_process(model, name, path):
            return ImportDamInitialStateProcess(
                model, KM.Parameters(json.dumps({
                    "model_part_name": name,
                    "input_file_name": str(path),
                    "import_stresses": True,
                    "import_temperature": True,
                    "import_reference_temperature": True,
                    "reference_temperature_policy": "preserve",
                    "transfer_method": "spatial_mapping"})))

        # Stage A: gravity on the source mesh at T = T_ref.
        a_model = KM.Model()
        a_part = build(a_model, "A", _hex8_mesh(2, 4, 2), 20.0, 20.0)
        _make_strategy(a_part).Solve()
        path_a = self._state_file("A.json")
        export(a_part, path_a)

        # Stage B: import A into a different mesh and solve gravity.
        b_model = KM.Model()
        b_part = build(b_model, "B", _hex8_mesh(3, 6, 3), 20.0, 20.0)
        b_strategy = _make_strategy(b_part)
        b_process = make_process(b_model, "B", path_a)
        b_process.PrevalidateInitialState()
        b_process.ExecuteBeforeSolutionLoop()
        for node in b_part.Nodes:
            for displacement in node.GetSolutionStepValue(KM.DISPLACEMENT):
                self.assertEqual(displacement, 0.0)
        b_strategy.Solve()

        # Stage B thermal phase: apply T = 30, T_ref = 20 and re-solve.
        for node in b_part.Nodes:
            node.SetSolutionStepValue(KM.TEMPERATURE, 30.0)
            node.SetSolutionStepValue(Dam.NODAL_REFERENCE_TEMPERATURE, 20.0)
        b_strategy.Solve()
        self.assertGreater(
            math.sqrt(sum(value * value for value in _displacements(b_part))),
            0.0)
        path_b = self._state_file("B.json")
        export(b_part, path_b)

        # Stage C: import the B thermal state into another mesh.
        c_model = KM.Model()
        c_part = build(c_model, "C", _hex8_mesh(2, 4, 2), 30.0, 20.0)
        c_strategy = _make_strategy(c_part)
        c_process = make_process(c_model, "C", path_b)
        c_process.PrevalidateInitialState()
        c_process.ExecuteBeforeSolutionLoop()
        for node in c_part.Nodes:
            for displacement in node.GetSolutionStepValue(KM.DISPLACEMENT):
                self.assertEqual(displacement, 0.0)
            self.assertAlmostEqual(
                node.GetSolutionStepValue(KM.TEMPERATURE), 30.0, places=10)
            self.assertAlmostEqual(
                node.GetSolutionStepValue(
                    Dam.NODAL_REFERENCE_TEMPERATURE), 20.0, places=10)
        initial = dam_initial_state_mapping.EvaluateMechanicalResidual(
            c_part, force_scale=weight)
        c_strategy.Solve()
        final = dam_initial_state_mapping.EvaluateMechanicalResidual(
            c_part, force_scale=weight)

        reference_model = KM.Model()
        reference = build(
            reference_model, "reference", _hex8_mesh(2, 4, 2), 30.0, 20.0)
        _make_strategy(reference).Solve()
        reference_stress = _stresses(reference)
        scale = max(
            abs(value) for points in reference_stress.values()
            for point in points for value in point)

        self.assertLess(initial["free_residual_norm"], 0.2 * weight)
        self.assertLess(final["free_residual_norm"], 1.e-6 * weight)
        self.assertAlmostEqual(final["reaction_sum"][1], weight, places=6)
        self.assertLessEqual(
            _max_error(_stresses(c_part), reference_stress), 0.25 * scale)

    def test_multipart_gravity_interface_3d(self):
        body = [0.0, -9.81, 0.0]
        concrete_density = 2400.0
        soil_density = 1800.0
        weight = 9.81 * (concrete_density + soil_density)

        def material_blocks():
            return [
                ("Parts_Concrete", lambda x, y, z: x < 0.5,
                 {"law": SMA.LinearElastic3DLaw(),
                  "young_modulus": 1.e7, "poisson_ratio": 0.2,
                  "density": concrete_density}),
                ("Parts_Soil", lambda x, y, z: x >= 0.5,
                 {"law": SMA.LinearElastic3DLaw(),
                  "young_modulus": 5.e6, "poisson_ratio": 0.25,
                  "density": soil_density})]

        # The source must be mechanically connected at the interface, exactly
        # like the intended destination model. Two disconnected block models
        # would describe a different physical problem at the interface.
        source_model = KM.Model()
        source, source_parts = _build_two_material_part(
            source_model, "source", _hex8_mesh(2, 2, 1),
            material_blocks(), body=body, fixed_bottom=True,
            temperature=30.0, reference_temperature=20.0)
        _make_strategy(source).Solve()
        source_paths = {}
        for name in ("Parts_Concrete", "Parts_Soil"):
            path = self._state_file(name + ".json")
            ExportDamFinalStateProcess(
                source_model, KM.Parameters(json.dumps({
                    "model_part_name": "source." + name,
                    "export_stresses": True,
                    "export_temperature": True,
                    "export_reference_temperature": True,
                    "output_file_name": str(path),
                    "echo_level": 0}))).ExecuteFinalize()
            source_paths[name] = path

        destination_model = KM.Model()
        destination, destination_parts = _build_two_material_part(
            destination_model, "destination", _hex8_mesh(2, 4, 2),
            material_blocks(), body=body, fixed_bottom=True,
            temperature=30.0, reference_temperature=20.0)
        strategy = _make_strategy(destination)
        processes = []
        for name in ("Parts_Concrete", "Parts_Soil"):
            processes.append(ImportDamInitialStateProcess(
                destination_model, KM.Parameters(json.dumps({
                    "model_part_name": "destination." + name,
                    "input_file_name": str(source_paths[name]),
                    "import_stresses": True,
                    "import_temperature": True,
                    "import_reference_temperature": True,
                    "transfer_method": "spatial_mapping"}))))
        for process in processes:
            process.PrevalidateInitialState()
        dam_initial_state_file.ValidateCombinedImports(processes)
        for process in processes:
            process.ExecuteBeforeSolutionLoop()

        reference_model = KM.Model()
        reference, _ = _build_two_material_part(
            reference_model, "reference", _hex8_mesh(2, 4, 2),
            material_blocks(), body=body, fixed_bottom=True,
            temperature=30.0, reference_temperature=20.0)
        _make_strategy(reference).Solve()
        reference_stress = _stresses(reference)

        interface_nodes = [
            node for node in destination.Nodes if abs(node.X0 - 0.5) < 1.e-14]
        self.assertTrue(interface_nodes)
        for node in interface_nodes:
            self.assertAlmostEqual(
                node.GetSolutionStepValue(KM.TEMPERATURE), 30.0, places=10)
            self.assertAlmostEqual(
                node.GetSolutionStepValue(
                    Dam.NODAL_REFERENCE_TEMPERATURE), 20.0, places=10)

        initial = dam_initial_state_mapping.EvaluateMechanicalResidual(
            destination, force_scale=weight)
        strategy.Solve()
        final = dam_initial_state_mapping.EvaluateMechanicalResidual(
            destination, force_scale=weight)
        self.assertLess(initial["normalized_free_residual_norm"], 0.25)
        self.assertLess(final["free_residual_norm"], 1.e-6 * weight)
        self.assertAlmostEqual(final["reaction_sum"][1], weight, places=6)
        scale = max(
            abs(value) for points in reference_stress.values()
            for point in points for value in point)
        self.assertLessEqual(
            _max_error(_stresses(destination), reference_stress),
            1.e-8 * scale)
        for node_id in [node.Id for node in interface_nodes]:
            node = destination.Nodes[node_id]
            if (node.IsFixed(KM.DISPLACEMENT_X) or
                    node.IsFixed(KM.DISPLACEMENT_Y) or
                    node.IsFixed(KM.DISPLACEMENT_Z)):
                continue
            residual = final["residual_by_node"][node_id]
            self.assertLess(
                math.sqrt(sum(value * value for value in residual)), 1.e-3)


if __name__ == "__main__":
    KratosUnittest.main()
