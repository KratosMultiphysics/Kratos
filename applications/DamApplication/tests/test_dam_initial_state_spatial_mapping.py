"""Spatial mapping of Dam initial state between nonmatching 2D meshes."""
import json
import math
import tempfile
from pathlib import Path

import KratosMultiphysics as KM
import KratosMultiphysics.DamApplication as Dam
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.ConvectionDiffusionApplication as CDA

Dam.__path__.insert(0, str(Path(__file__).resolve().parents[1] / "python_scripts"))

from KratosMultiphysics.DamApplication import dam_initial_state_file
from KratosMultiphysics.DamApplication import dam_initial_state_mapping
from KratosMultiphysics.DamApplication.export_dam_final_state_process import (
    ExportDamFinalStateProcess)
from KratosMultiphysics.DamApplication.import_dam_initial_state_process import (
    ImportDamInitialStateProcess)


def _nodes(nx, ny, lx=1.0, ly=2.0):
    nodes = {}
    for j in range(ny + 1):
        for i in range(nx + 1):
            nodes[j * (nx + 1) + i + 1] = [lx * i / nx, ly * j / ny, 0.0]
    return nodes


def _q4_mesh(nx, ny, lx=1.0, ly=2.0):
    nodes = _nodes(nx, ny, lx, ly)
    elements = []
    element_id = 1
    for j in range(ny):
        for i in range(nx):
            n0 = j * (nx + 1) + i + 1
            elements.append(
                (element_id, [n0, n0 + 1, n0 + nx + 2, n0 + nx + 1]))
            element_id += 1
    return {"nodes": nodes, "elements": elements, "geometry_type": 8,
            "element": "SmallDisplacementElement2D4N"}


def _t3_mesh(nx, ny, lx=1.0, ly=2.0):
    nodes = _nodes(nx, ny, lx, ly)
    elements = []
    element_id = 1
    for j in range(ny):
        for i in range(nx):
            n0 = j * (nx + 1) + i + 1
            n1 = n0 + 1
            n2 = n0 + nx + 2
            n3 = n0 + nx + 1
            elements.append((element_id, [n0, n1, n2]))
            element_id += 1
            elements.append((element_id, [n0, n2, n3]))
            element_id += 1
    return {"nodes": nodes, "elements": elements, "geometry_type": 16,
            "element": "SmallDisplacementElement2D3N"}


def _add_mid_nodes(nodes, a, b, mid_nodes, next_id):
    key = tuple(sorted((a, b)))
    if key not in mid_nodes:
        mid_nodes[key] = next_id
        nodes[next_id] = [
            0.5 * (nodes[a][0] + nodes[b][0]),
            0.5 * (nodes[a][1] + nodes[b][1]), 0.0]
        next_id += 1
    return mid_nodes[key], next_id


def _t6_mesh(nx, ny, lx=1.0, ly=2.0):
    nodes = _nodes(nx, ny, lx, ly)
    next_id = max(nodes) + 1
    mid_nodes = {}
    elements = []
    element_id = 1
    for j in range(ny):
        for i in range(nx):
            n0 = j * (nx + 1) + i + 1
            n1 = n0 + 1
            n2 = n0 + nx + 2
            n3 = n0 + nx + 1
            for a, b, c in ((n0, n1, n2), (n0, n2, n3)):
                m0, next_id = _add_mid_nodes(nodes, a, b, mid_nodes, next_id)
                m1, next_id = _add_mid_nodes(nodes, b, c, mid_nodes, next_id)
                m2, next_id = _add_mid_nodes(nodes, c, a, mid_nodes, next_id)
                elements.append((element_id, [a, b, c, m0, m1, m2]))
                element_id += 1
    return {"nodes": nodes, "elements": elements, "geometry_type": 17,
            "element": "SmallDisplacementElement2D6N"}


def _q8_mesh(nx, ny, lx=1.0, ly=2.0):
    nodes = _nodes(nx, ny, lx, ly)
    next_id = max(nodes) + 1
    mid_nodes = {}
    elements = []
    element_id = 1
    for j in range(ny):
        for i in range(nx):
            n0 = j * (nx + 1) + i + 1
            n1 = n0 + 1
            n2 = n0 + nx + 2
            n3 = n0 + nx + 1
            m0, next_id = _add_mid_nodes(nodes, n0, n1, mid_nodes, next_id)
            m1, next_id = _add_mid_nodes(nodes, n1, n2, mid_nodes, next_id)
            m2, next_id = _add_mid_nodes(nodes, n2, n3, mid_nodes, next_id)
            m3, next_id = _add_mid_nodes(nodes, n3, n0, mid_nodes, next_id)
            elements.append((element_id, [n0, n1, n2, n3, m0, m1, m2, m3]))
            element_id += 1
    return {"nodes": nodes, "elements": elements, "geometry_type": 9,
            "element": "SmallDisplacementElement2D8N"}


def _q9_mesh(nx, ny, lx=1.0, ly=2.0):
    mesh = _q8_mesh(nx, ny, lx, ly)
    nodes = mesh["nodes"]
    next_id = max(nodes) + 1
    elements = []
    for element_id, connectivity in mesh["elements"]:
        corners = connectivity[:4]
        center = next_id
        nodes[center] = [
            sum(nodes[node][0] for node in corners) / 4.0,
            sum(nodes[node][1] for node in corners) / 4.0, 0.0]
        next_id += 1
        elements.append((element_id, connectivity + [center]))
    return {"nodes": nodes, "elements": elements, "geometry_type": 10,
            "element": "SmallDisplacementElement2D9N"}


def _distort(mesh, amplitude, phase=0.0):
    nodes = {}
    for node_id, (x, y, z) in mesh["nodes"].items():
        if 0.0 < x < 1.0 and 0.0 < y < 2.0:
            x += amplitude * math.sin(phase + 7.0 * y + 3.0 * x)
            y += amplitude * math.sin(phase + 5.0 * x + 2.0 * y)
        nodes[node_id] = [x, y, z]
    return {**mesh, "nodes": nodes}


def _curved_q8_mesh():
    return {
        "nodes": {
            1: [0.0, 0.0, 0.0], 2: [1.0, 0.0, 0.0],
            3: [1.0, 1.0, 0.0], 4: [0.0, 1.0, 0.0],
            5: [0.5, -0.1, 0.0], 6: [1.1, 0.5, 0.0],
            7: [0.5, 1.1, 0.0], 8: [-0.1, 0.5, 0.0]},
        "elements": [(1, [1, 2, 3, 4, 5, 6, 7, 8])],
        "geometry_type": 9, "element": "SmallDisplacementElement2D8N"}


def _curved_q9_mesh():
    mesh = _curved_q8_mesh()
    mesh["nodes"][9] = [0.55, 0.45, 0.0]
    mesh["elements"] = [(1, list(range(1, 10)))]
    mesh["geometry_type"] = 10
    mesh["element"] = "SmallDisplacementElement2D9N"
    return mesh


def _curved_t6_mesh():
    return {
        "nodes": {
            1: [0.0, 0.0, 0.0], 2: [1.0, 0.0, 0.0], 3: [0.0, 1.0, 0.0],
            4: [0.5, -0.1, 0.0], 5: [0.6, 0.6, 0.0],
            6: [-0.1, 0.5, 0.0]},
        "elements": [(1, [1, 2, 3, 4, 5, 6])],
        "geometry_type": 17, "element": "SmallDisplacementElement2D6N"}


def _triangle_t3_mesh():
    return {
        "nodes": {1: [0.0, 0.0, 0.0], 2: [1.0, 0.0, 0.0],
                  3: [0.0, 1.0, 0.0]},
        "elements": [(1, [1, 2, 3])],
        "geometry_type": 16, "element": "SmallDisplacementElement2D3N"}


def _build_part(model, name, mesh, law, integration_order=None,
                thermal_expansion=None, body=None,
                temperature=20.0, reference_temperature=20.0,
                young_modulus=1.e4, poisson_ratio=0.2, density=2.0,
                fixed_bottom=False):
    part = model.CreateModelPart(name, 2)
    for variable in (KM.DISPLACEMENT, KM.REACTION, KM.VOLUME_ACCELERATION,
                     KM.CAUCHY_STRESS_VECTOR, KM.GREEN_LAGRANGE_STRAIN_VECTOR,
                     KM.TEMPERATURE, Dam.NODAL_REFERENCE_TEMPERATURE):
        part.AddNodalSolutionStepVariable(variable)
    part.ProcessInfo[KM.DOMAIN_SIZE] = 2
    part.ProcessInfo[KM.DELTA_TIME] = 1.0
    for node_id, (x, y, z) in mesh["nodes"].items():
        node = part.CreateNewNode(node_id, x, y, z)
        for dof, reaction in ((KM.DISPLACEMENT_X, KM.REACTION_X),
                              (KM.DISPLACEMENT_Y, KM.REACTION_Y),
                              (KM.DISPLACEMENT_Z, KM.REACTION_Z)):
            node.AddDof(dof, reaction)
        if fixed_bottom and abs(y) < 1.e-14:
            node.Fix(KM.DISPLACEMENT_X)
            node.Fix(KM.DISPLACEMENT_Y)
        if body is not None:
            node.SetSolutionStepValue(KM.VOLUME_ACCELERATION, body)
        node.SetSolutionStepValue(KM.TEMPERATURE, temperature)
        node.SetSolutionStepValue(
            Dam.NODAL_REFERENCE_TEMPERATURE, reference_temperature)
    properties = part.CreateNewProperties(1)
    properties[KM.YOUNG_MODULUS] = young_modulus
    properties[KM.POISSON_RATIO] = poisson_ratio
    properties[KM.DENSITY] = density
    properties[KM.THICKNESS] = 1.0
    if thermal_expansion is not None:
        properties[Dam.THERMAL_EXPANSION] = thermal_expansion
    if integration_order is not None:
        properties[KM.KratosGlobals.GetVariable("INTEGRATION_ORDER")] = \
            integration_order
    properties[KM.CONSTITUTIVE_LAW] = law
    for element_id, connectivity in mesh["elements"]:
        element = part.CreateNewElement(
            mesh["element"], element_id, connectivity, properties)
        element.Initialize(part.ProcessInfo)
    return part


def _integrations(part):
    result = {}
    for element in part.Elements:
        local_points = \
            Dam.DamInitialStateTransferUtility.GetIntegrationPointCoordinates(
                element)
        geometry = element.GetGeometry()
        coordinates = []
        for local in local_points:
            point = geometry.GlobalCoordinates(KM.Array3(local))
            coordinates.append([float(point[0]), float(point[1])])
        result[element.Id] = coordinates
    return result


def _assign_analytic_stresses(part, function, compensate=False):
    stresses = {}
    for element in part.Elements:
        values = []
        for x, y in _integrations(part)[element.Id]:
            values.append([float(value) for value in function(x, y)])
        stresses[element.Id] = values
    Dam.DamInitialStateTransferUtility.AssignInitialStresses(
        part, stresses, compensate)


def _make_strategy(part):
    builder = KM.ResidualBasedBlockBuilderAndSolver(
        KM.SkylineLUFactorizationSolver())
    scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
    strategy = KM.ResidualBasedLinearStrategy(
        part, scheme, builder, True, False, False, False)
    strategy.SetEchoLevel(0)
    strategy.Initialize()
    strategy.Check()
    return strategy


def _displacements(part):
    return [
        value for node in part.Nodes
        for value in list(node.GetSolutionStepValue(KM.DISPLACEMENT))[:2]]


def _shift_mesh(mesh, dx):
    return {
        "nodes": {
            node_id: [x + dx, y, z]
            for node_id, (x, y, z) in mesh["nodes"].items()},
        "elements": list(mesh["elements"]),
        "geometry_type": mesh["geometry_type"],
        "element": mesh["element"],
    }


def _scale_mesh(mesh, x0, x1, y0, y1):
    base = _q4_mesh(4, 4, lx=1.0, ly=1.0)
    nodes = {
        node_id: [x0 + (x1 - x0) * x, y0 + (y1 - y0) * y, z]
        for node_id, (x, y, z) in base["nodes"].items()}
    return {**base, "nodes": nodes}


def _renumber_mesh(mesh, node_offset, element_offset):
    return {
        "nodes": {
            node_id + node_offset: list(coordinates)
            for node_id, coordinates in mesh["nodes"].items()},
        "elements": [
            (element_id + element_offset,
             [node_id + node_offset for node_id in connectivity])
            for element_id, connectivity in mesh["elements"]],
        "geometry_type": mesh["geometry_type"],
        "element": mesh["element"],
    }


def _build_two_material_part(model, name, mesh, materials, body=None,
                             fixed_bottom=False, temperature=20.0,
                             reference_temperature=20.0):
    """materials: list of (subpart_name, predicate(x, y), properties_dict)."""
    part = model.CreateModelPart(name, 2)
    for variable in (KM.DISPLACEMENT, KM.REACTION, KM.VOLUME_ACCELERATION,
                     KM.CAUCHY_STRESS_VECTOR, KM.GREEN_LAGRANGE_STRAIN_VECTOR,
                     KM.TEMPERATURE, Dam.NODAL_REFERENCE_TEMPERATURE):
        part.AddNodalSolutionStepVariable(variable)
    part.ProcessInfo[KM.DOMAIN_SIZE] = 2
    part.ProcessInfo[KM.DELTA_TIME] = 1.0
    for node_id, (x, y, z) in mesh["nodes"].items():
        node = part.CreateNewNode(node_id, x, y, z)
        for dof, reaction in ((KM.DISPLACEMENT_X, KM.REACTION_X),
                              (KM.DISPLACEMENT_Y, KM.REACTION_Y),
                              (KM.DISPLACEMENT_Z, KM.REACTION_Z)):
            node.AddDof(dof, reaction)
        if fixed_bottom and abs(y) < 1.e-14:
            node.Fix(KM.DISPLACEMENT_X)
            node.Fix(KM.DISPLACEMENT_Y)
        if body is not None:
            node.SetSolutionStepValue(KM.VOLUME_ACCELERATION, body)
        node.SetSolutionStepValue(KM.TEMPERATURE, temperature)
        node.SetSolutionStepValue(
            Dam.NODAL_REFERENCE_TEMPERATURE, reference_temperature)
    subparts = {}
    for index, (subpart_name, predicate, values) in enumerate(materials, 1):
        selected = []
        for element_id, connectivity in mesh["elements"]:
            x = sum(mesh["nodes"][node][0] for node in connectivity) \
                / len(connectivity)
            y = sum(mesh["nodes"][node][1] for node in connectivity) \
                / len(connectivity)
            if predicate(x, y):
                selected.append((element_id, connectivity))
        if not selected:
            raise RuntimeError(
                "No destination element was assigned to '{}'."
                .format(subpart_name))
        subpart = part.CreateSubModelPart(subpart_name)
        node_ids = sorted({
            node_id for _, connectivity in selected
            for node_id in connectivity})
        subpart.AddNodes(node_ids)
        properties = subpart.CreateNewProperties(index)
        properties[KM.YOUNG_MODULUS] = values["young_modulus"]
        properties[KM.POISSON_RATIO] = values["poisson_ratio"]
        properties[KM.DENSITY] = values.get("density", 2.0)
        properties[KM.THICKNESS] = 1.0
        properties[KM.CONSTITUTIVE_LAW] = values["law"]
        for element_id, connectivity in selected:
            element = part.CreateNewElement(
                mesh["element"], element_id, connectivity, properties)
            subpart.AddElements([element_id])
            element.Initialize(part.ProcessInfo)
        subparts[subpart_name] = subpart
    return part, subparts


def _stresses(part):
    return {
        element.Id: [
            [float(value) for value in point]
            for point in element.CalculateOnIntegrationPoints(
                KM.CAUCHY_STRESS_VECTOR, part.ProcessInfo)]
        for element in part.Elements}


def _write_state(path, part, stresses=True, temperature=False, reference=False):
    data = dam_initial_state_file.BuildFileData(
        part, stresses, temperature, reference)
    dam_initial_state_file.WriteFile(path, data)
    return data


def _import(model, part_name, path, stresses=True, temperature=False,
            reference=False, policy="preserve"):
    settings = KM.Parameters(json.dumps({
        "model_part_name": part_name,
        "input_file_name": str(path),
        "import_stresses": stresses,
        "import_temperature": temperature,
        "import_reference_temperature": reference,
        "reference_temperature_policy": policy,
        "transfer_method": "spatial_mapping",
        "echo_level": 0}))
    process = ImportDamInitialStateProcess(model, settings)
    process.PrevalidateInitialState()
    process.ExecuteBeforeSolutionLoop()
    return process


def _max_error(actual, expected):
    maximum = 0.0
    for element_id, points in expected.items():
        for point, expected_point in zip(actual[element_id], points):
            for value, expected_value in zip(point, expected_point):
                maximum = max(maximum, abs(value - expected_value))
    return maximum


def _constant_function(value):
    return lambda x, y: value


def _linear_function(x, y):
    return [100.0 + 10.0 * x + 20.0 * y,
            -50.0 + 5.0 * x - 7.0 * y,
            30.0 + 7.0 * x - 3.0 * y]


def _biquadratic_function(x, y):
    return [100.0 + 10.0 * x + 20.0 * y + 2.0 * x * x + 3.0 * x * y - y * y,
            -50.0 + 5.0 * x - 7.0 * y + x * x - 2.0 * x * y + 4.0 * y * y,
            30.0 + 7.0 * x - 3.0 * y + 5.0 * x * x - x * y + 2.0 * y * y]


class TestDamInitialStateSpatialMapping(KratosUnittest.TestCase):

    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.directory.cleanup()

    def _state_file(self, name):
        return Path(self.directory.name) / name

    def _roundtrip(self, source_mesh, source_law, destination_mesh,
                   destination_law, function, source_order=None,
                   destination_order=None, source_thermal_expansion=None,
                   tolerance=1.e-9):
        source_model = KM.Model()
        source = _build_part(
            source_model, "source", source_mesh, source_law,
            integration_order=source_order,
            thermal_expansion=source_thermal_expansion)
        _assign_analytic_stresses(source, function)
        path = self._state_file("state.json")
        _write_state(path, source, stresses=True)

        destination_model = KM.Model()
        destination = _build_part(
            destination_model, "destination", destination_mesh,
            destination_law,
            integration_order=destination_order,
            thermal_expansion=source_thermal_expansion)
        _import(destination_model, "destination", path, stresses=True)
        actual = _stresses(destination)
        expected = {
            element.Id: [list(function(x, y))
                         for x, y in _integrations(destination)[element.Id]]
            for element in destination.Elements}
        error = _max_error(actual, expected)
        self.assertLessEqual(error, tolerance)
        for node in destination.Nodes:
            for value in node.GetSolutionStepValue(KM.DISPLACEMENT):
                self.assertEqual(value, 0.0)

    def test_constant_stress_q4_coarse_to_fine(self):
        self._roundtrip(
            _q4_mesh(1, 2), SMA.LinearElasticPlaneStrain2DLaw(),
            _q4_mesh(4, 8), SMA.LinearElasticPlaneStrain2DLaw(),
            _constant_function([100.0, -50.0, 30.0]))

    def test_constant_stress_q4_fine_to_coarse(self):
        self._roundtrip(
            _q4_mesh(4, 8), SMA.LinearElasticPlaneStrain2DLaw(),
            _q4_mesh(1, 2), SMA.LinearElasticPlaneStrain2DLaw(),
            _constant_function([100.0, -50.0, 30.0]))

    def test_different_destination_ids(self):
        self._roundtrip(
            _q4_mesh(1, 2), SMA.LinearElasticPlaneStrain2DLaw(),
            _renumber_mesh(_q4_mesh(4, 8), 1000, 2000),
            SMA.LinearElasticPlaneStrain2DLaw(),
            _linear_function)

    def test_linear_stress_q4_coarse_to_fine(self):
        self._roundtrip(
            _q4_mesh(1, 2), SMA.LinearElasticPlaneStrain2DLaw(),
            _q4_mesh(4, 8), SMA.LinearElasticPlaneStrain2DLaw(),
            _linear_function)

    def test_linear_stress_q4_fine_to_coarse(self):
        self._roundtrip(
            _q4_mesh(4, 8), SMA.LinearElasticPlaneStrain2DLaw(),
            _q4_mesh(1, 2), SMA.LinearElasticPlaneStrain2DLaw(),
            _linear_function)

    def test_constant_stress_t3_to_q4(self):
        self._roundtrip(
            _t3_mesh(1, 2), SMA.LinearElasticPlaneStrain2DLaw(),
            _q4_mesh(4, 8), SMA.LinearElasticPlaneStrain2DLaw(),
            _constant_function([100.0, -50.0, 30.0]))

    def test_linear_stress_q4_to_t3(self):
        # The destination T3 Gauss points are inside source Q4 elements; the
        # Q4 reconstruction has enough information for a linear field.
        self._roundtrip(
            _q4_mesh(4, 8), SMA.LinearElasticPlaneStrain2DLaw(),
            _t3_mesh(1, 2), SMA.LinearElasticPlaneStrain2DLaw(),
            _linear_function)

    def test_linear_stress_t6_to_q4(self):
        self._roundtrip(
            _t6_mesh(1, 2), SMA.LinearElasticPlaneStrain2DLaw(),
            _q4_mesh(4, 8), SMA.LinearElasticPlaneStrain2DLaw(),
            _linear_function)

    def test_biquadratic_stress_q9_to_q4(self):
        # Q9 has nine Gauss points and exercises the Dam bounding-box search
        # configuration, because Quadrilateral2D9 has no HasIntersection.
        self._roundtrip(
            _q9_mesh(1, 2), SMA.LinearElasticPlaneStrain2DLaw(),
            _q4_mesh(4, 8), SMA.LinearElasticPlaneStrain2DLaw(),
            _biquadratic_function, tolerance=1.e-8)

    def test_linear_stress_q8_to_q4(self):
        self._roundtrip(
            _q8_mesh(1, 2), SMA.LinearElasticPlaneStrain2DLaw(),
            _q4_mesh(4, 8), SMA.LinearElasticPlaneStrain2DLaw(),
            _linear_function, tolerance=1.e-8)

    def test_linear_stress_distorted_q4(self):
        self._roundtrip(
            _distort(_q4_mesh(2, 4), 0.03, 0.0),
            SMA.LinearElasticPlaneStrain2DLaw(),
            _distort(_q4_mesh(4, 8), 0.03, 1.0),
            SMA.LinearElasticPlaneStrain2DLaw(),
            _linear_function, tolerance=1.e-8)

    def test_linear_stress_distorted_t6(self):
        self._roundtrip(
            _distort(_t6_mesh(2, 4), 0.02, 0.0),
            SMA.LinearElasticPlaneStrain2DLaw(),
            _distort(_t6_mesh(4, 8), 0.02, 1.0),
            SMA.LinearElasticPlaneStrain2DLaw(),
            _linear_function, source_order=4, tolerance=1.e-8)

    def test_linear_stress_distorted_q9(self):
        self._roundtrip(
            _distort(_q9_mesh(2, 4), 0.02, 0.0),
            SMA.LinearElasticPlaneStrain2DLaw(),
            _distort(_q9_mesh(4, 8), 0.02, 1.0),
            SMA.LinearElasticPlaneStrain2DLaw(),
            _linear_function, tolerance=1.e-8)

    def test_linear_stress_curved_q8(self):
        self._roundtrip(
            _curved_q8_mesh(), SMA.LinearElasticPlaneStrain2DLaw(),
            _scale_mesh(None, 0.2, 0.8, 0.2, 0.8),
            SMA.LinearElasticPlaneStrain2DLaw(),
            _linear_function, tolerance=1.e-8)

    def test_linear_stress_curved_q9(self):
        self._roundtrip(
            _curved_q9_mesh(), SMA.LinearElasticPlaneStrain2DLaw(),
            _scale_mesh(None, 0.2, 0.8, 0.2, 0.8),
            SMA.LinearElasticPlaneStrain2DLaw(),
            _linear_function, tolerance=1.e-8)

    def test_linear_stress_curved_t6(self):
        self._roundtrip(
            _curved_t6_mesh(), SMA.LinearElasticPlaneStrain2DLaw(),
            _triangle_t3_mesh(), SMA.LinearElasticPlaneStrain2DLaw(),
            _linear_function, source_order=4, tolerance=1.e-8)

    def test_linear_stress_q4_to_t6(self):
        self._roundtrip(
            _q4_mesh(4, 8), SMA.LinearElasticPlaneStrain2DLaw(),
            _t6_mesh(1, 2), SMA.LinearElasticPlaneStrain2DLaw(),
            _linear_function, tolerance=1.e-8)

    def test_linear_stress_q4_to_q8(self):
        self._roundtrip(
            _q4_mesh(4, 8), SMA.LinearElasticPlaneStrain2DLaw(),
            _q8_mesh(1, 2), SMA.LinearElasticPlaneStrain2DLaw(),
            _linear_function, tolerance=1.e-8)

    def test_linear_stress_q4_to_q9(self):
        self._roundtrip(
            _q4_mesh(4, 8), SMA.LinearElasticPlaneStrain2DLaw(),
            _q9_mesh(1, 2), SMA.LinearElasticPlaneStrain2DLaw(),
            _linear_function, tolerance=1.e-8)

    def test_boundary_nodes_within_tolerance(self):
        source_model = KM.Model()
        source = _build_part(
            source_model, "source", _q4_mesh(1, 2),
            SMA.LinearElasticPlaneStrain2DLaw(), temperature=20.0)
        path = self._state_file("boundary.json")
        _write_state(path, source, stresses=False, temperature=True)
        destination_mesh = _q4_mesh(2, 4)
        shifted = {
            node_id: [x + 2.e-10, y, z]
            for node_id, (x, y, z) in destination_mesh["nodes"].items()}
        destination_mesh = {**destination_mesh, "nodes": shifted}
        destination_model = KM.Model()
        destination = _build_part(
            destination_model, "destination", destination_mesh,
            SMA.LinearElasticPlaneStrain2DLaw())
        _import(destination_model, "destination", path, stresses=False,
                temperature=True)
        for node in destination.Nodes:
            self.assertAlmostEqual(
                node.GetSolutionStepValue(KM.TEMPERATURE), 20.0, places=10)

    def test_clearly_outside_boundary_is_rejected(self):
        source_model = KM.Model()
        source = _build_part(
            source_model, "source", _q4_mesh(1, 2),
            SMA.LinearElasticPlaneStrain2DLaw(), temperature=20.0)
        path = self._state_file("boundary.json")
        _write_state(path, source, stresses=False, temperature=True)
        destination_mesh = _q4_mesh(2, 4)
        shifted = {
            node_id: [x + 1.e-6, y, z]
            for node_id, (x, y, z) in destination_mesh["nodes"].items()}
        destination_mesh = {**destination_mesh, "nodes": shifted}
        destination_model = KM.Model()
        _build_part(destination_model, "destination", destination_mesh,
                    SMA.LinearElasticPlaneStrain2DLaw())
        with self.assertRaisesRegex(RuntimeError, "outside the source domain"):
            _import(destination_model, "destination", path, stresses=False,
                    temperature=True)

    def test_reduced_q4_integration_constant(self):
        self._roundtrip(
            _q4_mesh(1, 2), SMA.LinearElasticPlaneStrain2DLaw(),
            _q4_mesh(4, 8), SMA.LinearElasticPlaneStrain2DLaw(),
            _constant_function([100.0, -50.0, 30.0]), source_order=1)

    def test_q4_overintegration_linear(self):
        self._roundtrip(
            _q4_mesh(1, 2), SMA.LinearElasticPlaneStrain2DLaw(),
            _q4_mesh(4, 8), SMA.LinearElasticPlaneStrain2DLaw(),
            _linear_function, source_order=3, tolerance=1.e-8)

    def test_uniform_temperature(self):
        source_model = KM.Model()
        source = _build_part(
            source_model, "source", _q4_mesh(1, 2),
            SMA.LinearElasticPlaneStrain2DLaw(),
            temperature=25.0, reference_temperature=15.0)
        path = self._state_file("temperature.json")
        _write_state(path, source, stresses=False, temperature=True, reference=True)

        destination_model = KM.Model()
        destination = _build_part(
            destination_model, "destination", _q4_mesh(4, 8),
            SMA.LinearElasticPlaneStrain2DLaw(),
            temperature=0.0, reference_temperature=0.0)
        _import(destination_model, "destination", path, stresses=False,
                temperature=True, reference=True)
        for node in destination.Nodes:
            self.assertAlmostEqual(
                node.GetSolutionStepValue(KM.TEMPERATURE), 25.0, places=10)
            self.assertAlmostEqual(
                node.GetSolutionStepValue(
                    Dam.NODAL_REFERENCE_TEMPERATURE), 15.0, places=10)

    def test_initial_temperature_policy(self):
        source_model = KM.Model()
        source = _build_part(
            source_model, "source", _q4_mesh(1, 2),
            SMA.LinearElasticPlaneStrain2DLaw(),
            temperature=25.0, reference_temperature=15.0)
        path = self._state_file("temperature_policy.json")
        _write_state(path, source, stresses=False, temperature=True)

        destination_model = KM.Model()
        destination = _build_part(
            destination_model, "destination", _q4_mesh(4, 8),
            SMA.LinearElasticPlaneStrain2DLaw(),
            temperature=0.0, reference_temperature=0.0)
        _import(destination_model, "destination", path, stresses=False,
                temperature=True, policy="initial_temperature")
        for node in destination.Nodes:
            self.assertAlmostEqual(
                node.GetSolutionStepValue(KM.TEMPERATURE), 25.0, places=10)
            self.assertAlmostEqual(
                node.GetSolutionStepValue(
                    Dam.NODAL_REFERENCE_TEMPERATURE), 25.0, places=10)

    def test_linear_independent_temperatures(self):
        source_model = KM.Model()
        source = _build_part(
            source_model, "source", _q4_mesh(2, 4),
            SMA.LinearElasticPlaneStrain2DLaw())
        for node in source.Nodes:
            node.SetSolutionStepValue(
                KM.TEMPERATURE, 10.0 + 3.0 * node.X - 2.0 * node.Y)
            node.SetSolutionStepValue(
                Dam.NODAL_REFERENCE_TEMPERATURE,
                5.0 + node.X + node.Y)
        path = self._state_file("temperatures.json")
        _write_state(path, source, stresses=False, temperature=True, reference=True)

        destination_model = KM.Model()
        destination = _build_part(
            destination_model, "destination", _q4_mesh(5, 10),
            SMA.LinearElasticPlaneStrain2DLaw())
        _import(destination_model, "destination", path, stresses=False,
                temperature=True, reference=True)
        differing = 0
        for node in destination.Nodes:
            temperature = 10.0 + 3.0 * node.X - 2.0 * node.Y
            reference = 5.0 + node.X + node.Y
            self.assertAlmostEqual(
                node.GetSolutionStepValue(KM.TEMPERATURE), temperature, places=10)
            self.assertAlmostEqual(
                node.GetSolutionStepValue(
                    Dam.NODAL_REFERENCE_TEMPERATURE), reference, places=10)
            if abs(temperature - reference) > 1.e-12:
                differing += 1
        self.assertGreater(differing, 0)

    def test_outside_domain_is_rejected(self):
        source_model = KM.Model()
        source = _build_part(
            source_model, "source", _q4_mesh(1, 2),
            SMA.LinearElasticPlaneStrain2DLaw())
        _assign_analytic_stresses(
            source, _constant_function([100.0, -50.0, 30.0]))
        path = self._state_file("state.json")
        _write_state(path, source, stresses=True)

        destination_model = KM.Model()
        destination = _build_part(
            destination_model, "destination", _q4_mesh(2, 4, lx=1.0, ly=3.0),
            SMA.LinearElasticPlaneStrain2DLaw())
        with self.assertRaisesRegex(RuntimeError, "outside the source domain"):
            _import(destination_model, "destination", path, stresses=True)
        for node in destination.Nodes:
            self.assertEqual(
                node.GetSolutionStepValue(KM.TEMPERATURE), 20.0)
            for value in node.GetSolutionStepValue(KM.DISPLACEMENT):
                self.assertEqual(value, 0.0)

    def test_thermal_expansion_is_irrelevant_for_mechanical_law(self):
        source_model = KM.Model()
        source = _build_part(
            source_model, "source", _q4_mesh(1, 2),
            SMA.LinearElasticPlaneStrain2DLaw(),
            thermal_expansion=1.e-5)
        _assign_analytic_stresses(
            source, _constant_function([100.0, -50.0, 30.0]))
        path = self._state_file("state.json")
        _write_state(path, source, stresses=True)

        destination_model = KM.Model()
        destination = _build_part(
            destination_model, "destination", _q4_mesh(2, 4),
            SMA.LinearElasticPlaneStrain2DLaw())
        _import(destination_model, "destination", path, stresses=True)
        for element in destination.Elements:
            for point in element.CalculateOnIntegrationPoints(
                    KM.CAUCHY_STRESS_VECTOR, destination.ProcessInfo):
                self.assertAlmostEqual(point[0], 100.0, places=8)
                self.assertAlmostEqual(point[1], -50.0, places=8)
                self.assertAlmostEqual(point[2], 30.0, places=8)

    def test_material_incompatibility_is_rejected(self):
        source_model = KM.Model()
        source = _build_part(
            source_model, "source", _q4_mesh(1, 2),
            SMA.LinearElasticPlaneStrain2DLaw(),
            thermal_expansion=1.e-5)
        _assign_analytic_stresses(
            source, _constant_function([100.0, -50.0, 30.0]))
        path = self._state_file("state.json")
        _write_state(path, source, stresses=True)

        destination_model = KM.Model()
        destination = _build_part(
            destination_model, "destination", _q4_mesh(2, 4),
            SMA.LinearElasticPlaneStress2DLaw(),
            thermal_expansion=1.e-5)
        with self.assertRaisesRegex(RuntimeError, "material"):
            _import(destination_model, "destination", path, stresses=True)

    def test_destination_reduced_integration_linear(self):
        self._roundtrip(
            _q4_mesh(1, 2), SMA.LinearElasticPlaneStrain2DLaw(),
            _q4_mesh(4, 8), SMA.LinearElasticPlaneStrain2DLaw(),
            _linear_function, destination_order=1, tolerance=1.e-8)

    def test_thermal_stress_with_independent_temperatures(self):
        law = Dam.ThermalLinearElastic2DPlaneStrain()
        source_model = KM.Model()
        source = _build_part(
            source_model, "source", _q4_mesh(1, 2), law,
            thermal_expansion=1.e-5, temperature=30.0,
            reference_temperature=20.0)
        _assign_analytic_stresses(
            source, _constant_function([100.0, -50.0, 30.0]),
            compensate=True)
        path = self._state_file("thermal.json")
        _write_state(path, source, stresses=True, temperature=True,
                     reference=True)

        destination_model = KM.Model()
        destination = _build_part(
            destination_model, "destination", _q4_mesh(4, 8), law,
            thermal_expansion=1.e-5, temperature=5.0,
            reference_temperature=5.0)
        _import(destination_model, "destination", path, stresses=True,
                temperature=True, reference=True)
        for element in destination.Elements:
            for point in element.CalculateOnIntegrationPoints(
                    KM.CAUCHY_STRESS_VECTOR, destination.ProcessInfo):
                self.assertAlmostEqual(point[0], 100.0, places=8)
                self.assertAlmostEqual(point[1], -50.0, places=8)
                self.assertAlmostEqual(point[2], 30.0, places=8)
        for node in destination.Nodes:
            self.assertAlmostEqual(
                node.GetSolutionStepValue(KM.TEMPERATURE), 30.0, places=10)
            self.assertAlmostEqual(
                node.GetSolutionStepValue(
                    Dam.NODAL_REFERENCE_TEMPERATURE), 20.0, places=10)
            for value in node.GetSolutionStepValue(KM.DISPLACEMENT):
                self.assertEqual(value, 0.0)

    def _gravity_case(self, source_n, dest_n=8):
        body = [0.0, -10.0, 0.0]
        source_model = KM.Model()
        source = _build_part(
            source_model, "source", _q4_mesh(source_n, 2 * source_n),
            SMA.LinearElasticPlaneStrain2DLaw(), body=body,
            fixed_bottom=True, temperature=30.0,
            reference_temperature=20.0)
        _make_strategy(source).Solve()
        path = self._state_file("gravity_{}.json".format(source_n))
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
            reference_model, "reference", _q4_mesh(dest_n, 2 * dest_n),
            SMA.LinearElasticPlaneStrain2DLaw(), body=body,
            fixed_bottom=True)
        _make_strategy(reference).Solve()
        reference_stress = _stresses(reference)
        reference_displacement = _displacements(reference)
        reference_norm = math.sqrt(
            sum(value * value for value in reference_displacement))

        destination_model = KM.Model()
        destination = _build_part(
            destination_model, "destination", _q4_mesh(dest_n, 2 * dest_n),
            SMA.LinearElasticPlaneStrain2DLaw(), body=body,
            fixed_bottom=True)
        destination_strategy = _make_strategy(destination)
        _import(destination_model, "destination", path, stresses=True,
                temperature=True, reference=True)

        for node in destination.Nodes:
            for value in node.GetSolutionStepValue(KM.DISPLACEMENT):
                self.assertEqual(value, 0.0)
            self.assertAlmostEqual(
                node.GetSolutionStepValue(KM.TEMPERATURE), 30.0, places=10)
            self.assertAlmostEqual(
                node.GetSolutionStepValue(
                    Dam.NODAL_REFERENCE_TEMPERATURE), 20.0, places=10)

        diagnostic = dam_initial_state_mapping.EvaluateMechanicalResidual(
            destination, force_scale=40.0)
        free_residual = diagnostic["free_residual_norm"]

        destination_strategy.Solve()
        correction = _displacements(destination)
        correction_norm = math.sqrt(
            sum(value * value for value in correction))
        stress_error = _max_error(_stresses(destination), reference_stress)

        diagnostic = dam_initial_state_mapping.EvaluateMechanicalResidual(
            destination, force_scale=40.0)
        return {
            "free_residual": free_residual,
            "correction": correction_norm,
            "reference_displacement": reference_norm,
            "stress_error": stress_error,
            "reaction_y": diagnostic["reaction_sum"][1],
        }

    def test_gravity_equilibrium_and_first_step(self):
        metrics = self._gravity_case(4, 8)
        self.assertLess(metrics["free_residual"], 0.2 * 40.0)
        self.assertLess(metrics["correction"], 0.05 * metrics["reference_displacement"])
        self.assertLess(metrics["stress_error"], 1.e-6)
        self.assertAlmostEqual(metrics["reaction_y"], 40.0, places=6)

    def test_gravity_converges_with_source_refinement(self):
        coarse = self._gravity_case(2, 8)
        fine = self._gravity_case(4, 8)
        self.assertLess(fine["free_residual"], coarse["free_residual"])
        self.assertLess(fine["correction"], coarse["correction"])
        self.assertAlmostEqual(coarse["reaction_y"], 40.0, places=6)
        self.assertAlmostEqual(fine["reaction_y"], 40.0, places=6)

    def test_equilibrium_diagnostic_includes_conditions(self):
        model = KM.Model()
        part = _build_part(
            model, "part", _q4_mesh(1, 2),
            SMA.LinearElasticPlaneStrain2DLaw())
        loaded_node = max(part.Nodes, key=lambda node: node.Y)
        condition = part.CreateNewCondition(
            "PointLoadCondition2D1N", 1, [loaded_node.Id],
            part.GetProperties(1))
        condition.SetValue(SMA.POINT_LOAD, [3.0, 7.0, 0.0])
        diagnostic = dam_initial_state_mapping.EvaluateMechanicalResidual(part)
        self.assertAlmostEqual(
            max(abs(value) for value in diagnostic["free_residual"]),
            7.0, places=10)

    def test_equilibrium_diagnostic_filters_thermal_conditions(self):
        model = KM.Model()
        part = _build_part(
            model, "part", _q4_mesh(1, 2),
            SMA.LinearElasticPlaneStrain2DLaw())
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
        top_nodes = [node.Id for node in part.Nodes if abs(node.Y - 2.0) < 1.e-14]
        thermal = part.CreateNewCondition(
            "ThermalFace2D2N", 1, sorted(top_nodes), properties)
        thermal_rhs = KM.Vector()
        thermal.CalculateRightHandSide(thermal_rhs, part.ProcessInfo)
        self.assertGreater(max(abs(value) for value in thermal_rhs), 1.0)

        line_load = part.CreateNewCondition(
            "LineLoadCondition2D2N", 2, sorted(top_nodes), properties)
        line_load.SetValue(SMA.LINE_LOAD, [0.0, -5.0, 0.0])
        diagnostic = dam_initial_state_mapping.EvaluateMechanicalResidual(part)
        self.assertAlmostEqual(
            diagnostic["max_free_residual"], 2.5, places=10)
        for node_id in top_nodes:
            self.assertAlmostEqual(
                diagnostic["residual_by_node"][node_id][1], -2.5, places=10)

    def test_mechanical_dof_row_filter_handles_mixed_order(self):
        model = KM.Model()
        part = _build_part(
            model, "part", _q4_mesh(1, 1),
            SMA.LinearElasticPlaneStrain2DLaw())
        node = next(iter(part.Nodes))
        node.AddDof(KM.TEMPERATURE)

        class _Entity:
            def __init__(self, dofs):
                self.dofs = dofs

            def GetDofList(self, process_info):
                return self.dofs

        entity = _Entity([
            node.GetDof(KM.DISPLACEMENT_X),
            node.GetDof(KM.TEMPERATURE),
            node.GetDof(KM.DISPLACEMENT_Y)])
        rows, number_of_dofs = \
            dam_initial_state_mapping._MechanicalDofRows(
                entity, part.ProcessInfo)
        self.assertEqual(number_of_dofs, 3)
        self.assertEqual(rows, [(0, node.Id, 0), (2, node.Id, 1)])

    def test_equilibrium_diagnostic_is_read_only(self):
        model = KM.Model()
        part = _build_part(
            model, "part", _q4_mesh(2, 4),
            SMA.LinearElasticPlaneStrain2DLaw())
        _assign_analytic_stresses(
            part, _constant_function([100.0, -50.0, 30.0]))
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

    def test_equilibrium_diagnostic_traverses_subparts(self):
        model = KM.Model()
        root = model.CreateModelPart("root", 2)
        for variable in (KM.DISPLACEMENT, KM.REACTION, KM.VOLUME_ACCELERATION):
            root.AddNodalSolutionStepVariable(variable)
        root.ProcessInfo[KM.DOMAIN_SIZE] = 2
        mesh = _q4_mesh(1, 2)
        for node_id, (x, y, z) in mesh["nodes"].items():
            node = root.CreateNewNode(node_id, x, y, z)
            for dof, reaction in ((KM.DISPLACEMENT_X, KM.REACTION_X),
                                  (KM.DISPLACEMENT_Y, KM.REACTION_Y),
                                  (KM.DISPLACEMENT_Z, KM.REACTION_Z)):
                node.AddDof(dof, reaction)
            if y == 0.0:
                node.Fix(KM.DISPLACEMENT_X)
                node.Fix(KM.DISPLACEMENT_Y)
            node.SetSolutionStepValue(KM.VOLUME_ACCELERATION, [0.0, -10.0, 0.0])
        subpart = root.CreateSubModelPart("body")
        subpart.AddNodes(sorted(mesh["nodes"]))
        properties = subpart.CreateNewProperties(1)
        properties[KM.YOUNG_MODULUS] = 1.e4
        properties[KM.POISSON_RATIO] = 0.2
        properties[KM.DENSITY] = 2.0
        properties[KM.THICKNESS] = 1.0
        properties[KM.CONSTITUTIVE_LAW] = SMA.LinearElasticPlaneStrain2DLaw()
        for element_id, connectivity in mesh["elements"]:
            element = subpart.CreateNewElement(
                "SmallDisplacementElement2D4N", element_id, connectivity,
                properties)
            element.Initialize(root.ProcessInfo)
        diagnostic = dam_initial_state_mapping.EvaluateMechanicalResidual(root)
        self.assertGreater(diagnostic["free_residual_norm"], 1.0)

    def test_multipart_material_interface(self):
        concrete_mesh = _q4_mesh(2, 4)
        soil_mesh = _shift_mesh(_q4_mesh(2, 4), 1.0)
        concrete_value = [100.0, -50.0, 30.0]
        soil_value = [-20.0, 40.0, 10.0]

        concrete_model = KM.Model()
        concrete = _build_part(
            concrete_model, "Concrete", concrete_mesh,
            SMA.LinearElasticPlaneStrain2DLaw(),
            young_modulus=1.e4, poisson_ratio=0.2)
        _assign_analytic_stresses(
            concrete, _constant_function(concrete_value))
        concrete_path = self._state_file("concrete.json")
        _write_state(concrete_path, concrete, stresses=True)

        soil_model = KM.Model()
        soil = _build_part(
            soil_model, "Soil", soil_mesh,
            SMA.LinearElasticPlaneStrain2DLaw(),
            young_modulus=5.e3, poisson_ratio=0.25)
        _assign_analytic_stresses(soil, _constant_function(soil_value))
        soil_path = self._state_file("soil.json")
        _write_state(soil_path, soil, stresses=True)

        destination_mesh = _q4_mesh(16, 16, lx=2.0, ly=2.0)
        destination_model = KM.Model()
        _, destination_parts = _build_two_material_part(
            destination_model, "destination", destination_mesh,
            [("Parts_Concrete", lambda x, y: x < 1.0,
              {"law": SMA.LinearElasticPlaneStrain2DLaw(),
               "young_modulus": 1.e4, "poisson_ratio": 0.2}),
             ("Parts_Soil", lambda x, y: x >= 1.0,
              {"law": SMA.LinearElasticPlaneStrain2DLaw(),
               "young_modulus": 5.e3, "poisson_ratio": 0.25})])

        concrete_process = ImportDamInitialStateProcess(
            destination_model, KM.Parameters(json.dumps({
                "model_part_name": "destination.Parts_Concrete",
                "input_file_name": str(concrete_path),
                "import_stresses": True,
                "transfer_method": "spatial_mapping"})))
        soil_process = ImportDamInitialStateProcess(
            destination_model, KM.Parameters(json.dumps({
                "model_part_name": "destination.Parts_Soil",
                "input_file_name": str(soil_path),
                "import_stresses": True,
                "transfer_method": "spatial_mapping"})))
        concrete_process.PrevalidateInitialState()
        soil_process.PrevalidateInitialState()
        dam_initial_state_file.ValidateCombinedImports(
            [concrete_process, soil_process])
        concrete_process.ExecuteBeforeSolutionLoop()
        soil_process.ExecuteBeforeSolutionLoop()

        process_info = destination_model["destination"].ProcessInfo
        for element in destination_parts["Parts_Concrete"].Elements:
            for point in element.CalculateOnIntegrationPoints(
                    KM.CAUCHY_STRESS_VECTOR, process_info):
                for actual, expected in zip(point, concrete_value):
                    self.assertAlmostEqual(actual, expected, places=8)
        for element in destination_parts["Parts_Soil"].Elements:
            for point in element.CalculateOnIntegrationPoints(
                    KM.CAUCHY_STRESS_VECTOR, process_info):
                for actual, expected in zip(point, soil_value):
                    self.assertAlmostEqual(actual, expected, places=8)

    def test_two_material_gravity_spatial_mapping(self):
        body = [0.0, -10.0, 0.0]
        source_mesh = _q4_mesh(8, 16, lx=2.0, ly=2.0)
        destination_mesh = _q4_mesh(16, 32, lx=2.0, ly=2.0)
        materials = [
            ("Parts_Concrete", lambda x, y: x < 1.0,
             {"law": SMA.LinearElasticPlaneStrain2DLaw(),
              "young_modulus": 1.e4, "poisson_ratio": 0.2}),
            ("Parts_Soil", lambda x, y: x >= 1.0,
             {"law": SMA.LinearElasticPlaneStrain2DLaw(),
              "young_modulus": 5.e3, "poisson_ratio": 0.25}),
        ]
        source_model = KM.Model()
        source_root, source_parts = _build_two_material_part(
            source_model, "source", source_mesh, materials, body=body,
            fixed_bottom=True)
        _make_strategy(source_root).Solve()
        paths = {}
        for name, _, _ in materials:
            paths[name] = self._state_file(name + "_gravity.json.gz")
            ExportDamFinalStateProcess(
                source_model, KM.Parameters(json.dumps({
                    "model_part_name": "source." + name,
                    "export_stresses": True,
                    "output_file_name": str(paths[name]),
                    "echo_level": 0}))).ExecuteFinalize()

        reference_model = KM.Model()
        reference_root, _ = _build_two_material_part(
            reference_model, "reference", destination_mesh, materials,
            body=body, fixed_bottom=True)
        _make_strategy(reference_root).Solve()
        reference_stress = _stresses(reference_root)
        reference_u = _displacements(reference_root)
        reference_norm = math.sqrt(
            sum(value * value for value in reference_u))

        destination_model = KM.Model()
        destination_root, destination_parts = _build_two_material_part(
            destination_model, "destination", destination_mesh, materials,
            body=body, fixed_bottom=True)
        destination_strategy = _make_strategy(destination_root)
        processes = []
        for name, _, _ in materials:
            processes.append(ImportDamInitialStateProcess(
                destination_model, KM.Parameters(json.dumps({
                    "model_part_name": "destination." + name,
                    "input_file_name": str(paths[name]),
                    "import_stresses": True,
                    "transfer_method": "spatial_mapping"}))))
        for process in processes:
            process.PrevalidateInitialState()
        dam_initial_state_file.ValidateCombinedImports(processes)
        for process in processes:
            process.ExecuteBeforeSolutionLoop()
        for node in destination_root.Nodes:
            for value in node.GetSolutionStepValue(KM.DISPLACEMENT):
                self.assertEqual(value, 0.0)

        initial = dam_initial_state_mapping.EvaluateMechanicalResidual(
            destination_root, force_scale=80.0)
        before = _displacements(destination_root)
        destination_strategy.Solve()
        after = _displacements(destination_root)
        correction = [
            after[i] - before[i] for i in range(len(before))]
        correction_norm = math.sqrt(
            sum(value * value for value in correction))
        self.assertLess(correction_norm, 0.05 * reference_norm)
        self.assertLess(initial["free_residual_norm"], 0.2 * 80.0)
        self.assertLess(_max_error(_stresses(destination_root),
                                   reference_stress), 1.e-6)
        final = dam_initial_state_mapping.EvaluateMechanicalResidual(
            destination_root, force_scale=80.0)
        self.assertAlmostEqual(final["reaction_sum"][1], 80.0, places=6)
        interface_nodes = [
            nid for nid, (x, y, z) in destination_mesh["nodes"].items()
            if abs(x - 1.0) < 1.e-14 and y > 1.e-14]
        interface = max(
            math.hypot(*final["residual_by_node"][nid])
            for nid in interface_nodes)
        self.assertLess(interface, 1.e-8 * 80.0)

    def test_shared_node_temperature_consistency(self):
        self._shared_node_case(20.0, 20.0, expect_conflict=False)
        self._shared_node_case(20.0, 25.0, expect_conflict=True)

    def _shared_node_case(self, concrete_temperature, soil_temperature,
                          expect_conflict):
        concrete_model = KM.Model()
        concrete = _build_part(
            concrete_model, "Concrete", _q4_mesh(2, 4),
            SMA.LinearElasticPlaneStrain2DLaw(),
            temperature=concrete_temperature)
        concrete_path = self._state_file("concrete_temperature.json")
        _write_state(concrete_path, concrete, stresses=False, temperature=True)

        soil_model = KM.Model()
        soil = _build_part(
            soil_model, "Soil", _shift_mesh(_q4_mesh(2, 4), 1.0),
            SMA.LinearElasticPlaneStrain2DLaw(),
            young_modulus=5.e3, poisson_ratio=0.25,
            temperature=soil_temperature)
        soil_path = self._state_file("soil_temperature.json")
        _write_state(soil_path, soil, stresses=False, temperature=True)

        destination_model = KM.Model()
        destination_mesh = _q4_mesh(16, 16, lx=2.0, ly=2.0)
        _, destination_parts = _build_two_material_part(
            destination_model, "destination", destination_mesh,
            [("Parts_Concrete", lambda x, y: x < 1.0,
              {"law": SMA.LinearElasticPlaneStrain2DLaw(),
               "young_modulus": 1.e4, "poisson_ratio": 0.2}),
             ("Parts_Soil", lambda x, y: x >= 1.0,
              {"law": SMA.LinearElasticPlaneStrain2DLaw(),
               "young_modulus": 5.e3, "poisson_ratio": 0.25})])
        concrete_process = ImportDamInitialStateProcess(
            destination_model, KM.Parameters(json.dumps({
                "model_part_name": "destination.Parts_Concrete",
                "input_file_name": str(concrete_path),
                "import_stresses": False,
                "import_temperature": True,
                "transfer_method": "spatial_mapping"})))
        soil_process = ImportDamInitialStateProcess(
            destination_model, KM.Parameters(json.dumps({
                "model_part_name": "destination.Parts_Soil",
                "input_file_name": str(soil_path),
                "import_stresses": False,
                "import_temperature": True,
                "transfer_method": "spatial_mapping"})))
        concrete_process.PrevalidateInitialState()
        soil_process.PrevalidateInitialState()
        if expect_conflict:
            with self.assertRaisesRegex(RuntimeError, "shared node"):
                dam_initial_state_file.ValidateCombinedImports(
                    [concrete_process, soil_process])
        else:
            dam_initial_state_file.ValidateCombinedImports(
                [concrete_process, soil_process])
            concrete_process.ExecuteBeforeSolutionLoop()
            soil_process.ExecuteBeforeSolutionLoop()
            for node in destination_parts["Parts_Concrete"].Nodes:
                self.assertAlmostEqual(
                    node.GetSolutionStepValue(KM.TEMPERATURE), 20.0, places=10)

    def test_shared_node_linear_temperature_nonmatching_interfaces(self):
        concrete_model = KM.Model()
        concrete = _build_part(
            concrete_model, "Concrete", _q4_mesh(2, 4),
            SMA.LinearElasticPlaneStrain2DLaw())
        for node in concrete.Nodes:
            node.SetSolutionStepValue(KM.TEMPERATURE, 20.0 + 0.1 * node.Y)
        concrete_path = self._state_file("concrete_linear_temperature.json")
        _write_state(concrete_path, concrete, stresses=False, temperature=True)

        soil_model = KM.Model()
        soil = _build_part(
            soil_model, "Soil", _shift_mesh(_q4_mesh(2, 4), 1.0),
            SMA.LinearElasticPlaneStrain2DLaw(),
            young_modulus=5.e3, poisson_ratio=0.25)
        for node in soil.Nodes:
            node.SetSolutionStepValue(KM.TEMPERATURE, 20.0 + 0.1 * node.Y)
        soil_path = self._state_file("soil_linear_temperature.json")
        _write_state(soil_path, soil, stresses=False, temperature=True)

        destination_model = KM.Model()
        destination_mesh = _q4_mesh(16, 16, lx=2.0, ly=2.0)
        _build_two_material_part(
            destination_model, "destination", destination_mesh,
            [("Parts_Concrete", lambda x, y: x < 1.0,
              {"law": SMA.LinearElasticPlaneStrain2DLaw(),
               "young_modulus": 1.e4, "poisson_ratio": 0.2}),
             ("Parts_Soil", lambda x, y: x >= 1.0,
              {"law": SMA.LinearElasticPlaneStrain2DLaw(),
               "young_modulus": 5.e3, "poisson_ratio": 0.25})])
        concrete_process = ImportDamInitialStateProcess(
            destination_model, KM.Parameters(json.dumps({
                "model_part_name": "destination.Parts_Concrete",
                "input_file_name": str(concrete_path),
                "import_stresses": False,
                "import_temperature": True,
                "transfer_method": "spatial_mapping"})))
        soil_process = ImportDamInitialStateProcess(
            destination_model, KM.Parameters(json.dumps({
                "model_part_name": "destination.Parts_Soil",
                "input_file_name": str(soil_path),
                "import_stresses": False,
                "import_temperature": True,
                "transfer_method": "spatial_mapping"})))
        concrete_process.PrevalidateInitialState()
        soil_process.PrevalidateInitialState()
        dam_initial_state_file.ValidateCombinedImports(
            [concrete_process, soil_process])
        concrete_process.ExecuteBeforeSolutionLoop()
        soil_process.ExecuteBeforeSolutionLoop()
        for node in destination_model["destination"].Nodes:
            if abs(node.X - 1.0) < 1.e-14:
                self.assertAlmostEqual(
                    node.GetSolutionStepValue(KM.TEMPERATURE),
                    20.0 + 0.1 * node.Y, places=9)

    def test_nonzero_destination_displacement_is_rejected(self):
        source_model = KM.Model()
        source = _build_part(
            source_model, "source", _q4_mesh(1, 2),
            SMA.LinearElasticPlaneStrain2DLaw())
        _assign_analytic_stresses(
            source, _constant_function([100.0, -50.0, 30.0]))
        path = self._state_file("state.json")
        _write_state(path, source, stresses=True)

        destination_model = KM.Model()
        destination = _build_part(
            destination_model, "destination", _q4_mesh(2, 4),
            SMA.LinearElasticPlaneStrain2DLaw())
        next(iter(destination.Nodes)).SetSolutionStepValue(
            KM.DISPLACEMENT, [1.e-3, 0.0, 0.0])
        with self.assertRaisesRegex(RuntimeError, "nonzero displacement"):
            _import(destination_model, "destination", path, stresses=True)

    def test_duplicate_initial_state_is_rejected(self):
        source_model = KM.Model()
        source = _build_part(
            source_model, "source", _q4_mesh(1, 2),
            SMA.LinearElasticPlaneStrain2DLaw())
        _assign_analytic_stresses(
            source, _constant_function([100.0, -50.0, 30.0]))
        path = self._state_file("state.json")
        _write_state(path, source, stresses=True)

        destination_model = KM.Model()
        _build_part(destination_model, "destination", _q4_mesh(2, 4),
                    SMA.LinearElasticPlaneStrain2DLaw())
        _import(destination_model, "destination", path, stresses=True)
        second = ImportDamInitialStateProcess(
            destination_model, KM.Parameters(json.dumps({
                "model_part_name": "destination",
                "input_file_name": str(path),
                "import_stresses": True,
                "transfer_method": "spatial_mapping"})))
        with self.assertRaisesRegex(RuntimeError, "already has an initial state"):
            second.PrevalidateInitialState()

    def test_spatial_mapping_rejects_3d_source(self):
        source_model = KM.Model()
        source = _build_part(
            source_model, "source", _q4_mesh(1, 2),
            SMA.LinearElasticPlaneStrain2DLaw())
        _assign_analytic_stresses(
            source, _constant_function([100.0, -50.0, 30.0]))
        path = self._state_file("state.json")
        data = _write_state(path, source, stresses=True)
        data["dimension"] = 3
        dam_initial_state_file.WriteFile(path, data)

        destination_model = KM.Model()
        _build_part(destination_model, "destination", _q4_mesh(2, 4),
                    SMA.LinearElasticPlaneStrain2DLaw())
        with self.assertRaisesRegex(RuntimeError, "2D source files only"):
            _import(destination_model, "destination", path, stresses=True)

    def test_incomplete_file_is_rejected(self):
        source_model = KM.Model()
        source = _build_part(
            source_model, "source", _q4_mesh(1, 2),
            SMA.LinearElasticPlaneStrain2DLaw())
        _assign_analytic_stresses(
            source, _constant_function([100.0, -50.0, 30.0]))
        path = self._state_file("state.json")
        data = _write_state(path, source, stresses=True)
        del data["integration_rules"]
        dam_initial_state_file.WriteFile(path, data)

        destination_model = KM.Model()
        _build_part(destination_model, "destination", _q4_mesh(2, 4),
                    SMA.LinearElasticPlaneStrain2DLaw())
        with self.assertRaisesRegex(RuntimeError, "integration rules"):
            _import(destination_model, "destination", path, stresses=True)

    def test_matching_mesh_default_rejects_nonmatching_meshes(self):
        source_model = KM.Model()
        source = _build_part(
            source_model, "source", _q4_mesh(1, 2),
            SMA.LinearElasticPlaneStrain2DLaw())
        _assign_analytic_stresses(
            source, _constant_function([100.0, -50.0, 30.0]))
        path = self._state_file("state.json")
        _write_state(path, source, stresses=True)

        destination_model = KM.Model()
        _build_part(destination_model, "destination", _q4_mesh(2, 4),
                    SMA.LinearElasticPlaneStrain2DLaw())
        settings = KM.Parameters(json.dumps({
            "model_part_name": "destination",
            "input_file_name": str(path),
            "import_stresses": True,
            "echo_level": 0}))
        with self.assertRaisesRegex(RuntimeError, "node counts differ"):
            ImportDamInitialStateProcess(
                destination_model, settings).PrevalidateInitialState()

    def test_unsupported_transfer_method_is_rejected(self):
        model = KM.Model()
        _build_part(model, "destination", _q4_mesh(1, 2),
                    SMA.LinearElasticPlaneStrain2DLaw())
        settings = KM.Parameters(json.dumps({
            "model_part_name": "destination",
            "input_file_name": "unused.json",
            "transfer_method": "not_a_method"}))
        with self.assertRaisesRegex(Exception, "unsupported transfer_method"):
            ImportDamInitialStateProcess(model, settings)


if __name__ == "__main__":
    KratosUnittest.main()
