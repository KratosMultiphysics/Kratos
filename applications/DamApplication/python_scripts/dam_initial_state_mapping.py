"""Spatial mapping of Dam initial-state fields between nonmatching 2D or 3D meshes.

The source state is the existing JSON v2 file. Source geometry is rebuilt from
the stored nodes, connectivity, geometry types and integration rules. The Core
point locator is used through the Dam bounding-box search configuration, so
higher-order geometries that do not implement Geometry::HasIntersection are
supported.

The same algorithm is used in 2D and 3D. The only dimension-specific parts are
the geometry tables, the polynomial reconstruction basis, the locator wrapper
and the number of displacement components in the equilibrium diagnostic.

Stresses are reconstructed element-locally from the source Gauss-point values
with a small rank-checked least-squares polynomial fit. No nodal smoothing is
performed, so discontinuities at material interfaces are preserved. Nodal
temperatures are interpolated with the source element shape functions.
"""
import json
import math

import numpy as np

import KratosMultiphysics as KM
import KratosMultiphysics.DamApplication as KratosDam
import KratosMultiphysics.StructuralMechanicsApplication as SMA

from KratosMultiphysics.DamApplication import dam_initial_state_file


SUPPORTED_GEOMETRIES = {
    2: {
        16,  # Triangle2D3
        17,  # Triangle2D6
        8,   # Quadrilateral2D4
        9,   # Quadrilateral2D8
        10,  # Quadrilateral2D9
    },
    3: {
        15,  # Tetrahedra3D4
        14,  # Tetrahedra3D10
        5,   # Prism3D6
        4,   # Prism3D15
        3,   # Hexahedra3D8
        1,   # Hexahedra3D20
        2,   # Hexahedra3D27
    },
}

GEOMETRY_TYPE_TO_NODES = {
    16: 3,
    17: 6,
    8: 4,
    9: 8,
    10: 9,
    15: 4,
    14: 10,
    5: 6,
    4: 15,
    3: 8,
    1: 20,
    2: 27,
}

GEOMETRY_TYPE_TO_ELEMENT = {
    16: "SmallDisplacementElement2D3N",
    17: "SmallDisplacementElement2D6N",
    8: "SmallDisplacementElement2D4N",
    9: "SmallDisplacementElement2D8N",
    10: "SmallDisplacementElement2D9N",
    15: "SmallDisplacementElement3D4N",
    14: "SmallDisplacementElement3D10N",
    5: "SmallDisplacementElement3D6N",
    4: "SmallDisplacementElement3D15N",
    3: "SmallDisplacementElement3D8N",
    1: "SmallDisplacementElement3D20N",
    2: "SmallDisplacementElement3D27N",
}

SUPPORTED_LAWS = {
    "LinearElasticPlaneStrain2DLaw",
    "LinearElasticPlaneStress2DLaw",
    "LinearElastic3DLaw",
    "ThermalLinearElastic2DPlaneStrain",
    "ThermalLinearElastic2DPlaneStress",
    "ThermalLinearElastic3DLaw",
}

SUPPORTED_INTEGRATION_METHODS = {
    "GI_GAUSS_1", "GI_GAUSS_2", "GI_GAUSS_3", "GI_GAUSS_4", "GI_GAUSS_5",
}

SEARCH_TOLERANCE = 1.e-9
MAX_RECONSTRUCTION_CONDITION = 1.e12
MAX_REPORTED_MISSING_POINTS = 5


def _law_instance(name):
    if name == "LinearElasticPlaneStrain2DLaw":
        return SMA.LinearElasticPlaneStrain2DLaw()
    if name == "LinearElasticPlaneStress2DLaw":
        return SMA.LinearElasticPlaneStress2DLaw()
    if name == "LinearElastic3DLaw":
        return SMA.LinearElastic3DLaw()
    if name == "ThermalLinearElastic2DPlaneStrain":
        return KratosDam.ThermalLinearElastic2DPlaneStrain()
    if name == "ThermalLinearElastic2DPlaneStress":
        return KratosDam.ThermalLinearElastic2DPlaneStress()
    if name == "ThermalLinearElastic3DLaw":
        return KratosDam.ThermalLinearElastic3DLaw()
    raise RuntimeError(
        "Spatial mapping supports only linear elastic and thermal linear "
        "elastic laws; '{}' is not supported.".format(name))


def _canonical_material(record):
    normalized = {}
    for key in sorted(record):
        value = record[key]
        if isinstance(value, float):
            normalized[key] = float("{:.12g}".format(value))
        else:
            normalized[key] = value
    # THERMAL_EXPANSION does not participate in a non-thermal elastic law.
    if normalized.get("law") not in dam_initial_state_file.THERMAL_LAWS:
        normalized.pop("thermal_expansion", None)
    # THICKNESS is only a 2D section property. It must not prevent a match when
    # a 3D destination happens to carry a legacy THICKNESS value.
    if normalized.get("working_dimension") != 2:
        normalized.pop("thickness", None)
    return json.dumps(normalized, sort_keys=True)


def _source_material_key(record):
    return _canonical_material(record)


def _destination_material_record(element):
    properties = element.Properties
    if not properties.Has(KM.CONSTITUTIVE_LAW):
        raise RuntimeError(
            "Destination element {} has no constitutive law.".format(element.Id))
    law = properties[KM.CONSTITUTIVE_LAW]
    working_dimension = int(law.WorkingSpaceDimension())
    required = [KM.YOUNG_MODULUS, KM.POISSON_RATIO, KM.DENSITY]
    if working_dimension == 2:
        required.append(KM.THICKNESS)
    for variable in required:
        if not properties.Has(variable):
            raise RuntimeError(
                "Destination element {} does not define {}, which is required "
                "for {} spatial mapping.".format(
                    element.Id, variable.Name(), working_dimension))
    record = {
        "law": type(law).__name__,
        "young_modulus": float(properties[KM.YOUNG_MODULUS]),
        "poisson_ratio": float(properties[KM.POISSON_RATIO]),
        "density": float(properties[KM.DENSITY]),
        "strain_size": int(law.GetStrainSize()),
        "working_dimension": working_dimension,
    }
    if working_dimension == 2:
        record["thickness"] = float(properties[KM.THICKNESS])
    if type(law).__name__ in dam_initial_state_file.THERMAL_LAWS:
        if not properties.Has(KratosDam.THERMAL_EXPANSION):
            raise RuntimeError(
                "Thermal destination element {} does not define "
                "THERMAL_EXPANSION.".format(element.Id))
        record["thermal_expansion"] = float(
            properties[KratosDam.THERMAL_EXPANSION])
    elif properties.Has(KratosDam.THERMAL_EXPANSION):
        record["thermal_expansion"] = float(
            properties[KratosDam.THERMAL_EXPANSION])
    return record


def _destination_material_key(element):
    return _canonical_material(_destination_material_record(element))


class _SourceGeometry:
    """Temporary source model part reconstructed from the JSON v2 metadata."""

    def __init__(self, data):
        self.data = data
        self.dimension = int(data["dimension"])
        self.element_records = {
            element["id"]: element for element in data["elements"]}
        self.stress_by_element_id = {
            int(key): values
            for key, values in data.get("initial_stresses", {}).items()}
        self.model = KM.Model()
        self.root = self.model.CreateModelPart("dam_initial_state_source", 2)
        self.root.ProcessInfo[KM.DOMAIN_SIZE] = self.dimension

        for node in data["nodes"]:
            self.root.CreateNewNode(node[0], node[1], node[2], node[3])

        grouped = {}
        for element in data["elements"]:
            material = data["materials"][element["material"]]
            key = _source_material_key(material)
            grouped.setdefault(key, {
                "material": material,
                "elements": [],
            })["elements"].append(element)

        self.groups = grouped
        self.group_index = {}
        self.locators = {}
        for index, (key, group) in enumerate(self.groups.items()):
            self.group_index[key] = index
            part = self.root.CreateSubModelPart("material_group_{}".format(index))
            node_ids = sorted({
                node_id
                for element in group["elements"]
                for node_id in element["connectivity"]})
            part.AddNodes(node_ids)

            properties = part.CreateNewProperties(index + 1)
            material = group["material"]
            properties[KM.YOUNG_MODULUS] = material["young_modulus"]
            properties[KM.POISSON_RATIO] = material["poisson_ratio"]
            properties[KM.DENSITY] = material["density"]
            if "thickness" in material:
                properties[KM.THICKNESS] = material["thickness"]
            if "thermal_expansion" in material:
                properties[KratosDam.THERMAL_EXPANSION] = material["thermal_expansion"]
            properties[KM.CONSTITUTIVE_LAW] = _law_instance(material["law"])

            for element in group["elements"]:
                # These temporary elements provide geometry and integration-rule
                # identity only. Their constitutive state is never used; source
                # stresses come from the JSON file.
                element_name = GEOMETRY_TYPE_TO_ELEMENT[element["geometry_type"]]
                part.CreateNewElement(
                    element_name, element["id"], element["connectivity"],
                    properties)
            group["part"] = part

        self.reconstruction_cache = {}
        self.coefficient_cache = {}
        self.group_key_by_element_id = {
            element["id"]: key
            for key, group in self.groups.items()
            for element in group["elements"]}

    def get_locator(self, material_key):
        if material_key not in self.groups:
            raise RuntimeError(
                "The source file has no material compatible with the "
                "destination material '{}'.".format(material_key))
        if material_key not in self.locators:
            if self.dimension == 3:
                locator_class = KratosDam.DamInitialStateSpatialSearch3D
            else:
                locator_class = KratosDam.DamInitialStateSpatialSearch2D
            locator = locator_class(self.groups[material_key]["part"])
            locator.UpdateSearchDatabase()
            self.locators[material_key] = locator
        return self.locators[material_key]

    def reconstruction_operator(self, rule_key):
        if rule_key not in self.reconstruction_cache:
            points = np.asarray(
                self.data["integration_rules"][rule_key], dtype=float)
            self.reconstruction_cache[rule_key] = \
                _build_reconstruction_operator(points, self.dimension)
        return self.reconstruction_cache[rule_key]

    def stress_at(self, element_record, local_coordinates):
        rule_key = element_record["integration_rule"]
        terms, operator = self.reconstruction_operator(rule_key)
        element_id = element_record["id"]
        if element_id not in self.coefficient_cache:
            values = np.asarray(
                self.stress_by_element_id[element_id], dtype=float)
            self.coefficient_cache[element_id] = operator @ values
        coefficients = self.coefficient_cache[element_id]
        basis = np.asarray([
            term(*local_coordinates[:self.dimension])
            for term in terms])
        return basis @ coefficients


def _term_constant(x, y, z=0.0):
    return 1.0


def _term_x(x, y, z=0.0):
    return x


def _term_y(x, y, z=0.0):
    return y


def _term_z(x, y, z=0.0):
    return z


def _term_xy(x, y, z=0.0):
    return x * y


def _term_xz(x, y, z=0.0):
    return x * z


def _term_yz(x, y, z=0.0):
    return y * z


def _term_xyz(x, y, z=0.0):
    return x * y * z


def _term_xx(x, y, z=0.0):
    return x * x


def _term_yy(x, y, z=0.0):
    return y * y


def _term_zz(x, y, z=0.0):
    return z * z


def _term_xxy(x, y, z=0.0):
    return x * x * y


def _term_xxz(x, y, z=0.0):
    return x * x * z


def _term_xyy(x, y, z=0.0):
    return x * y * y


def _term_yyz(x, y, z=0.0):
    return y * y * z


def _term_xzz(x, y, z=0.0):
    return x * z * z


def _term_yzz(x, y, z=0.0):
    return y * z * z


def _term_xxyy(x, y, z=0.0):
    return x * x * y * y


def _term_xxzz(x, y, z=0.0):
    return x * x * z * z


def _term_yyzz(x, y, z=0.0):
    return y * y * z * z


def _term_xxyz(x, y, z=0.0):
    return x * x * y * z


def _term_xyyz(x, y, z=0.0):
    return x * y * y * z


def _term_xyzz(x, y, z=0.0):
    return x * y * z * z


def _term_xxyyz(x, y, z=0.0):
    return x * x * y * y * z


def _term_xxyzz(x, y, z=0.0):
    return x * x * y * z * z


def _term_xyyzz(x, y, z=0.0):
    return x * y * y * z * z


def _term_xxyyzz(x, y, z=0.0):
    return x * x * y * y * z * z


# The 2D ordering is kept unchanged to preserve the validated 2D behaviour.
_TERMS_2D = (
    _term_constant,
    _term_x,
    _term_y,
    _term_xy,
    _term_xx,
    _term_yy,
    _term_xxy,
    _term_xyy,
    _term_xxyy,
)

# Hierarchical 3D ordering. Bilinear terms precede the pure quadratics and the
# trilinear term follows them, so a distorted isoparametric hexahedron can
# still reproduce a field that is affine in physical coordinates, because its
# geometric mapping is trilinear. All 27 tri-quadratic monomials are available
# for 27-point rules such as the H20/H27 default.
_TERMS_3D = (
    _term_constant,
    _term_x,
    _term_y,
    _term_z,
    _term_xy,
    _term_xz,
    _term_yz,
    _term_xyz,
    _term_xx,
    _term_yy,
    _term_zz,
    _term_xxy,
    _term_xxz,
    _term_xyy,
    _term_yyz,
    _term_xzz,
    _term_yzz,
    _term_xxyy,
    _term_xxzz,
    _term_yyzz,
    _term_xxyz,
    _term_xyyz,
    _term_xyzz,
    _term_xxyyz,
    _term_xxyzz,
    _term_xyyzz,
    _term_xxyyzz,
)


def _build_reconstruction_operator(points, dimension):
    """Builds a rank-checked local least-squares reconstruction operator.

    The number of polynomial terms is never larger than the number of source
    Gauss points. This is essential: an underdetermined fit does not even
    reproduce a constant field. Terms are added greedily while they increase
    the rank and keep the Vandermonde matrix reasonably conditioned.

    The same procedure is used in 2D and 3D; only the ordered term set and the
    number of local coordinates differ.
    """
    terms = _TERMS_2D if dimension == 2 else _TERMS_3D
    coordinates = [point[:dimension] for point in points]
    number_of_points = len(coordinates)
    selected = []
    for term in terms:
        if len(selected) >= number_of_points:
            break
        trial = selected + [term]
        vander = np.asarray([
            [candidate(*point) for candidate in trial]
            for point in coordinates])
        if np.linalg.matrix_rank(vander) < len(trial):
            continue
        if np.linalg.cond(vander) > MAX_RECONSTRUCTION_CONDITION:
            continue
        selected = trial
    if not selected:
        raise RuntimeError(
            "Unable to build a stress reconstruction operator from the source "
            "integration points.")
    vander = np.asarray([
        [term(*point) for term in selected]
        for point in coordinates])
    operator = np.linalg.pinv(vander, rcond=1.e-12)
    return tuple(selected), operator


def _destination_gauss_point_coordinates(element):
    method_name = KratosDam.DamInitialStateTransferUtility.GetIntegrationMethodName(
        element)
    if method_name == "UNSUPPORTED":
        raise RuntimeError(
            "Destination element {} uses an unsupported integration method."
            .format(element.Id))
    local_points = \
        KratosDam.DamInitialStateTransferUtility.GetIntegrationPointCoordinates(
            element)
    geometry = element.GetGeometry()
    coordinates = []
    for local in local_points:
        point = geometry.GlobalCoordinates(KM.Array3(local))
        coordinates.append([float(point[0]), float(point[1]), float(point[2])])
    return coordinates


def _validate_source_metadata(data):
    dimension = int(data["dimension"])
    if dimension not in SUPPORTED_GEOMETRIES:
        raise RuntimeError(
            "Spatial mapping supports 2D and 3D source files only; "
            "the file dimension is {}.".format(dimension))

    supported_geometries = SUPPORTED_GEOMETRIES[dimension]
    for element in data["elements"]:
        geometry_type = element["geometry_type"]
        if geometry_type not in supported_geometries:
            raise RuntimeError(
                "Source element {} has unsupported {}-dimensional geometry "
                "type {}.".format(element["id"], dimension, geometry_type))
        if len(element["connectivity"]) != GEOMETRY_TYPE_TO_NODES[geometry_type]:
            raise RuntimeError(
                "Source element {} has {} connectivity entries but geometry "
                "type {} requires {}.".format(
                    element["id"], len(element["connectivity"]), geometry_type,
                    GEOMETRY_TYPE_TO_NODES[geometry_type]))
        method = element["integration_method"]
        if method not in SUPPORTED_INTEGRATION_METHODS:
            raise RuntimeError(
                "Source element {} has unsupported integration method '{}'."
                .format(element["id"], method))
        expected_rule = "{}|{}".format(geometry_type, method)
        if element["integration_rule"] != expected_rule:
            raise RuntimeError(
                "Source element {} has an integration rule that does not match "
                "its geometry and method.".format(element["id"]))
        points = data["integration_rules"][element["integration_rule"]]
        if not points:
            raise RuntimeError(
                "Source element {} has an empty integration rule."
                .format(element["id"]))
        if "initial_stresses" in data:
            stresses = data["initial_stresses"].get(str(element["id"]))
            if stresses is not None and len(points) != len(stresses):
                raise RuntimeError(
                    "Source element {} has {} integration points but {} stress "
                    "vectors.".format(
                        element["id"], len(points), len(stresses)))

    expected_strain_size = 3 if dimension == 2 else 6
    for material in data["materials"]:
        if material["law"] not in SUPPORTED_LAWS:
            raise RuntimeError(
                "Source material law '{}' is not supported by spatial "
                "mapping.".format(material["law"]))
        if (material["strain_size"] != expected_strain_size or
                material["working_dimension"] != dimension):
            raise RuntimeError(
                "Source material '{}' is not a {}-dimensional material."
                .format(material["law"], dimension))


def _validate_destination_metadata(part):
    dimension = int(part.ProcessInfo[KM.DOMAIN_SIZE])
    supported_geometries = SUPPORTED_GEOMETRIES[dimension]
    expected_strain_size = 3 if dimension == 2 else 6
    for element in part.Elements:
        geometry_type = int(element.GetGeometry().GetGeometryType())
        if geometry_type not in supported_geometries:
            raise RuntimeError(
                "Destination element {} has unsupported {}-dimensional "
                "geometry type {}.".format(
                    element.Id, dimension, geometry_type))
        if not element.Properties.Has(KM.CONSTITUTIVE_LAW):
            raise RuntimeError(
                "Destination element {} has no constitutive law."
                .format(element.Id))
        law = element.Properties[KM.CONSTITUTIVE_LAW]
        if type(law).__name__ not in SUPPORTED_LAWS:
            raise RuntimeError(
                "Destination element {} uses unsupported constitutive law '{}' "
                "for spatial mapping.".format(
                    element.Id, type(law).__name__))
        if (law.GetStrainSize() != expected_strain_size or
                law.WorkingSpaceDimension() != dimension):
            raise RuntimeError(
                "Destination element {} is not a {}-dimensional material."
                .format(element.Id, dimension))
        if dimension == 2 and not element.Properties.Has(KM.THICKNESS):
            raise RuntimeError(
                "Destination element {} does not define THICKNESS, which is "
                "required for 2D spatial mapping.".format(element.Id))


def _bounds(nodes):
    values = np.asarray(nodes, dtype=float)
    return values.min(axis=0), values.max(axis=0)


def _validate_domain_overlap(part, data):
    dimension = int(data["dimension"])
    source_nodes = [[node[1], node[2], node[3]] for node in data["nodes"]]
    destination_nodes = [[node.X0, node.Y0, node.Z0] for node in part.Nodes]
    source_low, source_high = _bounds(source_nodes)
    destination_low, destination_high = _bounds(destination_nodes)
    scale = max(
        float(np.max(source_high[:dimension] - source_low[:dimension])),
        float(np.max(destination_high[:dimension] -
                     destination_low[:dimension])), 1.0)
    tolerance = 1.e-9 * scale
    for axis in range(dimension):
        if (destination_high[axis] < source_low[axis] - tolerance or
                destination_low[axis] > source_high[axis] + tolerance):
            raise RuntimeError(
                "Source and destination domains do not overlap in the {} "
                "coordinate.".format("xyz"[axis]))


def ValidateSpatialDestination(part, data, import_stresses, import_temperature,
                               import_reference_temperature,
                               reference_temperature_policy):
    _validate_source_metadata(data)
    dam_initial_state_file.ValidateImportSelection(
        part, data, import_stresses, import_temperature,
        import_reference_temperature, reference_temperature_policy)
    dam_initial_state_file.ValidateDestinationDimension(part, data)
    _validate_destination_metadata(part)
    if part.NumberOfNodes() == 0 or part.NumberOfElements() == 0:
        raise RuntimeError(
            "Destination model part has no nodes or no elements.")
    _validate_domain_overlap(part, data)

    dam_initial_state_file.ValidateDestinationZeroDisplacement(part)
    dam_initial_state_file.ValidateThermalStressImport(
        part, data, import_stresses)


def _interpolate_nodal_value(data, field, source_element, shape_functions):
    values = data[field]
    result = 0.0
    geometry = source_element.GetGeometry()
    for index, node in enumerate(geometry):
        key = str(node.Id)
        if key not in values:
            raise RuntimeError(
                "Source field '{}' has no value for node {}."
                .format(field, node.Id))
        result += shape_functions[index] * float(values[key])
    return result


def _agree(first, second):
    tolerance = dam_initial_state_file.SHARED_NODE_TOLERANCE
    scale = max(1.0, abs(first), abs(second))
    return abs(first - second) <= tolerance * scale


def _map_nodal_fields(part, data, source, material_keys, import_temperature,
                      import_reference_temperature,
                      reference_temperature_policy):
    fields = []
    if import_temperature or reference_temperature_policy == "initial_temperature":
        fields.append("temperature")
    if import_reference_temperature:
        fields.append("reference_temperature")
    if not fields:
        return {}

    node_groups = {}
    for element in part.Elements:
        key = material_keys[element.Id]
        if key not in source.groups:
            raise RuntimeError(
                "Destination element {} material is not present in the source "
                "file.".format(element.Id))
        for node in element.GetGeometry():
            node_groups.setdefault(node.Id, []).append((element.Id, key))
    for assignments in node_groups.values():
        assignments.sort()

    mapped = {field: {} for field in fields}
    missing = []
    for node in part.Nodes:
        assignments = node_groups.get(node.Id)
        if not assignments:
            raise RuntimeError(
                "Destination node {} is not attached to any destination element."
                .format(node.Id))
        coordinates = KM.Array3([node.X0, node.Y0, node.Z0])
        found = {field: [] for field in fields}
        adjacent_material_keys = []
        for _, material_key in assignments:
            if material_key not in adjacent_material_keys:
                adjacent_material_keys.append(material_key)
        for material_key in adjacent_material_keys:
            locator = source.get_locator(material_key)
            result = locator.FindPoint(coordinates, SEARCH_TOLERANCE)
            if not result[0]:
                continue
            source_element = result[2]
            shape_functions = result[1]
            for field in fields:
                value = _interpolate_nodal_value(
                    data, field, source_element, shape_functions)
                found[field].append(value)
        if not any(found.values()):
            missing.append((node.Id, [node.X0, node.Y0, node.Z0]))
            continue
        for field in fields:
            values = found[field]
            if not values:
                missing.append((node.Id, [node.X0, node.Y0, node.Z0]))
                continue
            for other in values[1:]:
                if not _agree(values[0], other):
                    raise RuntimeError(
                        "Destination node {} receives conflicting {} values "
                        "from adjacent source material groups: {} and {}."
                        .format(node.Id, field, values[0], other))
            mapped[field][str(node.Id)] = values[0]

    if missing:
        details = ", ".join(
            "{} at ({}, {}, {})".format(node_id, *coordinates)
            for node_id, coordinates in missing[:MAX_REPORTED_MISSING_POINTS])
        raise RuntimeError(
            "{} destination node(s) are outside the source domain or could not "
            "be located. First missing: {}.".format(len(missing), details))
    return mapped


def _map_stresses(part, data, source, material_keys):
    mapped = {}
    missing = []
    for element in part.Elements:
        material_key = material_keys[element.Id]
        if material_key not in source.groups:
            raise RuntimeError(
                "Destination element {} material is not present in the source "
                "file.".format(element.Id))
        locator = source.get_locator(material_key)
        stresses = []
        for coordinates in _destination_gauss_point_coordinates(element):
            result = locator.FindPoint(
                KM.Array3(coordinates), SEARCH_TOLERANCE)
            if not result[0]:
                missing.append((element.Id, coordinates))
                continue
            source_element = result[2]
            source_record = source.element_records[source_element.Id]
            if source.group_key_by_element_id[source_element.Id] != material_key:
                raise RuntimeError(
                    "Destination element {} was mapped to a source element of "
                    "a different material.".format(element.Id))
            stresses.append([
                float(value) for value in source.stress_at(
                    source_record, result[3])])
        mapped[str(element.Id)] = stresses

    if missing:
        details = ", ".join(
            "element {} at ({}, {}, {})".format(element_id, *coordinates)
            for element_id, coordinates in missing[:MAX_REPORTED_MISSING_POINTS])
        raise RuntimeError(
            "{} destination integration point(s) are outside the source domain "
            "or could not be located. First missing: {}.".format(
                len(missing), details))
    return mapped


def _CollectModelPartTree(part):
    parts = [part]
    for name in part.GetSubModelPartNames():
        parts.extend(_CollectModelPartTree(part.GetSubModelPart(name)))
    return parts


_DISPLACEMENT_COMPONENTS = (
    KM.DISPLACEMENT_X,
    KM.DISPLACEMENT_Y,
    KM.DISPLACEMENT_Z,
)


def _MechanicalDofRows(entity, process_info, dimension=2):
    """Returns (local row, node ID, component) for displacement DOFs.

    The entity's own DOF list is authoritative. A thermal entity has only
    TEMPERATURE DOFs and therefore yields no rows, even when its nodes also
    carry displacement DOFs. A mixed DOF list contributes only its displacement
    rows, using the same local ordering as CalculateLocalSystem. Components
    outside the working dimension are ignored.
    """
    dofs = entity.GetDofList(process_info)
    component_by_name = {
        variable.Name(): component
        for component, variable in enumerate(_DISPLACEMENT_COMPONENTS)}
    rows = []
    for local_index, dof in enumerate(dofs):
        component = component_by_name.get(dof.GetVariable().Name())
        if component is not None and component < dimension:
            rows.append((local_index, dof.Id(), component))
    return rows, len(dofs)


def _EntityResidual(entity, process_info):
    left_hand_side, right_hand_side = KM.Matrix(), KM.Vector()
    entity.CalculateLocalSystem(left_hand_side, right_hand_side, process_info)
    if right_hand_side.Size() == 0:
        entity.CalculateRightHandSide(right_hand_side, process_info)
    return right_hand_side


def EvaluateMechanicalResidual(part, force_scale=None):
    """Evaluates a mechanical residual diagnostic in 2D or 3D at the current state.

    The function is read-only. It must be called after mechanical loads,
    supports and imported states have been initialized; the import process
    deliberately does not call it because process ordering is analysis-dependent.

    The supplied model part and all its sub-model-parts are traversed. Mechanical
    elements and conditions are assembled into one nodal residual by node ID.
    Conditions are included, so pressure, surface and point-load conditions are
    represented. Nodal loads that are not represented by a condition are not
    included. Master-slave constraints are not eliminated; the raw element and
    condition residual is reported. Shared interface nodes must be the same node
    objects (same ID); a shared ID with inconsistent coordinates or fixity is
    rejected.

    The assembled element/condition RHS is external minus internal force. Free
    entries measure the unbalanced force; nodal REACTION values are summed where
    available. The result is only physically meaningful when the model part
    actually represents the intended mechanical connection between its parts.
    """
    if part.ProcessInfo.Has(KM.DOMAIN_SIZE):
        dimension = int(part.ProcessInfo[KM.DOMAIN_SIZE])
    else:
        dimension = 2
    if dimension not in (2, 3):
        raise RuntimeError(
            "EvaluateMechanicalResidual: only 2D and 3D model parts are "
            "supported; DOMAIN_SIZE is {}.".format(dimension))
    components = _DISPLACEMENT_COMPONENTS[:dimension]

    parts = _CollectModelPartTree(part)
    nodes_by_id = {}
    elements = {}
    conditions = {}
    for current in parts:
        for element in current.Elements:
            elements.setdefault(element.Id, element)
        for condition in current.Conditions:
            conditions.setdefault(condition.Id, condition)
        for node in current.Nodes:
            if node.Id not in nodes_by_id:
                nodes_by_id[node.Id] = node
                continue
            other = nodes_by_id[node.Id]
            scale = max(
                1.0, abs(node.X0), abs(node.Y0), abs(node.Z0),
                abs(other.X0), abs(other.Y0), abs(other.Z0))
            inconsistent = any(
                abs(first - second) > 1.e-12 * scale
                for first, second in (
                    (node.X0, other.X0), (node.Y0, other.Y0),
                    (node.Z0, other.Z0)))
            for variable in components:
                inconsistent = inconsistent or (
                    node.IsFixed(variable) != other.IsFixed(variable))
            if inconsistent:
                raise RuntimeError(
                    "EvaluateMechanicalResidual: node {} has inconsistent "
                    "coordinates or fixity across model parts.".format(node.Id))

    mechanical_nodes = [
        node for node in nodes_by_id.values()
        if all(node.HasDofFor(variable) for variable in components)]
    if not mechanical_nodes:
        raise RuntimeError(
            "EvaluateMechanicalResidual: no mechanical displacement nodes found.")
    index = {node.Id: i for i, node in enumerate(mechanical_nodes)}
    residual = [0.0] * (dimension * len(mechanical_nodes))

    for element in elements.values():
        rows, number_of_dofs = _MechanicalDofRows(
            element, part.ProcessInfo, dimension)
        if not rows:
            continue
        right_hand_side = _EntityResidual(element, part.ProcessInfo)
        if right_hand_side.Size() == 0:
            continue
        if right_hand_side.Size() != number_of_dofs:
            raise RuntimeError(
                "EvaluateMechanicalResidual: element {} has {} DOFs but an "
                "RHS of size {}.".format(
                    element.Id, number_of_dofs, right_hand_side.Size()))
        for local_index, node_id, component in rows:
            residual[dimension * index[node_id] + component] += \
                right_hand_side[local_index]

    for condition in conditions.values():
        rows, number_of_dofs = _MechanicalDofRows(
            condition, part.ProcessInfo, dimension)
        if not rows:
            continue
        right_hand_side = _EntityResidual(condition, part.ProcessInfo)
        if right_hand_side.Size() == 0:
            continue
        if right_hand_side.Size() != number_of_dofs:
            raise RuntimeError(
                "EvaluateMechanicalResidual: condition {} has {} DOFs but an "
                "RHS of size {}.".format(
                    condition.Id, number_of_dofs, right_hand_side.Size()))
        for local_index, node_id, component in rows:
            residual[dimension * index[node_id] + component] += \
                right_hand_side[local_index]

    free_residual = []
    reaction = [0.0] * dimension
    for i, node in enumerate(mechanical_nodes):
        for component_index, variable in enumerate(components):
            if not node.IsFixed(variable):
                free_residual.append(
                    residual[dimension * i + component_index])
        if node.SolutionStepsDataHas(KM.REACTION):
            value = node.GetSolutionStepValue(KM.REACTION)
            for component_index in range(dimension):
                reaction[component_index] += value[component_index]

    free_residual_norm = math.sqrt(
        sum(value * value for value in free_residual))
    result = {
        "free_residual": free_residual,
        "max_free_residual": max(
            (abs(value) for value in free_residual), default=0.0),
        "free_residual_norm": free_residual_norm,
        "reaction_sum": reaction,
        "residual_by_node": {
            node.Id: [residual[dimension * i + component]
                      for component in range(dimension)]
            for i, node in enumerate(mechanical_nodes)},
    }
    if force_scale is not None and force_scale > 0.0:
        result["normalized_free_residual_norm"] = (
            free_residual_norm / force_scale)
        result["normalized_max_free_residual"] = (
            result["max_free_residual"] / force_scale)
    return result


def BuildSpatialMappedData(part, data, import_stresses, import_temperature,
                           import_reference_temperature,
                           reference_temperature_policy):
    """Maps the selected source fields to destination IDs without modifying the model."""
    source = _SourceGeometry(data)
    material_keys = {
        element.Id: _destination_material_key(element)
        for element in part.Elements}
    mapped = {"fields": list(data["fields"])}

    nodal = _map_nodal_fields(
        part, data, source, material_keys, import_temperature,
        import_reference_temperature, reference_temperature_policy)
    for field, values in nodal.items():
        mapped[field] = values

    if import_stresses:
        mapped["initial_stresses"] = _map_stresses(
            part, data, source, material_keys)
        if len(mapped["initial_stresses"]) != part.NumberOfElements():
            raise RuntimeError(
                "Spatial mapping did not produce stresses for every destination "
                "element.")
        KratosDam.DamInitialStateTransferUtility.ValidateInitialStresses(
            part,
            {int(key): value
             for key, value in mapped["initial_stresses"].items()},
            True)

    for field in ("temperature", "reference_temperature"):
        if field in mapped and len(mapped[field]) != part.NumberOfNodes():
            raise RuntimeError(
                "Spatial mapping did not produce '{}' for every destination "
                "node.".format(field))
    return mapped
