"""Spatial mapping of Dam initial-state fields between nonmatching 2D meshes.

The source state is the existing JSON v2 file. Source geometry is rebuilt from
the stored nodes, connectivity, geometry types and integration rules. The Core
point locator is used through the Dam bounding-box search configuration, so
higher-order geometries that do not implement Geometry::HasIntersection are
supported.

Stresses are reconstructed element-locally from the source Gauss-point values
with a small least-squares polynomial fit. No nodal smoothing is performed, so
discontinuities at material interfaces are preserved. Nodal temperatures are
interpolated with the source element shape functions.
"""
import json
import math

import numpy as np

import KratosMultiphysics as KM
import KratosMultiphysics.DamApplication as KratosDam
import KratosMultiphysics.StructuralMechanicsApplication as SMA

from KratosMultiphysics.DamApplication import dam_initial_state_file


SUPPORTED_2D_GEOMETRIES = {
    16,  # Triangle2D3
    17,  # Triangle2D6
    8,   # Quadrilateral2D4
    9,   # Quadrilateral2D8
    10,  # Quadrilateral2D9
}

GEOMETRY_TYPE_TO_NODES = {
    16: 3,
    17: 6,
    8: 4,
    9: 8,
    10: 9,
}

GEOMETRY_TYPE_TO_ELEMENT = {
    16: "SmallDisplacementElement2D3N",
    17: "SmallDisplacementElement2D6N",
    8: "SmallDisplacementElement2D4N",
    9: "SmallDisplacementElement2D8N",
    10: "SmallDisplacementElement2D9N",
}

SUPPORTED_2D_LAWS = {
    "LinearElasticPlaneStrain2DLaw",
    "LinearElasticPlaneStress2DLaw",
    "ThermalLinearElastic2DPlaneStrain",
    "ThermalLinearElastic2DPlaneStress",
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
    if name == "ThermalLinearElastic2DPlaneStrain":
        return KratosDam.ThermalLinearElastic2DPlaneStrain()
    if name == "ThermalLinearElastic2DPlaneStress":
        return KratosDam.ThermalLinearElastic2DPlaneStress()
    raise RuntimeError(
        "Spatial mapping supports only 2D linear elastic and thermal linear "
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
    return json.dumps(normalized, sort_keys=True)


def _source_material_key(record):
    return _canonical_material(record)


def _destination_material_record(element):
    properties = element.Properties
    if not properties.Has(KM.CONSTITUTIVE_LAW):
        raise RuntimeError(
            "Destination element {} has no constitutive law.".format(element.Id))
    law = properties[KM.CONSTITUTIVE_LAW]
    for variable in (KM.YOUNG_MODULUS, KM.POISSON_RATIO, KM.DENSITY, KM.THICKNESS):
        if not properties.Has(variable):
            raise RuntimeError(
                "Destination element {} does not define {}, which is required "
                "for 2D spatial mapping.".format(element.Id, variable.Name()))
    record = {
        "law": type(law).__name__,
        "young_modulus": float(properties[KM.YOUNG_MODULUS]),
        "poisson_ratio": float(properties[KM.POISSON_RATIO]),
        "density": float(properties[KM.DENSITY]),
        "strain_size": int(law.GetStrainSize()),
        "working_dimension": int(law.WorkingSpaceDimension()),
        "thickness": float(properties[KM.THICKNESS]),
    }
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
        self.element_records = {
            element["id"]: element for element in data["elements"]}
        self.stress_by_element_id = {
            int(key): values
            for key, values in data.get("initial_stresses", {}).items()}
        self.model = KM.Model()
        self.root = self.model.CreateModelPart("dam_initial_state_source", 2)
        self.root.ProcessInfo[KM.DOMAIN_SIZE] = 2

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
            locator = KratosDam.DamInitialStateSpatialSearch2D(
                self.groups[material_key]["part"])
            locator.UpdateSearchDatabase()
            self.locators[material_key] = locator
        return self.locators[material_key]

    def reconstruction_operator(self, rule_key):
        if rule_key not in self.reconstruction_cache:
            points = np.asarray(
                self.data["integration_rules"][rule_key], dtype=float)
            self.reconstruction_cache[rule_key] = _build_reconstruction_operator(
                points)
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
            term(local_coordinates[0], local_coordinates[1])
            for term in terms])
        return basis @ coefficients


def _term_constant(x, y):
    return 1.0


def _term_x(x, y):
    return x


def _term_y(x, y):
    return y


def _term_xy(x, y):
    return x * y


def _term_xx(x, y):
    return x * x


def _term_yy(x, y):
    return y * y


def _term_xxy(x, y):
    return x * x * y


def _term_xyy(x, y):
    return x * y * y


def _term_xxyy(x, y):
    return x * x * y * y


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


def _build_reconstruction_operator(points):
    """Builds a rank-checked local least-squares reconstruction operator.

    The number of polynomial terms is never larger than the number of source
    Gauss points. This is essential: an underdetermined fit does not even
    reproduce a constant field. Terms are added greedily while they increase
    the rank and keep the Vandermonde matrix reasonably conditioned.
    """
    number_of_points = len(points)
    selected = []
    for term in _TERMS_2D:
        if len(selected) >= number_of_points:
            break
        trial = selected + [term]
        vander = np.asarray([
            [candidate(point[0], point[1]) for candidate in trial]
            for point in points])
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
        [term(point[0], point[1]) for term in selected]
        for point in points])
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
    if data["dimension"] != 2:
        raise RuntimeError(
            "Spatial mapping currently supports 2D source files only; "
            "the file dimension is {}.".format(data["dimension"]))

    for element in data["elements"]:
        geometry_type = element["geometry_type"]
        if geometry_type not in SUPPORTED_2D_GEOMETRIES:
            raise RuntimeError(
                "Source element {} has unsupported 2D geometry type {}."
                .format(element["id"], geometry_type))
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

    for material in data["materials"]:
        if material["law"] not in SUPPORTED_2D_LAWS:
            raise RuntimeError(
                "Source material law '{}' is not supported by 2D spatial "
                "mapping.".format(material["law"]))
        if material["strain_size"] != 3 or material["working_dimension"] != 2:
            raise RuntimeError(
                "Source material '{}' is not a 2D material."
                .format(material["law"]))


def _validate_destination_metadata(part):
    for element in part.Elements:
        geometry_type = int(element.GetGeometry().GetGeometryType())
        if geometry_type not in SUPPORTED_2D_GEOMETRIES:
            raise RuntimeError(
                "Destination element {} has unsupported 2D geometry type {}."
                .format(element.Id, geometry_type))
        if not element.Properties.Has(KM.CONSTITUTIVE_LAW):
            raise RuntimeError(
                "Destination element {} has no constitutive law."
                .format(element.Id))
        law = element.Properties[KM.CONSTITUTIVE_LAW]
        if type(law).__name__ not in SUPPORTED_2D_LAWS:
            raise RuntimeError(
                "Destination element {} uses unsupported constitutive law '{}' "
                "for 2D spatial mapping.".format(
                    element.Id, type(law).__name__))
        if law.GetStrainSize() != 3 or law.WorkingSpaceDimension() != 2:
            raise RuntimeError(
                "Destination element {} is not a 2D material."
                .format(element.Id))
        if not element.Properties.Has(KM.THICKNESS):
            raise RuntimeError(
                "Destination element {} does not define THICKNESS, which is "
                "required for 2D spatial mapping.".format(element.Id))


def _bounds(nodes):
    values = np.asarray(nodes, dtype=float)
    return values.min(axis=0), values.max(axis=0)


def _validate_domain_overlap(part, data):
    source_nodes = [[node[1], node[2], node[3]] for node in data["nodes"]]
    destination_nodes = [[node.X0, node.Y0, node.Z0] for node in part.Nodes]
    source_low, source_high = _bounds(source_nodes)
    destination_low, destination_high = _bounds(destination_nodes)
    scale = max(
        float(np.max(source_high[:2] - source_low[:2])),
        float(np.max(destination_high[:2] - destination_low[:2])), 1.0)
    tolerance = 1.e-9 * scale
    for axis in (0, 1):
        if (destination_high[axis] < source_low[axis] - tolerance or
                destination_low[axis] > source_high[axis] + tolerance):
            raise RuntimeError(
                "Source and destination domains do not overlap in the {} "
                "coordinate.".format("xy"[axis]))


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
        material_keys = []
        for _, material_key in assignments:
            if material_key not in material_keys:
                material_keys.append(material_key)
        for material_key in material_keys:
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


def _HasDisplacementDofs(node):
    return (node.HasDofFor(KM.DISPLACEMENT_X) and
            node.HasDofFor(KM.DISPLACEMENT_Y))


def _IsMechanicalEntity(entity):
    return all(_HasDisplacementDofs(node) for node in entity.GetGeometry())


def _EntityResidual(entity, process_info):
    left_hand_side, right_hand_side = KM.Matrix(), KM.Vector()
    entity.CalculateLocalSystem(left_hand_side, right_hand_side, process_info)
    if right_hand_side.Size() == 0:
        entity.CalculateRightHandSide(right_hand_side, process_info)
    return right_hand_side


def EvaluateMechanicalResidual(part, force_scale=None):
    """Evaluates a 2D mechanical residual diagnostic at the current state.

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
            scale = max(1.0, abs(node.X0), abs(node.Y0), abs(other.X0), abs(other.Y0))
            if (abs(node.X0 - other.X0) > 1.e-12 * scale or
                    abs(node.Y0 - other.Y0) > 1.e-12 * scale or
                    node.IsFixed(KM.DISPLACEMENT_X) !=
                    other.IsFixed(KM.DISPLACEMENT_X) or
                    node.IsFixed(KM.DISPLACEMENT_Y) !=
                    other.IsFixed(KM.DISPLACEMENT_Y)):
                raise RuntimeError(
                    "EvaluateMechanicalResidual: node {} has inconsistent "
                    "coordinates or fixity across model parts.".format(node.Id))

    mechanical_nodes = [
        node for node in nodes_by_id.values() if _HasDisplacementDofs(node)]
    if not mechanical_nodes:
        raise RuntimeError(
            "EvaluateMechanicalResidual: no mechanical displacement nodes found.")
    index = {node.Id: i for i, node in enumerate(mechanical_nodes)}
    residual = [0.0] * (2 * len(mechanical_nodes))

    for element in elements.values():
        if not _IsMechanicalEntity(element):
            continue
        right_hand_side = _EntityResidual(element, part.ProcessInfo)
        if right_hand_side.Size() < 2 * element.GetGeometry().PointsNumber():
            continue
        for i, node in enumerate(element.GetGeometry()):
            residual[2 * index[node.Id]] += right_hand_side[2 * i]
            residual[2 * index[node.Id] + 1] += right_hand_side[2 * i + 1]

    for condition in conditions.values():
        if not _IsMechanicalEntity(condition):
            continue
        right_hand_side = _EntityResidual(condition, part.ProcessInfo)
        if right_hand_side.Size() < 2 * condition.GetGeometry().PointsNumber():
            continue
        for i, node in enumerate(condition.GetGeometry()):
            residual[2 * index[node.Id]] += right_hand_side[2 * i]
            residual[2 * index[node.Id] + 1] += right_hand_side[2 * i + 1]

    free_residual = []
    reaction = [0.0, 0.0]
    for i, node in enumerate(mechanical_nodes):
        if not node.IsFixed(KM.DISPLACEMENT_X):
            free_residual.append(residual[2 * i])
        if not node.IsFixed(KM.DISPLACEMENT_Y):
            free_residual.append(residual[2 * i + 1])
        if node.SolutionStepsDataHas(KM.REACTION):
            value = node.GetSolutionStepValue(KM.REACTION)
            reaction[0] += value[0]
            reaction[1] += value[1]

    free_residual_norm = math.sqrt(
        sum(value * value for value in free_residual))
    result = {
        "free_residual": free_residual,
        "max_free_residual": max(
            (abs(value) for value in free_residual), default=0.0),
        "free_residual_norm": free_residual_norm,
        "reaction_sum": reaction,
        "residual_by_node": {
            node.Id: [residual[2 * i], residual[2 * i + 1]]
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
