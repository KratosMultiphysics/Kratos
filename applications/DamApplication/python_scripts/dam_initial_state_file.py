"""Lightweight Dam initial-state file format and validation.

The file stores only initial Gauss-point stresses and optional nodal temperature
fields. It is deliberately not a Kratos restart file.
"""
import gzip
import json
import math
import os
from pathlib import Path

import KratosMultiphysics as KM
import KratosMultiphysics.DamApplication as KratosDam


FILE_FORMAT = "dam_initial_state"
FILE_VERSION = 1
SUPPORTED_LAWS = {
    "LinearElasticPlaneStrain2DLaw",
    "LinearElasticPlaneStress2DLaw",
    "LinearElastic3DLaw",
    "ThermalLinearElastic2DPlaneStrain",
    "ThermalLinearElastic2DPlaneStress",
    "ThermalLinearElastic3DLaw",
}
THERMAL_LAWS = {
    "ThermalLinearElastic2DPlaneStrain",
    "ThermalLinearElastic2DPlaneStress",
    "ThermalLinearElastic3DLaw",
}
STRESS_SOURCES = ("current_cauchy", "stored_initial_stress")


def WriteFile(file_name, data):
    path = Path(file_name)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    if str(path).endswith(".gz"):
        with gzip.open(temporary, "wt", encoding="utf-8") as stream:
            json.dump(data, stream, separators=(",", ":"), allow_nan=False)
    else:
        with open(temporary, "w", encoding="utf-8") as stream:
            json.dump(data, stream, separators=(",", ":"), allow_nan=False)
    os.replace(temporary, path)


def ReadFile(file_name):
    path = Path(file_name)
    if not path.is_file():
        raise FileNotFoundError(
            "Dam initial-state file not found: {}".format(file_name))
    if str(path).endswith(".gz"):
        with gzip.open(path, "rt", encoding="utf-8") as stream:
            data = json.load(stream)
    else:
        with open(path, "r", encoding="utf-8") as stream:
            data = json.load(stream)
    ValidateFile(data)
    return data


def _law_name(properties):
    return type(properties[KM.CONSTITUTIVE_LAW]).__name__


def _integration_info(element):
    properties = element.Properties
    order_variable = KM.KratosGlobals.GetVariable("INTEGRATION_ORDER")
    if properties.Has(order_variable):
        order = int(properties[order_variable])
        if order < 1 or order > 5:
            raise RuntimeError(
                "Unsupported INTEGRATION_ORDER {} on element {}.".format(
                    order, element.Id))
    else:
        order = 0
    method_name = KratosDam.SelectiveInitialStressUtility.GetIntegrationMethodName(element)
    if method_name == "UNSUPPORTED":
        raise RuntimeError(
            "Unsupported integration method on element {}.".format(element.Id))
    geometry_type = int(element.GetGeometry().GetGeometryType())
    key = "{}|{}".format(geometry_type, method_name)
    points = [
        list(point) for point in
        KratosDam.SelectiveInitialStressUtility.GetIntegrationPointCoordinates(element)]
    return key, order, method_name, points


def _material_record(element, name):
    properties = element.Properties
    law = properties[KM.CONSTITUTIVE_LAW]
    record = {
        "law": name,
        "young_modulus": float(properties[KM.YOUNG_MODULUS]),
        "poisson_ratio": float(properties[KM.POISSON_RATIO]),
        "density": float(properties[KM.DENSITY]),
        "strain_size": int(law.GetStrainSize()),
        "working_dimension": int(law.WorkingSpaceDimension()),
    }
    if properties.Has(KM.THICKNESS):
        record["thickness"] = float(properties[KM.THICKNESS])
    if properties.Has(KratosDam.THERMAL_EXPANSION):
        record["thermal_expansion"] = float(properties[KratosDam.THERMAL_EXPANSION])
    return record


def BuildFileData(part, save_stresses, save_temperature, save_reference_temperature,
                  stress_source="current_cauchy"):
    if not (save_stresses or save_temperature or save_reference_temperature):
        raise RuntimeError("At least one initial-state field must be selected for export.")
    if save_stresses and stress_source not in STRESS_SOURCES:
        raise RuntimeError(
            "Unsupported stress_source '{}'. Use one of {}.".format(
                stress_source, STRESS_SOURCES))
    if save_temperature and not part.HasNodalSolutionStepVariable(KM.TEMPERATURE):
        raise RuntimeError("TEMPERATURE is not defined on the exported model part.")
    if (save_reference_temperature and
            not part.HasNodalSolutionStepVariable(KratosDam.NODAL_REFERENCE_TEMPERATURE)):
        raise RuntimeError(
            "NODAL_REFERENCE_TEMPERATURE is not defined on the exported model part.")

    if not part.ProcessInfo.Has(KM.DOMAIN_SIZE):
        raise RuntimeError("DOMAIN_SIZE is not defined on the exported model part.")
    dimension = int(part.ProcessInfo[KM.DOMAIN_SIZE])
    fields = []
    if save_stresses:
        fields.append("initial_stresses")
    if save_temperature:
        fields.append("temperature")
    if save_reference_temperature:
        fields.append("reference_temperature")

    nodes = []
    for node in part.Nodes:
        nodes.append([node.Id, node.X0, node.Y0, node.Z0])

    materials = []
    material_index = {}
    elements = []
    integration_rules = {}
    for element in part.Elements:
        law_name = _law_name(element.Properties)
        if law_name not in SUPPORTED_LAWS:
            raise RuntimeError(
                "Unsupported constitutive law '{}' on element {}.".format(
                    law_name, element.Id))
        record = _material_record(element, law_name)
        if dimension == 2 and "thickness" not in record:
            raise RuntimeError(
                "2D element {} does not define THICKNESS.".format(element.Id))
        key = json.dumps(record, sort_keys=True)
        if key not in material_index:
            material_index[key] = len(materials)
            materials.append(record)
        rule_key, order, method_name, points = _integration_info(element)
        if rule_key in integration_rules:
            if len(integration_rules[rule_key]) != len(points):
                raise RuntimeError(
                    "Inconsistent integration rule for element {}.".format(element.Id))
        else:
            integration_rules[rule_key] = points
        elements.append({
            "id": element.Id,
            "element_type": type(element).__name__,
            "geometry_type": int(element.GetGeometry().GetGeometryType()),
            "integration_rule": rule_key,
            "integration_order": order,
            "integration_method": method_name,
            "connectivity": [node.Id for node in element.GetGeometry()],
            "material": material_index[key],
        })

    data = {
        "format": FILE_FORMAT,
        "version": FILE_VERSION,
        "model_part_name": part.Name,
        "dimension": dimension,
        "fields": fields,
        "nodes": nodes,
        "elements": elements,
        "materials": materials,
        "integration_rules": integration_rules,
    }
    if save_stresses:
        data["stress_source"] = stress_source
        if stress_source == "current_cauchy":
            stresses = KratosDam.SelectiveInitialStressUtility.GetCurrentCauchyStresses(part)
        else:
            stresses = KratosDam.SelectiveInitialStressUtility.GetStoredInitialStresses(part)
        data["initial_stresses"] = {
            str(element_id): [list(value) for value in values]
            for element_id, values in stresses.items()}
    if save_temperature:
        data["temperature"] = {
            str(node.Id): node.GetSolutionStepValue(KM.TEMPERATURE)
            for node in part.Nodes}
    if save_reference_temperature:
        data["reference_temperature"] = {
            str(node.Id): node.GetSolutionStepValue(KratosDam.NODAL_REFERENCE_TEMPERATURE)
            for node in part.Nodes}
    ValidateFile(data)
    return data


def _require(condition, message):
    if not condition:
        raise RuntimeError(message)


def _check_finite(value, message):
    _require(isinstance(value, (int, float)) and math.isfinite(value), message)


def _integer_keys(values, label):
    result = set()
    for key in values:
        try:
            integer_key = int(key)
        except (TypeError, ValueError):
            raise RuntimeError(
                "{} contains a non-integer ID key '{}'.".format(label, key))
        if str(integer_key) != str(key):
            raise RuntimeError(
                "{} contains a non-canonical integer ID key '{}'.".format(
                    label, key))
        result.add(integer_key)
    return result


def ValidateFile(data):
    _require(isinstance(data, dict), "Dam initial-state file must contain a JSON object.")
    _require(data.get("format") == FILE_FORMAT,
             "Unknown Dam initial-state format '{}'.".format(data.get("format")))
    _require(data.get("version") == FILE_VERSION,
             "Unsupported Dam initial-state version '{}'.".format(data.get("version")))
    _require(data.get("dimension") in (2, 3),
             "Dam initial-state dimension must be 2 or 3.")

    fields = data.get("fields")
    _require(isinstance(fields, list) and fields,
             "Dam initial-state 'fields' must be a nonempty list.")
    allowed_fields = {"initial_stresses", "temperature", "reference_temperature"}
    _require(set(fields).issubset(allowed_fields),
             "Dam initial-state file contains an unsupported field selection.")

    nodes = data.get("nodes")
    _require(isinstance(nodes, list) and nodes, "Dam initial-state file has no nodes.")
    node_ids = set()
    for node in nodes:
        _require(isinstance(node, list) and len(node) == 4,
                 "Every node entry must be [id, x, y, z].")
        _require(node[0] not in node_ids, "Duplicate node ID in Dam initial-state file.")
        node_ids.add(node[0])
        for coordinate in node[1:]:
            _check_finite(coordinate, "Non-finite reference coordinate in initial-state file.")

    materials = data.get("materials")
    _require(isinstance(materials, list) and materials,
             "Dam initial-state file has no materials.")
    for material in materials:
        _require(isinstance(material, dict), "Material entry must be an object.")
        for key in ("law", "young_modulus", "poisson_ratio", "density",
                    "strain_size", "working_dimension"):
            _require(key in material, "Material entry is missing '{}'.".format(key))
        _require(material["law"] in SUPPORTED_LAWS,
                 "Unsupported constitutive law '{}' in initial-state file.".format(
                     material["law"]))
        _require(material["strain_size"] in (3, 6),
                 "Invalid strain size in initial-state file.")
        _require(material["working_dimension"] in (2, 3),
                 "Invalid working dimension in initial-state file.")
        for key in ("young_modulus", "poisson_ratio", "density"):
            _check_finite(float(material[key]), "Invalid material value '{}'.".format(key))
        if "thickness" in material:
            _check_finite(float(material["thickness"]), "Invalid thickness.")
        if "thermal_expansion" in material:
            _check_finite(float(material["thermal_expansion"]), "Invalid thermal expansion.")
        if material["law"] in THERMAL_LAWS:
            _require("thermal_expansion" in material,
                     "Thermal law material is missing thermal_expansion.")

    rules = data.get("integration_rules")
    _require(isinstance(rules, dict) and rules,
             "Dam initial-state file has no integration rules.")
    for key, points in rules.items():
        _require(isinstance(points, list) and points,
                 "Integration rule '{}' is empty.".format(key))
        for point in points:
            _require(isinstance(point, list) and len(point) == 3,
                     "Integration point in rule '{}' must have three coordinates.".format(key))
            for coordinate in point:
                _check_finite(coordinate, "Non-finite integration-point coordinate.")

    elements = data.get("elements")
    _require(isinstance(elements, list) and elements,
             "Dam initial-state file has no elements.")
    element_ids = set()
    element_points = {}
    element_material = {}
    for element in elements:
        _require(isinstance(element, dict), "Element entry must be an object.")
        for key in ("id", "element_type", "geometry_type", "integration_rule",
                    "integration_order", "integration_method", "connectivity",
                    "material"):
            _require(key in element, "Element entry is missing '{}'.".format(key))
        element_id = element["id"]
        _require(element_id not in element_ids,
                 "Duplicate element ID in Dam initial-state file.")
        element_ids.add(element_id)
        _require(element["integration_rule"] in rules,
                 "Element {} references an unknown integration rule.".format(element_id))
        _require(0 <= element["material"] < len(materials),
                 "Element {} references an invalid material.".format(element_id))
        _require(isinstance(element["connectivity"], list) and element["connectivity"],
                 "Element {} has invalid connectivity.".format(element_id))
        for node_id in element["connectivity"]:
            _require(node_id in node_ids,
                     "Element {} references missing node {}.".format(element_id, node_id))
        element_points[element_id] = len(rules[element["integration_rule"]])
        element_material[element_id] = element["material"]

    if "initial_stresses" in fields:
        _require(data.get("stress_source") in STRESS_SOURCES,
                 "Invalid or missing stress_source.")
        stresses = data.get("initial_stresses")
        _require(isinstance(stresses, dict), "initial_stresses must be an object.")
        _require(_integer_keys(stresses, "initial_stresses") == element_ids,
                 "initial_stresses must contain exactly the file's element IDs.")
        for element in elements:
            element_id = element["id"]
            values = stresses[str(element_id)]
            _require(isinstance(values, list) and len(values) == element_points[element_id],
                     "Element {} has the wrong number of stress vectors.".format(element_id))
            expected_size = materials[element_material[element_id]]["strain_size"]
            for value in values:
                _require(isinstance(value, list) and len(value) == expected_size,
                         "Element {} has an invalid stress-vector size.".format(element_id))
                for component in value:
                    _check_finite(component, "Non-finite initial stress.")

    for field, variable_name in (
            ("temperature", "TEMPERATURE"),
            ("reference_temperature", "NODAL_REFERENCE_TEMPERATURE")):
        if field in fields:
            values = data.get(field)
            _require(isinstance(values, dict), "{} must be an object.".format(field))
            _require(_integer_keys(values, field) == node_ids,
                     "{} must contain exactly the file's node IDs.".format(field))
            for value in values.values():
                _check_finite(value, "Non-finite {} value.".format(variable_name))


def _file_node_data(data):
    return {node[0]: (node[1], node[2], node[3]) for node in data["nodes"]}


def _file_element_data(data):
    return {element["id"]: element for element in data["elements"]}


def _field_values(data, field):
    return {int(key): value for key, value in data[field].items()}


def _stress_values(data):
    return {
        int(key): [list(value) for value in values]
        for key, values in data["initial_stresses"].items()}


def _compare_value(actual, expected, label, tolerance=1.e-12):
    scale = 1.0 + abs(expected)
    if abs(actual - expected) > tolerance * scale:
        raise RuntimeError(
            "Dam initial-state mismatch for {}: {} vs {}.".format(
                label, actual, expected))


def ValidateDestination(part, data, transfer_stresses, transfer_temperature,
                        transfer_reference_temperature, reference_policy):
    if not (transfer_stresses or transfer_temperature or transfer_reference_temperature):
        raise RuntimeError("At least one initial-state field must be selected for import.")
    if reference_policy not in ("preserve", "initial_temperature"):
        raise RuntimeError(
            "Unsupported reference_temperature_policy '{}'.".format(reference_policy))
    if reference_policy == "initial_temperature":
        if not transfer_temperature:
            raise RuntimeError(
                "reference_temperature_policy='initial_temperature' requires "
                "transfer_initial_temperature=true.")
        if transfer_reference_temperature:
            raise RuntimeError(
                "Cannot both import reference_temperature and replace it with the "
                "imported initial temperature.")
    if transfer_stresses and "initial_stresses" not in data["fields"]:
        raise RuntimeError("The file does not contain initial_stresses.")
    if transfer_temperature and "temperature" not in data["fields"]:
        raise RuntimeError("The file does not contain temperature.")
    if (transfer_reference_temperature and
            "reference_temperature" not in data["fields"]):
        raise RuntimeError("The file does not contain reference_temperature.")
    if transfer_temperature and not part.HasNodalSolutionStepVariable(KM.TEMPERATURE):
        raise RuntimeError("TEMPERATURE is not defined on the destination model part.")
    if (transfer_reference_temperature and
            not part.HasNodalSolutionStepVariable(KratosDam.NODAL_REFERENCE_TEMPERATURE)):
        raise RuntimeError(
            "NODAL_REFERENCE_TEMPERATURE is not defined on the destination model part.")

    if not part.ProcessInfo.Has(KM.DOMAIN_SIZE):
        raise RuntimeError("DOMAIN_SIZE is not defined on the destination model part.")
    if int(part.ProcessInfo[KM.DOMAIN_SIZE]) != data["dimension"]:
        raise RuntimeError("Source and destination dimensions differ.")
    if part.NumberOfNodes() != len(data["nodes"]):
        raise RuntimeError("Source and destination node counts differ.")
    if part.NumberOfElements() != len(data["elements"]):
        raise RuntimeError("Source and destination element counts differ.")

    file_nodes = _file_node_data(data)
    for node in part.Nodes:
        if node.Id not in file_nodes:
            raise RuntimeError("Destination node {} is missing from the file.".format(node.Id))
        expected = file_nodes[node.Id]
        actual = (node.X0, node.Y0, node.Z0)
        for component, expected_component in zip(actual, expected):
            _compare_value(component, expected_component,
                           "reference coordinate of node {}".format(node.Id))

    file_elements = _file_element_data(data)
    for element in part.Elements:
        if element.Id not in file_elements:
            raise RuntimeError("Destination element {} is missing from the file.".format(element.Id))
        expected = file_elements[element.Id]
        if type(element).__name__ != expected["element_type"]:
            raise RuntimeError(
                "Element {} type mismatch: {} vs {}.".format(
                    element.Id, type(element).__name__, expected["element_type"]))
        if int(element.GetGeometry().GetGeometryType()) != expected["geometry_type"]:
            raise RuntimeError("Element {} geometry type mismatch.".format(element.Id))
        actual_connectivity = [node.Id for node in element.GetGeometry()]
        if actual_connectivity != expected["connectivity"]:
            raise RuntimeError("Element {} connectivity/order mismatch.".format(element.Id))
        rule_key, order, method_name, points = _integration_info(element)
        if (rule_key != expected["integration_rule"] or
                method_name != expected["integration_method"]):
            raise RuntimeError(
                "Element {} integration rule mismatch: {} vs {}.".format(
                    element.Id, rule_key, expected["integration_rule"]))
        file_points = data["integration_rules"][expected["integration_rule"]]
        if len(points) != len(file_points):
            raise RuntimeError("Element {} integration-point count mismatch.".format(element.Id))
        for actual_point, expected_point in zip(points, file_points):
            for actual_coordinate, expected_coordinate in zip(actual_point, expected_point):
                _compare_value(actual_coordinate, expected_coordinate,
                               "integration point of element {}".format(element.Id))

        material = data["materials"][expected["material"]]
        properties = element.Properties
        law = properties[KM.CONSTITUTIVE_LAW]
        actual_law = type(law).__name__
        if actual_law != material["law"]:
            raise RuntimeError(
                "Element {} constitutive law mismatch: {} vs {}.".format(
                    element.Id, actual_law, material["law"]))
        if int(law.GetStrainSize()) != material["strain_size"]:
            raise RuntimeError("Element {} strain size mismatch.".format(element.Id))
        if int(law.WorkingSpaceDimension()) != material["working_dimension"]:
            raise RuntimeError("Element {} working dimension mismatch.".format(element.Id))
        for key, variable in (("young_modulus", KM.YOUNG_MODULUS),
                              ("poisson_ratio", KM.POISSON_RATIO),
                              ("density", KM.DENSITY)):
            if not properties.Has(variable):
                raise RuntimeError(
                    "Element {} properties do not define {}.".format(
                        element.Id, variable.Name()))
            _compare_value(float(properties[variable]), float(material[key]),
                           "{} of element {}".format(variable.Name(), element.Id))
        if data["dimension"] == 2:
            if "thickness" not in material:
                raise RuntimeError("2D material entry is missing THICKNESS.")
            if not properties.Has(KM.THICKNESS):
                raise RuntimeError("Element {} does not define THICKNESS.".format(element.Id))
            _compare_value(float(properties[KM.THICKNESS]), float(material["thickness"]),
                           "THICKNESS of element {}".format(element.Id))
        elif "thickness" in material and properties.Has(KM.THICKNESS):
            _compare_value(float(properties[KM.THICKNESS]), float(material["thickness"]),
                           "THICKNESS of element {}".format(element.Id))
        if "thermal_expansion" in material:
            if not properties.Has(KratosDam.THERMAL_EXPANSION):
                raise RuntimeError(
                    "Element {} does not define THERMAL_EXPANSION.".format(element.Id))
            _compare_value(float(properties[KratosDam.THERMAL_EXPANSION]),
                           float(material["thermal_expansion"]),
                           "THERMAL_EXPANSION of element {}".format(element.Id))
        elif material["law"] in THERMAL_LAWS:
            raise RuntimeError("Thermal material entry is missing thermal_expansion.")

    if not part.HasNodalSolutionStepVariable(KM.DISPLACEMENT):
        raise RuntimeError("DISPLACEMENT is not defined on the destination model part.")
    for node in part.Nodes:
        for step in range(node.GetBufferSize()):
            displacement = node.GetSolutionStepValue(KM.DISPLACEMENT, step)
            if sum(value * value for value in displacement) > 1.e-24:
                raise RuntimeError(
                    "Destination node {} has a nonzero displacement in buffer step {}."
                    .format(node.Id, step))

    thermal_law_present = any(
        material["law"] in THERMAL_LAWS for material in data["materials"])
    if transfer_stresses and thermal_law_present:
        if not part.HasNodalSolutionStepVariable(KM.TEMPERATURE):
            raise RuntimeError(
                "TEMPERATURE is not defined on the destination model part; "
                "thermal-strain compensation for imported stresses is impossible.")
        if not part.HasNodalSolutionStepVariable(
                KratosDam.NODAL_REFERENCE_TEMPERATURE):
            raise RuntimeError(
                "NODAL_REFERENCE_TEMPERATURE is not defined on the destination "
                "model part; thermal-strain compensation for imported stresses "
                "is impossible.")

    if transfer_stresses:
        KratosDam.SelectiveInitialStressUtility.ValidateStoredInitialStresses(
            part, _stress_values(data), True)
    return data


SHARED_NODE_TOLERANCE = 1.e-9


def _shared_values_agree(first, second, tolerance=SHARED_NODE_TOLERANCE):
    scale = max(1.0, abs(first), abs(second))
    return abs(first - second) <= tolerance * scale


def _describe_import(process):
    return "'{}' (part '{}')".format(
        process.input_file_name, process.model_part.Name)


def ValidateCombinedImports(processes):
    """Validates multiple initial-state imports before any field is applied.

    Stress imports must target disjoint destination element sets. Nodal fields
    may overlap only when every shared node receives the same value (within a
    relative tolerance). All files and destinations are validated first; this
    function never modifies the model.
    """
    if not processes:
        return

    entries = []
    for process in processes:
        entries.append((process, process.PrevalidateInitialState()))

    stress_owner = {}
    temperature_values = {}
    reference_values = {}

    def check_shared_value(field, assignments, node_id, value, process):
        if node_id in assignments:
            other_value, other_process = assignments[node_id]
            if not _shared_values_agree(value, other_value):
                raise RuntimeError(
                    "Dam initial-state import conflict: shared node {} receives "
                    "{} = {} from {} and {} from {}. Shared-node imports must "
                    "agree within relative tolerance {}.".format(
                        node_id, field, other_value,
                        _describe_import(other_process), value,
                        _describe_import(process), SHARED_NODE_TOLERANCE))
        else:
            assignments[node_id] = (value, process)

    for process, data in entries:
        if process.transfer_initial_stresses:
            for element in process.model_part.Elements:
                element_id = element.Id
                if element_id in stress_owner:
                    raise RuntimeError(
                        "Dam initial-state import conflict: element {} is assigned "
                        "initial stresses by both {} and {}. Stress imports must "
                        "target disjoint element sets.".format(
                            element_id, _describe_import(stress_owner[element_id]),
                            _describe_import(process)))
                stress_owner[element_id] = process

        if process.transfer_initial_temperature:
            for node_id, value in data["temperature"].items():
                check_shared_value(
                    "TEMPERATURE", temperature_values, int(node_id),
                    float(value), process)

        if process.reference_temperature_policy == "initial_temperature":
            reference_source = data["temperature"]
        elif process.transfer_reference_temperature:
            reference_source = data["reference_temperature"]
        else:
            reference_source = None
        if reference_source is not None:
            for node_id, value in reference_source.items():
                check_shared_value(
                    "NODAL_REFERENCE_TEMPERATURE", reference_values,
                    int(node_id), float(value), process)


def ApplyTemperatures(part, data, transfer_temperature,
                      transfer_reference_temperature, reference_policy):
    temperature = _field_values(data, "temperature") if transfer_temperature else None
    reference = (_field_values(data, "reference_temperature")
                 if transfer_reference_temperature else None)
    for node in part.Nodes:
        if transfer_temperature:
            value = temperature[node.Id]
            for step in range(node.GetBufferSize()):
                node.SetSolutionStepValue(KM.TEMPERATURE, step, value)
        if reference_policy == "initial_temperature":
            value = temperature[node.Id]
            for step in range(node.GetBufferSize()):
                node.SetSolutionStepValue(
                    KratosDam.NODAL_REFERENCE_TEMPERATURE, step, value)
        elif transfer_reference_temperature:
            value = reference[node.Id]
            for step in range(node.GetBufferSize()):
                node.SetSolutionStepValue(
                    KratosDam.NODAL_REFERENCE_TEMPERATURE, step, value)


def ApplyStresses(part, data, compensate_initial_thermal_strain=True):
    KratosDam.SelectiveInitialStressUtility.AssignStoredInitialStresses(
        part, _stress_values(data), compensate_initial_thermal_strain)
