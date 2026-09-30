"""Single owner of the RomParameters.json format.

Two layouts exist:
- Version 1 (legacy): flat file mixing run-mode flags ('train_hrom', 'run_hrom', 'rom_manager'), solver settings,
  inline data ('nodal_modes', 'elements_and_weights') and HROM training settings. Numpy data files are found by naming convention.
- Version 2: a manifest describing the trained ROM as typed blocks. It contains a 'deployment' block with the default
  deployment options, a list of 'components' (one per solver), each with its 'projection' and 'manifold', an optional
  'coupling' block and a 'hyper_reduction' block. Numpy data files are referenced explicitly.
  The 'manifold' type sets the scope: 'global' (a single 'decoder') or 'local' (a 'selector' and a list of 'clusters',
  each with its own 'decoder'). The 'decoder' type sets how the reduced coordinates map to the full solution
  ('linear', 'ann_enhanced', ...). Scope and decoder are independent, so every decoder type is usable in both scopes.

Every read and write of RomParameters.json goes through this module. Files of either version can be read.
Files are still written as version 1 (see _OUTPUT_VERSION) until all the readers use the version 2 layout.
"""

import copy
import json
from pathlib import Path

ROM_PARAMETERS_VERSION = 2
_OUTPUT_VERSION = 1

# Conventional names of the numpy files stored next to a version 1 RomParameters.json
_RIGHT_BASIS_FILE = "RightBasisMatrix.npy"
_NODE_IDS_FILE = "NodeIds.npy"
_LEFT_BASIS_FILE = "LeftBasisMatrix.npy"
_HROM_FILES = {
    "element_ids" : "HROM_ElementIds.npy",
    "element_weights" : "HROM_ElementWeights.npy",
    "condition_ids" : "HROM_ConditionIds.npy",
    "condition_weights" : "HROM_ConditionWeights.npy"
}

# Values assumed by the readers when a version 1 key is missing
_V1_DEFAULTS = {
    "rom_manager" : False,
    "train_hrom" : False,
    "run_hrom" : False,
    "projection_strategy" : "galerkin",
    "assembling_strategy" : "global",
    "rom_format" : "json",
    "rom_settings" : {"rom_bns_settings" : {}},
    "hrom_settings" : {},
    "nodal_modes" : {},
    "elements_and_weights" : {}
}

_ANN_SUFFIX = "_ann"


def GetRomParametersFilePath(rom_basis_output_folder, rom_basis_output_name):
    return Path(rom_basis_output_folder) / Path(rom_basis_output_name).with_suffix(".json")


def LoadRomParameters(rom_basis_output_folder, rom_basis_output_name):
    """Returns the content of the RomParameters file in the version 2 layout, whatever the version stored."""
    with GetRomParametersFilePath(rom_basis_output_folder, rom_basis_output_name).open('r') as f:
        data = json.load(f)
    return data if IsVersion2(data) else UpgradeToVersion2(data)


def ReadRomParametersAsVersion1(rom_basis_output_folder, rom_basis_output_name):
    """Returns the content of the RomParameters file in the version 1 layout, whatever the version stored.
    Transitional: used by the readers that have not been moved to the version 2 layout yet."""
    with GetRomParametersFilePath(rom_basis_output_folder, rom_basis_output_name).open('r') as f:
        data = json.load(f)
    return DowngradeToVersion1(data) if IsVersion2(data) else data


def WriteRomParameters(rom_basis_output_folder, rom_basis_output_name, data):
    """Writes the RomParameters file. 'data' can be given in either layout; it is stored in the current output version."""
    if _OUTPUT_VERSION == 1:
        data = DowngradeToVersion1(data) if IsVersion2(data) else data
    else:
        data = data if IsVersion2(data) else UpgradeToVersion2(data)

    file_path = GetRomParametersFilePath(rom_basis_output_folder, rom_basis_output_name)
    file_path.parent.mkdir(parents=True, exist_ok=True)
    with file_path.open('w') as f:
        json.dump(data, f, indent=4)


def IsVersion2(data):
    return data.get("rom_parameters_version", 1) == 2


def UpgradeToVersion2(v1_data):
    """Converts a version 1 RomParameters dictionary into the version 2 layout. No information is lost."""
    v1 = copy.deepcopy(v1_data)
    rom_settings = v1.pop("rom_settings", {})

    # Projection: the ANN-enhanced strategies are the plain strategies with a nonlinear decoder
    strategy = v1.pop("projection_strategy", _V1_DEFAULTS["projection_strategy"])
    decoder_type = "linear"
    if strategy.endswith(_ANN_SUFFIX):
        strategy = strategy[:-len(_ANN_SUFFIX)]
        decoder_type = "ann_enhanced"
    projection = {
        "strategy" : strategy,
        "assembling_strategy" : v1.pop("assembling_strategy", _V1_DEFAULTS["assembling_strategy"]),
        "bns_settings" : rom_settings.pop("rom_bns_settings", {})
    }

    # Right basis
    rom_format = v1.pop("rom_format", _V1_DEFAULTS["rom_format"])
    basis = {"format" : rom_format}
    if "number_of_rom_dofs" in rom_settings:
        basis["number_of_modes"] = rom_settings.pop("number_of_rom_dofs")
    if rom_format == "numpy":
        basis["file"] = _RIGHT_BASIS_FILE
        basis["node_ids"] = _NODE_IDS_FILE
    nodal_modes = v1.pop("nodal_modes", {})
    if nodal_modes or rom_format == "json":
        basis["nodal_modes"] = nodal_modes

    # Left basis (Petrov-Galerkin)
    left_basis = {}
    if "petrov_galerkin_number_of_rom_dofs" in rom_settings:
        left_basis["number_of_modes"] = rom_settings.pop("petrov_galerkin_number_of_rom_dofs")
    if "petrov_galerkin_nodal_modes" in v1:
        left_basis["nodal_modes"] = v1.pop("petrov_galerkin_nodal_modes")
    if rom_format == "numpy" and left_basis.get("number_of_modes", 0) > 0:
        left_basis["file"] = _LEFT_BASIS_FILE
    if left_basis:
        left_basis["format"] = rom_format
        projection["left_basis"] = left_basis

    component = {
        "name" : "",
        "model_part_name" : "",
        "nodal_unknowns" : rom_settings.pop("nodal_unknowns", []),
        "projection" : projection,
        "manifold" : {"type" : "global", "decoder" : {"type" : decoder_type, "basis" : basis}}
    }
    if rom_settings:
        component["legacy_rom_settings"] = rom_settings

    # Hyper-reduction
    hrom_settings = v1.pop("hrom_settings", {})
    hyper_reduction = {"type" : hrom_settings.get("element_selection_type", "empirical_cubature")}
    hrom_format = hrom_settings.pop("hrom_format", None)
    if hrom_format is not None:
        hyper_reduction["format"] = hrom_format
        if hrom_format == "numpy":
            hyper_reduction.update(_HROM_FILES)
    hyper_reduction["elements_and_weights"] = v1.pop("elements_and_weights", {})
    hyper_reduction["training_settings"] = hrom_settings

    v2 = {
        "rom_parameters_version" : ROM_PARAMETERS_VERSION,
        "deployment" : {
            "use_hyper_reduction" : v1.pop("run_hrom", _V1_DEFAULTS["run_hrom"]),
            "train_hyper_reduction" : v1.pop("train_hrom", _V1_DEFAULTS["train_hrom"]), # Deprecated: hyper-reduction training no longer happens inside the simulation
            "rom_manager" : v1.pop("rom_manager", _V1_DEFAULTS["rom_manager"])
        },
        "components" : [component],
        "hyper_reduction" : hyper_reduction
    }

    # Coupled solvers sharing a single basis. The order of the solvers defines their column of HROM weights
    if "coupled_solvers" in v1:
        v2["coupling"] = {"basis_layout" : "monolithic", "solvers" : v1.pop("coupled_solvers")}

    if v1:
        v2["legacy"] = v1
    return v2


def DowngradeToVersion1(v2_data):
    """Converts a version 2 RomParameters dictionary into the version 1 layout.
    Only the version 2 content expressible in version 1 is supported (single component, global manifold, linear or ANN-enhanced decoder)."""
    v2 = copy.deepcopy(v2_data)
    components = v2["components"]
    if len(components) != 1:
        raise Exception(f"RomParameters with {len(components)} components cannot be written in the version 1 layout.")
    component = components[0]
    projection = component["projection"]
    manifold = component["manifold"]
    if manifold["type"] != "global":
        raise Exception(f"'{manifold['type']}' manifold cannot be written in the version 1 layout, which only supports 'global' manifolds.")
    decoder = manifold["decoder"]
    if decoder["type"] not in ("linear", "ann_enhanced"):
        raise Exception(f"'{decoder['type']}' decoder cannot be written in the version 1 layout, which only supports 'linear' and 'ann_enhanced' decoders.")
    basis = decoder["basis"]
    _CheckConventionalFileNames(basis, {"file" : _RIGHT_BASIS_FILE, "node_ids" : _NODE_IDS_FILE})

    strategy = projection["strategy"]
    if decoder["type"] == "ann_enhanced":
        strategy += _ANN_SUFFIX

    rom_settings = {"rom_bns_settings" : projection.get("bns_settings", {})}
    rom_settings["nodal_unknowns"] = component["nodal_unknowns"]
    if "number_of_modes" in basis:
        rom_settings["number_of_rom_dofs"] = basis["number_of_modes"]
    left_basis = projection.get("left_basis", {})
    _CheckConventionalFileNames(left_basis, {"file" : _LEFT_BASIS_FILE})
    if "number_of_modes" in left_basis:
        rom_settings["petrov_galerkin_number_of_rom_dofs"] = left_basis["number_of_modes"]
    rom_settings.update(component.get("legacy_rom_settings", {}))

    hyper_reduction = v2.get("hyper_reduction", {})
    _CheckConventionalFileNames(hyper_reduction, _HROM_FILES)
    hrom_settings = {}
    if "format" in hyper_reduction:
        hrom_settings["hrom_format"] = hyper_reduction["format"]
    hrom_settings.update(hyper_reduction.get("training_settings", {}))

    deployment = v2.get("deployment", {})
    v1 = {
        "rom_manager" : deployment.get("rom_manager", _V1_DEFAULTS["rom_manager"]),
        "train_hrom" : deployment.get("train_hyper_reduction", _V1_DEFAULTS["train_hrom"]),
        "run_hrom" : deployment.get("use_hyper_reduction", _V1_DEFAULTS["run_hrom"]),
        "projection_strategy" : strategy,
        "assembling_strategy" : projection.get("assembling_strategy", _V1_DEFAULTS["assembling_strategy"]),
        "rom_format" : basis.get("format", _V1_DEFAULTS["rom_format"]),
        "rom_settings" : rom_settings,
        "hrom_settings" : hrom_settings,
        "nodal_modes" : basis.get("nodal_modes", {}),
        "elements_and_weights" : hyper_reduction.get("elements_and_weights", {})
    }
    if "nodal_modes" in left_basis:
        v1["petrov_galerkin_nodal_modes"] = left_basis["nodal_modes"]
    if "coupling" in v2:
        v1["coupled_solvers"] = v2["coupling"]["solvers"]
    v1.update(v2.get("legacy", {}))
    return v1


def NormalizeVersion1(v1_data):
    """Returns a copy of a version 1 dictionary with the defaults assumed by the readers filled in.
    Two version 1 dictionaries describe the same ROM if their normalized forms are equal."""
    v1 = copy.deepcopy(v1_data)
    for key, default in _V1_DEFAULTS.items():
        v1.setdefault(key, copy.deepcopy(default))
    v1["rom_settings"].setdefault("rom_bns_settings", {})
    return v1


def _CheckConventionalFileNames(block, conventional_names):
    for key, conventional_name in conventional_names.items():
        if key in block and block[key] != conventional_name:
            raise Exception(f"'{key}': '{block[key]}' cannot be written in the version 1 layout, which requires '{conventional_name}'.")
