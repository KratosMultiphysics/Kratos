---
title: Meshio Output Process
keywords: meshio meshioplusplus output process core xdmf vtu gmsh med transient
tags: [process meshio output process]
sidebar: kratos_core_processes
summary: This document details the Meshio Output Process, which writes simulation results in any of the ~35 mesh formats supported by the bundled meshio++ library, including transient XDMF time series and multi-format file series.
---

# Meshio Output Process

## Overview 📝

The **Meshio Output Process** writes simulation results in any of the mesh file formats supported by the **meshio++** library bundled with the Kratos core (`external_libraries/meshioplusplus`): `vtu`, `vtk`, `gmsh`, `med`, `xdmf`, `abaqus`, `medit`, `stl`, `obj`, `ply`, `su2`, `tecplot`, `unv` and many more.

Behind the scenes this Python process is a thin wrapper for the C++ `MeshioPlusPlusIO` class (`kratos/input_output/meshioplusplus_io.h`), which converts the Kratos `ModelPart` through the meshio++ Kratos bridge (`kratos_bridge.hpp`) in a single O(n) pass and dispatches to the requested format writer.

The full list of formats available in your build can be queried at runtime:

```python
import KratosMultiphysics
print(KratosMultiphysics.MeshioPlusPlusIO.GetSupportedWriteFormats())
print(KratosMultiphysics.MeshioPlusPlusIO.GetSupportedReadFormats())
```

> **Build dependencies**: the HDF5-backed formats (`med`, `cgns`, `h5m`, `hmf` and the `HDF` data path of `xdmf`) are only available when Kratos was configured with HDF5 (CMake option `KRATOS_MESHIOPLUSPLUS_HDF5`, `ON` by default when HDF5 is found); the `exodus` format requires netCDF (`KRATOS_MESHIOPLUSPLUS_NETCDF`). All other formats are self-contained and always available. Requesting a compiled-out format raises an error naming the missing dependency.

## Transient output ⏱️

The process is called once per output step, and the underlying IO **extends the current output instead of overwriting it**:

- **XDMF** (`.xdmf`/`.xmf`, with `"time_series": "automatic"`, the default): all steps go into a **single file**. The mesh (geometry + topology) is written once; every output step appends one grid to the XDMF *temporal collection*, holding the `<Time Value="...">` and the nodal results. On the first write the IO checks whether the file already holds a valid time series — if it does (e.g. after a restart), it is **extended** rather than overwritten. Heavy data is stored according to `"xdmf_data_format"`: inline `XML` text, sibling `.bin` files (`Binary`), or a companion `.h5` file (`HDF`, recommended, requires HDF5). The produced files are readable by ParaView and by the meshio `TimeSeriesReader`.
- **Every other format**: one file per output step is written as a *file series* `<output_name>_<label>.<ext>`, where the label is the `STEP` or the `TIME` (see `"output_control_type"`), like the VTK output does.
- With `"time_series": "single_file"` every call overwrites the same file (useful for writing a final mesh).

## Usage and parameters 📋

Add it to the `output_processes` of your `ProjectParameters.json`:

```json
{
    "output_processes" : {
        "meshio_output" : [{
            "python_module" : "meshio_output_process",
            "kratos_module" : "KratosMultiphysics",
            "Parameters"    : {
                "model_part_name"                    : "Structure",
                "output_name"                        : "results.xdmf",
                "output_path"                        : "meshio_output",
                "output_control_type"                : "step",
                "output_interval"                    : 1,
                "nodal_solution_step_data_variables" : ["DISPLACEMENT", "REACTION"]
            }
        }]
    }
}
```

Or with the registry-based form:

```json
{
    "output_processes" : {
        "meshio_output" : [{
            "name"       : "KratosMultiphysics.MeshioOutputProcess",
            "Parameters" : { "..." : "..." }
        }]
    }
}
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `model_part_name` | string | - | Name of the model part to write (**required**) |
| `output_name` | string | `""` | Output file name; its extension selects the format when `format` is `"auto"`. Empty: the model part name |
| `output_path` | string | `"meshio_output"` | Folder the files are written into |
| `save_output_files_in_folder` | bool | `true` | If `false`, `output_name` is used as-is without `output_path` |
| `format` | string | `"auto"` | Any meshio++ format name (`"vtu"`, `"gmsh"`, `"med"`, `"xdmf"`, ...) or `"auto"` (resolve from the extension) |
| `file_format` | string | `"default"` | `"ascii"` or `"binary"` for the formats that support the choice (see below); `"default"` keeps the meshio++ registry default |
| `time_series` | string | `"automatic"` | `"automatic"` (XDMF: in-file append; others: file series), `"file_series"`, `"single_file"` |
| `output_control_type` | string | `"step"` | `"step"` or `"time"`: drives when output happens **and** the file-series label / XDMF time value |
| `output_interval` | double | `1.0` | Output frequency (in steps or in time, according to `output_control_type`) |
| `output_precision` | int | `7` | Digits of inline floating point data (XDMF `"XML"` data format) |
| `label_precision` | int | `4` | Digits of the time label in file-series names |
| `custom_name_prefix` | string | `""` | Prepended to the output file name |
| `custom_name_postfix` | string | `""` | Appended to the output file stem (before the step/time label) |
| `entity_type` | string | `"automatic"` | `"automatic"` (elements **and** conditions - note the divergence from VtkOutput's either/or), `"element"` (elements only), `"condition"` (conditions only) |
| `output_sub_model_parts` | bool | `false` | Additionally writes each first-level sub model part as its own output (own file series or own XDMF time series), named `<stem>_<parent>_<smp>` |
| `write_deformed_configuration` | bool | `false` | `false`: initial node coordinates (`X0`); `true`: current coordinates (`X`) |
| `write_ids` | bool | `false` | Writes `KRATOS_NODE_ID` as point data and `KRATOS_ELEMENT_ID`/`KRATOS_CONDITION_ID` + `PROPERTIES_ID` as cell data |
| `xdmf_data_format` | string | `"auto"` | `"auto"` (`HDF` if available, else `Binary`), `"XML"`, `"Binary"`, `"HDF"` |
| `nodal_solution_step_data_variables` | array | `[]` | Historical nodal variables written as point data |
| `nodal_data_value_variables` | array | `[]` | Non-historical nodal variables written as point data |
| `nodal_flags` | array | `[]` | Nodal flags written as point data (`1`/`0`, `-1` when the flag is undefined on the entity) |
| `element_data_value_variables` | array | `[]` | Non-historical elemental variables written as cell data |
| `element_flags` | array | `[]` | Element flags written as cell data |
| `condition_data_value_variables` | array | `[]` | Non-historical condition variables written as cell data |
| `condition_flags` | array | `[]` | Condition flags written as cell data |
| `gauss_point_variables_extrapolated_to_nodes` | array | `[]` | Integration point results extrapolated to the nodes (via `IntegrationValuesExtrapolationToNodesProcess`) and written as non-historical point data |
| `gauss_point_variables_in_elements` | array | `[]` | Integration point results averaged over the gauss points and written as cell data (elements and conditions) |

**Supported variable types** (mirroring the VTK output): `double`, `int`, `bool`, `array_1d<double, 3/4/6/9>` and `Vector` (size taken from the first entity). Variables of other types (e.g. `Matrix`) are skipped with a warning. For the gauss point lists the supported set is `double`, `int`, `bool`, `array_1d<double, 3/6>` and `Vector`.

**`file_format` support matrix**: the ascii/binary choice is honored for `vtu`, `vtk`, `gmsh`, `stl`, `ply`, `ansys` and `flac3d` (the formats whose meshio++ writers expose a binary flag). For every other format the setting is ignored with a warning and the registry default is used. For `xdmf` use `xdmf_data_format` instead.

**Cell data and mixed meshes**: when both elements and conditions are written (`entity_type: "automatic"`), cell arrays cover all cells - element rows first, then condition rows - and the rows of the entity kind a variable does not apply to are zero-filled.

**Sub model part output**: with `output_sub_model_parts: true` each first-level sub model part is written to its own file(s), including one XDMF temporal-collection file per sub model part in transient XDMF mode. Sub model parts whose entities reference nodes the sub model part does not contain raise a descriptive error.

## Importing meshes: modeler and solver input 📥

The same IO also reads every supported format **as a model part** (highest-dimension cells become `Elements`, lower-dimension cells become `Conditions`, and integer tag arrays such as gmsh physical groups become `SubModelParts` named `<tag-key>_<value>`, e.g. `gmsh_physical_1`).

**Option 1 — the `MeshioInputModeler`** (recommended; runs before the solver imports the model part, so use `"input_type": "use_input_model_part"` in the solver):

```json
{
    "modelers" : [{
        "name"       : "KratosMultiphysics.MeshioInputModeler",
        "Parameters" : {
            "input_filename"  : "geometry.msh",
            "input_format"    : "auto",
            "model_part_name" : "Structure"
        }
    }]
}
```

**Option 2 — directly in the solver settings**: `model_part_import_settings.input_type` accepts, in addition to `"mdpa"`, any meshio++ format name as well as `"meshio"`/`"auto"` (format resolved from the extension of `input_filename`):

```json
{
    "solver_settings" : {
        "model_part_import_settings" : {
            "input_type"     : "gmsh",
            "input_filename" : "geometry.msh"
        }
    }
}
```

**Option 3 — pure Python**:

```python
import KratosMultiphysics as KM

model = KM.Model()
model_part = model.CreateModelPart("Main")
KM.MeshioPlusPlusIO("geometry.msh").ReadModelPart(model_part)
```

## Known limitations ⚠️

- Serial only: distributed (MPI) model parts are rejected with an explicit error.
- Written entities carry the generic Kratos names (`Element3D4N`, `SurfaceCondition3D3N`, ...) on read; assign the final element technology with e.g. the `ReplaceElementsAndConditionsProcess` or `ConnectivityPreserveModeler`.
- Sub model parts are **read** from tagged formats but not written into the mesh files themselves (the meshio++ staging drops the grouping on export); use `output_sub_model_parts` to obtain per-sub-model-part output files instead.
- Higher-order cell node ordering follows the meshio convention; for `wedge15`/`pyramid13` *element* creation on read there is no registered generic Kratos element, so reading meshes whose bulk cells are of those types fails.
- `Matrix` variables are not exported (skipped with a warning), matching the VTK output.
- The generic core elements/conditions return no integration point results, so the `gauss_point_*` settings only produce meaningful data with application entities (FE elements with constitutive laws).
