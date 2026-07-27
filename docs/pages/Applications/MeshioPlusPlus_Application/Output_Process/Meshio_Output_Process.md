---
title: Meshio Output Process
keywords: meshio meshioplusplus output process xdmf vtu gmsh med transient
tags: [Meshio_Output_Process.md]
sidebar: meshioplusplus_application
summary: Writes simulation results in any of the mesh formats meshio++ supports, including transient XDMF time series and multi-format file series.
---

## Overview

The **Meshio Output Process** (`meshio_output_process.MeshioOutputProcess`) writes simulation results in any of the formats listed in the [Overview](../General/Overview.html#supported-formats).

It is a thin Python wrapper around the C++ `MeshioPlusPlusIO` class (`custom_io/meshioplusplus_io.h`), which converts the Kratos `ModelPart` through the meshio++ Kratos bridge in a single O(n) pass and dispatches to the requested format writer.

## Transient output

The process is called once per output step, and the underlying IO **extends the current output instead of overwriting it**:

- **XDMF** (`.xdmf`/`.xmf`, `"time_series": "automatic"`, the default): all steps go into a **single file**. The static mesh (geometry + topology) is written once; every output step appends one grid to the XDMF *temporal collection*, holding `<Time Value="...">` and the nodal/cell results. On the first write of a run the IO checks whether the file already holds a valid time series (e.g. after a restart) and, if so, **continues** it instead of overwriting it. Heavy data is stored according to `"xdmf_data_format"`: inline `XML` text, sibling `.bin` files (`Binary`), or a companion `.h5` file (`HDF`, recommended, requires HDF5). The light data (`.xdmf`) is rewritten after every step by default (`"xdmf_auto_flush"`), so a run that is killed still leaves a file readable in ParaView.
- **Every other format**: one file per output step is written as a file series `<output_name>_<label>.<ext>`, where the label is the `STEP` or the `TIME` (see `"output_control_type"`), like the VTK output does.
- With `"time_series": "single_file"` every call overwrites the same file (useful for writing a single, final mesh).

## Discovering what was written

`MeshioPlusPlusIO.GetNumberOfTimeSteps()` / `.GetTimeValues()` / `.GetTimeStepIndex(time)` answer "how many steps, and at what times" without reading any mesh data, for either transient shape:

- A single transient file (XDMF, or any format meshio++ itself reports time values for) is summarized through its own header.
- A **file series** is detected from disk: the directory is scanned for every `<prefix><stem><postfix>_<label><ext>` entry matching the naming `ComposeOutputPath` uses (same `"custom_name_prefix"`/`"custom_name_postfix"` settings as the writer), each `<label>` is parsed as a number, and the results are returned sorted ascending — regardless of write order or filesystem enumeration order.

```python
meshio_io = KratosMultiphysics.MeshioPlusPlusApplication.MeshioPlusPlusIO("results.vtu")
meshio_io.GetNumberOfTimeSteps()   # 3, discovered from results_1.vtu, results_2.vtu, results_5.vtu
meshio_io.GetTimeValues()          # [1.0, 2.0, 5.0]
meshio_io.GetTimeStepIndex(2.0)    # 1
```

A file series only carries as much time precision as `"label_precision"` kept when it was written, so a `GetTimeStepIndex` lookup finer than that will not match. A file that does not parse as numeric (an unrelated file sharing the extension, say) is silently skipped rather than failing the scan.

## Closing the output

The IO holds the time series open for its whole lifetime, and `ExecuteFinalize` closes it by calling `MeshioPlusPlusIO.CloseOutput()`. That finalizes the `.xdmf` light data and releases the file, so the series ends when the analysis ends rather than whenever the IO happens to be garbage collected.

> Call `CloseOutput()` before deleting the output of a transient run. A live writer finalizes on destruction, which **recreates** a file deleted underneath it — and because the series is opened in append mode, the *next* run then continues a series you believed deleted, so the mistake surfaces one run later as a wrong step count. A further `WriteModelPart` after `CloseOutput()` simply reopens the series and appends to what is on disk.

## Usage

```json
{
    "output_processes" : {
        "meshio_output" : [{
            "python_module" : "meshio_output_process",
            "kratos_module" : "KratosMultiphysics.MeshioPlusPlusApplication",
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

Or through the [Python registry](https://github.com/KratosMultiphysics/Kratos/wiki/Python-Orchestration-Classes):

```json
{
    "output_processes" : {
        "meshio_output" : [{
            "Parameters" : {
                "help"            : "meshio_output_process.MeshioOutputProcess",
                "model_part_name" : "Structure",
                "output_name"     : "results.xdmf"
            }
        }]
    }
}
```

## Parameters

| Setting | Default | Description |
|---|---|---|
| `model_part_name` | *(required)* | Model part to write. |
| `output_name` | `""` | Output file (name or path); the model part name when empty. |
| `output_path` | `"meshio_output"` | Folder the output is written into. |
| `save_output_files_in_folder` | `true` | Create `output_path` and write into it. |
| `format` | `"auto"` | Explicit format name, or `"auto"` to resolve it from `output_name`'s extension. |
| `file_format` | `"default"` | `"ascii"`/`"binary"` override, for the formats that support it. |
| `skin` | `true` | For surface-only formats (`stl`, `ply`): extract and write the boundary skin of a volume mesh. |
| `time_series` | `"automatic"` | `"automatic"` (see above), `"file_series"`, or `"single_file"`. |
| `output_control_type` | `"step"` | `"step"` or `"time"`, matching the VTK output process. |
| `output_interval` | `1.0` | How often to write, in the units of `output_control_type`. |
| `output_precision` | `7` | Digits for inline `XML` floating point data. |
| `label_precision` | `4` | Digits used for the file-series step/time label. |
| `custom_name_prefix` / `custom_name_postfix` | `""` | Extra text around the output file name. |
| `entity_type` | `"automatic"` | `"elements"`, `"conditions"`, or `"automatic"` (both). |
| `output_sub_model_parts` | `false` | Write one file per sub model part, keyed by its name. |
| `write_deformed_configuration` | `false` | Write current (`X()`) coordinates instead of the initial (`X0()`) ones. |
| `write_ids` | `false` | Attach node/element/condition ids as extra data arrays. |
| `xdmf_data_format` | `"auto"` | `"HDF"` (needs HDF5), `"XML"`, or `"Binary"`; `"auto"` picks `HDF` when available. |
| `xdmf_auto_flush` | `true` | Rewrite the `.xdmf` light data after every step, so a killed run stays readable. |
| `xdmf_gzip_level` | `-1` | gzip level for `HDF` datasets; `-1` means uncompressed. |
| `nodal_solution_step_data_variables` | `[]` | Historical nodal variables to write. |
| `nodal_data_value_variables` | `[]` | Non-historical nodal variables to write. |
| `nodal_flags` | `[]` | Nodal flags to write (1/0/-1: set/unset/undefined). |
| `element_data_value_variables` / `condition_data_value_variables` | `[]` | Non-historical element/condition variables to write. |
| `element_flags` / `condition_flags` | `[]` | Element/condition flags to write. |
| `gauss_point_variables_extrapolated_to_nodes` | `[]` | Gauss point variables extrapolated to nodes and written as nodal data. |
| `gauss_point_variables_in_elements` | `[]` | Gauss point variables averaged and written as cell data. |

See [`GetDefaultParameters`](../General/Overview.html) on `MeshioOutputProcess` for the authoritative, up-to-date list.
