---
title: Overview
keywords: meshio meshioplusplus mesh io formats xdmf vtu gmsh med mdpa operations
tags: [Overview.md]
sidebar: meshioplusplus_application
summary: What the Meshio++ Application is, how to build it, and where its features live.
---

## What is the Meshio++ Application?

The **Meshio++ Application** wraps the [meshio++](https://github.com/loumalouomega/meshioplusplus) library, bringing multi-format mesh input/output and a full mesh/data operations layer into *Kratos Multiphysics*.

meshio++ is consumed as a normal, **independently built external dependency** — nothing is vendored into Kratos. The only meshio++-shaped code that lives in this application is the bridge mapping Kratos entities to meshio++ ones (`custom_utilities/meshioplusplus_conversion_utilities.h`); everything else — the format readers/writers and the mesh operations — comes from the installed meshio++ C++ API.

## Features

- **Multi-format mesh input/output** through `MeshioPlusPlusIO` — see [Output Process](../Output_Process/Meshio_Output_Process.html) and [Modelers](../Modelers/Meshio_Input_Modeler.html).
  - 41 readable and 43 writable formats, including the Kratos native `mdpa`.
  - The format is resolved explicitly or inferred from the file extension.
  - Sub model parts map to meshio++ named regions, nesting included.
  - Registered entity names (`SmallDisplacementElement3D4N`, ...) are preserved across a round trip instead of degrading to the generic cell-type name.
  - `Properties` material data is carried in both directions, so a `mdpa` round trip keeps the values and not just the ids.
- **Transient output**: XDMF temporal collections in a single file, file series for every other format. Append mode continues an existing series (a restart does not destroy the previous run's output), and crash-safe flushing means a killed run still leaves a readable `.xdmf`. `GetNumberOfTimeSteps`/`GetTimeValues`/`GetTimeStepIndex` introspect either kind — a single transient file or a file series discovered from the output naming convention — without reading the mesh data.
- **Mesh and data operations** through `MeshioPlusPlusMeshOperations` — see [Utilities](../Utilities/Mesh_Operations.html): cleanup, transform, cell-type conversion, subsetting/extraction (split, crop, slice, isosurface, skin/surface extraction, merge), mesh improvement (refine, decimate, smooth, reorder), partitioning with ghost layers, diagnostics (stats, quality, diff), field-data manipulation (data_calc, data_condition, data_manage, data_info, point_data_to_cell_data, cell_data_to_point_data) and cross-mesh sampling (`Interpolate`, reachable through the dedicated `MeshioInterpolateModeler`).

## Building

Build and install meshio++ with the C++ API and the `KRATOS` mesh backend:

```bash
cmake -S <meshioplusplus> -B build \
  -DMESHIOPLUSPLUS_INSTALL_CPP=ON \
  -DMESHIOPLUSPLUS_INSTALL_CPP_BACKENDS="KRATOS" \
  -DMESHIOPLUSPLUS_MESH_BACKEND=KRATOS \
  -DMESHIOPLUSPLUS_BUILD_PYTHON=OFF \
  -DBUILD_SHARED_LIBS=ON
cmake --build build && cmake --install build --prefix <prefix>
```

then point the Kratos configure at it, and add the application:

```bash
-Dmeshioplusplus_DIR=<prefix>/lib/cmake/meshioplusplus
```

```python
add_app ${KRATOS_APP_DIR}/MeshioPlusPlusApplication
```

> **Version pin**: the application requires **meshio++ ABI 3 or ABI 4** — **v9.2.0 or newer** (ABI 3 covers v9.2.0-v9.4.1, ABI 4 covers v9.5.0 onwards). It pins `MESHIOPLUSPLUS_ABI_VERSION` rather than the release version — that counter moves only when the installed headers stop being compatible with an already-compiled consumer, so a release that cannot affect the application needs no rebuild:
>
> ```cmake
> find_package(meshioplusplus CONFIG REQUIRED COMPONENTS CXX)
> if(MESHIOPLUSPLUS_ABI_VERSION LESS 3 OR MESHIOPLUSPLUS_ABI_VERSION GREATER 4)
>     message(FATAL_ERROR "...")
> endif()
> ```
>
> Getting it wrong is not silent: the C++ variants' `SOVERSION` is the ABI version, and every consumer translation unit carries a link-time sentinel naming it, so a header/library skew is a link error rather than memory corruption. See meshio++'s [`doc/abi.md`](https://github.com/loumalouomega/meshioplusplus/blob/master/doc/abi.md) for the criterion that decides when the ABI version moves.

> **Upgrading from a meshio++ older than v9.4.0** requires relinking the application once: the C++ variants' `SOVERSION` changed from a flat `0` to the ABI version, so the needed library is now `libmeshioplusplus_core_kratos.so.3` rather than `.so.0`. Install the new meshio++ wholesale rather than part-upgrading — v9.4.0 headers reference a sentinel a `≤9.3.0` library does not define, which fails closed at link time.

> **Upgrading across the ABI 3 → ABI 4 boundary** (anything before v9.5.0 to v9.5.0+) likewise needs one relink, to `libmeshioplusplus_core_kratos.so.4`. ABI 4 moved because `RefineOptions` gained new data members for selective/adaptive refinement (a layout change); the application compiles around it with `MESHIOPLUSPLUS_VERSION_AT_LEAST(9, 5, 0)` rather than raising its ABI floor, so it keeps building against ABI 3 too — just without the new `refine` selectors, which are rejected with a clear error rather than silently ignored if set against an older build. See [Utilities](../Utilities/Mesh_Operations.html) for the settings.
>
> **meshio++ v9.6.0's MED improvements** (named regions via `FAS`/`GRO` groups, `med:num` global numbering, stricter MED major-version validation) needed no code changes here — `MeshioPlusPlusIO` treats MED like every other format through the generic region/format-table machinery.

## Supported formats

**Read (41):** `abaqus` `ansys` `ansysinp` `avsucd` `cgns` `dex` `dolfin` `ensight` `exodus` `flac3d` `flux` `freefem` `gmsh` `h5m` `hmf` `ip` `mdpa` `med` `medit` `mff` `mfm` `mphtxt` `nastran` `netgen` `obj` `off` `openfoam` `permas` `ply` `stl` `su2` `tecplot` `tetgen` `triangle` `ugrid` `unv` `vtk` `vtp` `vtu` `wkt` `xdmf`

**Write (43):** the same set except `openfoam` (read-only), plus `gmsh22`, `svg` and `tikz` (write-only).

`cgns`, `h5m`, `hmf`, `med` and the XDMF-HDF data path require an HDF5-enabled meshio++ build; `exodus` requires netCDF. A format compiled out is still resolved by extension and reports *why* it is unavailable rather than "unknown format". Query what your build supports at runtime:

```python
import KratosMultiphysics.MeshioPlusPlusApplication as KratosMeshioPlusPlus

print(KratosMeshioPlusPlus.MeshioPlusPlusIO.GetSupportedReadFormats())
print(KratosMeshioPlusPlus.MeshioPlusPlusIO.GetSupportedWriteFormats())
```

## Limitations

- **Serial only.** meshio++ has no MPI, no distributed reader or writer and no communicator. The intended distributed workflow is `partition` with ghost layers, feeding an MPI assembly.
- `Matrix`-valued variables are not written, and a `Matrix`-valued `Properties` entry is skipped with a warning — meshio++ has no representation for it.
- `Properties` entries that are not numeric are left to the materials file: a `Begin Table` curve, and text values such as a constitutive law name (instantiating one needs Kratos's own registry, which meshio++ deliberately does not link). Everything numeric — scalars, integers and vector-valued variables alike — round-trips.
- The `mdpa` reader throws by name on constructs it cannot represent (`Geometries`, `Mesh <id>`, `Constraints`, ...).
- Transient output is held open by the IO. Call `CloseOutput()` before deleting the output of a run: a live writer finalizes on destruction and would recreate the file, and since the series is opened in append mode the next run would then continue a series believed deleted. `MeshioOutputProcess` does this in `ExecuteFinalize`.

## Documentation

The meshio++ library documentation, including a page per format and per operation, is available at [loumalouomega.github.io/meshioplusplus](https://loumalouomega.github.io/meshioplusplus).
