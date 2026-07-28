# Meshio++ Application

<p align="center">
  <a href="https://github.com/loumalouomega/meshioplusplus"><img alt="meshio++" src="https://raw.githubusercontent.com/loumalouomega/meshioplusplus/master/doc/logo/logo-with-text.svg" width="100%"></a>
</p>

|            **Application**            |                                                                                                    **Description**                                                                                                    |                                              **Status**                                              |                                 **Authors**                                 |
|:-------------------------------------:|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|:----------------------------------------------------------------------------------------------------:|:---------------------------------------------------------------------------:|
| `MeshioPlusPlusApplication`           | The *Meshio++ Application* wraps the [meshio++](https://github.com/loumalouomega/meshioplusplus) library, bringing multi-format mesh input/output (41 readable, 43 writable formats) and a full mesh/data operations layer into *Kratos Multiphysics* | <img src="https://img.shields.io/badge/Status-%F0%9F%94%A7Maintained-blue"  width="300px"> | [*Vicente Mataix Ferrándiz*](mailto:vicente.mataix-ferrandiz@siemens.com) |

The application includes tests to check the proper functioning of the application.

## 😎 Features:

- **Multi-format mesh input/output through `MeshioPlusPlusIO`**

    * *41 readable and 43 writable formats, including the Kratos native `mdpa`*

    * *Format resolved explicitly or inferred from the file extension*

    * *Sub model parts mapped to meshio++ named regions, nesting included*

    * *Registered entity names (`SmallDisplacementElement3D4N`, ...) preserved across a round trip, instead of degrading to the generic cell-type name*

    * *`Properties` material data carried in both directions, so a `mdpa` round trip keeps the values and not just the ids*

- **Transient output**

    * *XDMF temporal collections written in a single file, file series for every other format*

    * *Append mode continues an existing series, so a restarted analysis does not destroy the previous run's output*

    * *Crash-safe flushing, so a killed run still leaves a readable `.xdmf`*

- **Mesh cleanup and transformation**

    * *`clean` — weld coincident nodes, drop degenerate and duplicate cells, remove orphan nodes*

    * *`transform` — translate, scale, rotate, with optional vector/tensor field rotation*

    * *`convert_cells` — linearize, simplexify or elevate (linear ⇄ quadratic element conversion)*

- **Subsetting and extraction**

    * *`split` — by cell type, connected component, integer tag or named region*

    * *`crop` — axis-aligned bounding box or half-space*

    * *`slice` and `isosurface` — planar cross-sections and level sets*

    * *`extract_surface` and `extract_skin` — boundary extraction*

    * *`merge` — concatenate meshes, optionally welding coincident nodes*

- **Mesh improvement and partitioning**

    * *`refine` — uniform subdivision, and `decimate` — quadric-error edge collapse*

    * *`smooth` — Laplacian and Taubin coordinate smoothing*

    * *`reorder` — reverse Cuthill–McKee bandwidth reduction, Morton and Hilbert space-filling curves*

    * *`partition` — SFC or KaHIP domain decomposition, with shared-node ghost layers*

- **Diagnostics and data operations**

    * *`stats` and `quality` — bounding box, volumes, scaled Jacobian, aspect ratio, inverted-cell detection*

    * *`diff` — structured mesh comparison with tolerances*

    * *`data_calc` — evaluate an expression over one or more arrays into a new one*

    * *`data_condition` — clamp, normalize or standardize an array's values*

    * *`data_manage` — keep, drop or rename arrays; `data_info` reports their shape and dtype*

    * *`point_data_to_cell_data` and `cell_data_to_point_data` — averaging transfer between nodes and cells*

    * *`interpolate` (`MeshioPlusPlusMeshOperations.Interpolate`, or the `MeshioInterpolateModeler`) — sample one mesh's field data onto another's geometry (nearest/barycentric, with extrapolation and conflict handling)*

    Field data (nodal/elemental/conditional variables, flags and ids) reaches all of the above through the same `nodal_solution_step_data_variables` / `nodal_data_value_variables` / `nodal_flags` / `element_data_value_variables` / `element_flags` / `condition_data_value_variables` / `condition_flags` / `gauss_point_variables_in_elements` / `write_ids` settings `MeshioPlusPlusIO` uses — pass them to `operation_settings` and the selected variables are staged into the meshio++ mesh before the operation runs and, where the result carries a matching registered `Variable` name, written back onto the output model part afterwards.

## 🛠️ Building:

meshio++ is consumed as a normal external dependency; **nothing is vendored into Kratos**. Build and install it with the C++ API and the `KRATOS` mesh backend:

```bash
cmake -S <meshioplusplus> -B build               \
  -DMESHIOPLUSPLUS_INSTALL_CPP=ON                \
  -DMESHIOPLUSPLUS_INSTALL_CPP_BACKENDS="KRATOS" \
  -DMESHIOPLUSPLUS_MESH_BACKEND=KRATOS           \
  -DMESHIOPLUSPLUS_BUILD_PYTHON=OFF              \
  -DBUILD_SHARED_LIBS=ON
cmake --build build && cmake --install build --prefix <prefix>
```

then point the Kratos configure at it:

```bash
-Dmeshioplusplus_DIR=<prefix>/lib/cmake/meshioplusplus
```

> [!IMPORTANT]
> The application requires **meshio++ ABI 3**, which is **v9.2.0 or newer**. The pin is on `MESHIOPLUSPLUS_ABI_VERSION` rather than the release version: that counter moves only when the installed headers stop being compatible with an already-compiled consumer, so a release that cannot affect this application needs no rebuild. Getting it wrong is not silent — the C++ variants' `SOVERSION` is the ABI version and every translation unit carries a link-time sentinel naming it, so a skew is a link error rather than memory corruption. See meshio++'s [`doc/abi.md`](https://github.com/loumalouomega/meshioplusplus/blob/master/doc/abi.md) for the criterion.

> [!NOTE]
> Upgrading from a meshio++ older than v9.4.0 requires relinking this application once: the C++ variants' `SOVERSION` changed from a flat `0` to the ABI version, so the needed library is now `libmeshioplusplus_core_kratos.so.3` rather than `.so.0`. Install the new meshio++ wholesale rather than part-upgrading — v9.4.0 headers reference a sentinel a `≤9.3.0` library does not define, which fails closed at link time.

## 📖 Usage:

Reading a mesh in any supported format, through the modeler:

```json
{
    "modeler_name" : "KratosMultiphysics.MeshioPlusPlusApplication.modelers.meshio_input_modeler.MeshioInputModeler",
    "Parameters"   : {
        "model_part_name" : "Structure",
        "input_filename"  : "bracket.msh"
    }
}
```

or directly from the solver settings, where any supported format name works as an `input_type`:

```json
"model_import_settings" : {
    "input_type"     : "vtu",
    "input_filename" : "bracket.vtu"
}
```

Writing results, through the output process:

```json
{
    "python_module" : "meshio_output_process",
    "kratos_module" : "KratosMultiphysics.MeshioPlusPlusApplication",
    "Parameters"    : {
        "model_part_name" : "Structure",
        "output_name"     : "results.xdmf",
        "nodal_solution_step_data_variables" : ["DISPLACEMENT"]
    }
}
```

Query what the current build supports:

```python
import KratosMultiphysics.MeshioPlusPlusApplication as KratosMeshioPlusPlus

print(KratosMeshioPlusPlus.MeshioPlusPlusIO.GetSupportedReadFormats())
print(KratosMeshioPlusPlus.MeshioPlusPlusIO.GetSupportedWriteFormats())
```

## 📁 Supported formats:

**Read (41):** `abaqus` `ansys` `ansysinp` `avsucd` `cgns` `dex` `dolfin` `ensight` `exodus` `flac3d` `flux` `freefem` `gmsh` `h5m` `hmf` `ip` `mdpa` `med` `medit` `mff` `mfm` `mphtxt` `nastran` `netgen` `obj` `off` `openfoam` `permas` `ply` `stl` `su2` `tecplot` `tetgen` `triangle` `ugrid` `unv` `vtk` `vtp` `vtu` `wkt` `xdmf`

**Write (43):** the same set except `openfoam` (read-only), plus `gmsh22`, `svg` and `tikz` (write-only).

`cgns`, `h5m`, `hmf`, `med` and the XDMF-HDF data path require an HDF5-enabled meshio++ build; `exodus` requires netCDF. A format compiled out is still resolved by extension and reports *why* it is unavailable rather than "unknown format".

## ⚠️ Limitations:

- **Serial only.** meshio++ has no MPI, no distributed reader or writer and no communicator. The intended distributed workflow is `partition` with ghost layers, feeding an MPI assembly.
- `Matrix`-valued variables are not written, and a `Matrix`-valued `Properties` entry is skipped with a warning — meshio++ has no representation for it.
- `Properties` entries that are not numeric are left to the materials file: a `Begin Table` curve, and text values such as a constitutive law name (instantiating one needs Kratos's own registry, which meshio++ deliberately does not link). Everything numeric — scalars, integers and vector-valued variables alike — round-trips.
- The C++ `mdpa` reader throws by name on constructs it cannot represent (`Geometries`, `Mesh <id>`, `Constraints`, ...); `ReadOptions::mLenient` downgrades those to a warning and a skip.
- Transient output is held open by the IO. Call `CloseOutput()` before deleting the output of a run: a live writer finalizes on destruction and would recreate the file, and since the series is opened in append mode the next run would then continue a series believed deleted. `MeshioOutputProcess` does this in `ExecuteFinalize`.
- The operations layer only writes a resulting array back onto the output model part when its name matches a registered `Variable<T>` with the right component count: Kratos stores non-historical/historical data keyed by `Variable` objects, not arbitrary strings, so an operation's own invented array names (`attach_quality`'s `"quality:scaled_jacobian"`, for instance) are computed but cannot be retrieved from the output model part — point the operation's own naming setting (`data_calc`'s `"output"`, `data_manage`'s renames, ...) at an existing variable name to get the result back.

## 🗎 Documentation:

The meshio++ library documentation, including a page per format and per operation, is available at [loumalouomega.github.io/meshioplusplus](https://loumalouomega.github.io/meshioplusplus).
