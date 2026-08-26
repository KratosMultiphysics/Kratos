# Meshio++ Application

<p align="center">
  <a href="https://github.com/loumalouomega/meshioplusplus"><img alt="meshio++" src="https://raw.githubusercontent.com/loumalouomega/meshioplusplus/master/doc/logo/logo-with-text.svg" width="100%"></a>
</p>

|            **Application**            |                                                                                                    **Description**                                                                                                    |                                              **Status**                                              |                                 **Authors**                                 |
|:-------------------------------------:|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|:----------------------------------------------------------------------------------------------------:|:---------------------------------------------------------------------------:|
| `MeshioPlusPlusApplication`           | The *Meshio++ Application* wraps the [meshio++](https://github.com/loumalouomega/meshioplusplus) library, bringing multi-format mesh input/output (43 readable, 46 writable formats) and a full mesh/data operations layer into *Kratos Multiphysics* | <img src="https://img.shields.io/badge/Status-%F0%9F%94%A7Maintained-blue"  width="300px"> | [*Vicente Mataix Ferrándiz*](mailto:vicente.mataix-ferrandiz@siemens.com) |

The application includes tests to check the proper functioning of the application.

## 😎 Features:

- **Multi-format mesh input/output through `MeshioPlusPlusIO`**

    * *43 readable and 46 writable formats, including the Kratos native `mdpa` and the GiD postprocess format*

    * *Format resolved explicitly or inferred from the file extension*

    * *Sub model parts mapped to meshio++ named regions, nesting included*

    * *Registered entity names (`SmallDisplacementElement3D4N`, ...) preserved across a round trip, instead of degrading to the generic cell-type name*

    * *`Properties` material data carried in both directions, so a `mdpa` round trip keeps the values and not just the ids*

- **Transient output**

    * *XDMF temporal collections written in a single file, GiD multi-step `.post` series, file series for every other format*

    * *Append mode continues an existing series, so a restarted analysis does not destroy the previous run's output*

    * *Crash-safe flushing, so a killed run still leaves a readable `.xdmf`*

- **Mesh cleanup and transformation**

    * *`clean` — weld coincident nodes, drop degenerate and duplicate cells, remove orphan nodes*

    * *`transform` — translate, scale, rotate, with optional vector/tensor field rotation*

    * *`convert_cells` — linearize, simplexify or elevate (linear ⇄ quadratic element conversion)*

- **Subsetting and extraction**

    * *`split` — by cell type, connected component, integer tag or named region*

    * *`crop` — axis-aligned bounding box, half-space, or a `cell_data` predicate (`crop_predicate`, sharing `refine`'s comparison vocabulary and its evaluator)*

    * *`slice` and `isosurface` — planar cross-sections and level sets*

    * *`extract_surface` and `extract_skin` — boundary extraction*

    * *`merge` — concatenate meshes, optionally welding coincident nodes*

- **Mesh improvement and partitioning**

    * *`refine` — uniform subdivision, or selective/adaptive subdivision (explicit cell list, named region, or a cell-data predicate, with a choice of conforming closure); `record_hierarchy` attaches the persistent parent/child arrays a multigrid caller resolves across the sequence of meshes it keeps. `subdivide` is the unconditional one-level sibling that needs no closure because it leaves no hanging node, and `decimate` the reverse — quadric-error edge collapse*

    * *`smooth` — Laplacian, Taubin and ODT coordinate smoothing*

    * *`reorder` — reverse Cuthill–McKee bandwidth reduction, Morton and Hilbert space-filling curves*

    * *`partition` — SFC or KaHIP domain decomposition, with shared-node ghost layers*

    * *`agglomerate` — group cells into coarser ones, and `undo_green` (`MeshioPlusPlusMeshOperations.UndoGreen`) — take back the green closure a selective `refine` added, which keeps a mesh conforming but carries no refinement information and degrades element quality*

- **Remeshing**

    * *`remesh` — ACVD surface remeshing of a triangle mesh, with isotropic, quadric or anisotropic metric, curvature gradation and boundary preservation*

    * *`remesh_volume` — a closed surface in, a tetrahedral volume out: a lattice cut against the surface's signed distance, so it shares the whole lattice/SDF settings block with `voxelize` and `compute_sdf`*

    * *`optimize_volume` — improves a tetrahedral mesh without changing its cell budget (vertex relocation plus 2-3/3-2 flips), with a monotone worst-element quality guarantee; `decimate_volume` is the tetrahedral counterpart of `decimate`*

- **Diagnostics and data operations**

    * *`stats` and `quality` — bounding box, volumes, scaled Jacobian, aspect ratio, inverted-cell detection*

    * *`diff` — structured mesh comparison with tolerances*

    * *`gradient` — gradient, divergence and curl of a nodal field (Green-Gauss or least-squares), exact for a linear field; `hessian` is its second-derivative sibling*

    * *`estimate_error` — Zienkiewicz-Zhu recovery-based error indicator with absolute, fraction or Dörfler marking: the array it marks is exactly what selective `refine`'s `predicate_array` was built to consume, closing the adaptive loop*

    * *`data_integrate` — the cell-measure-weighted integral and mean of an array, over the whole mesh and over every named region*

    * *`data_calc` — evaluate an expression over one or more arrays into a new one*

    * *`data_condition` — clamp, normalize or standardize an array's values*

    * *`data_manage` — keep, drop or rename arrays; `data_info` reports their shape and dtype*

    * *`point_data_to_cell_data` and `cell_data_to_point_data` — averaging transfer between nodes and cells*

    * *`interpolate` (`MeshioPlusPlusMeshOperations.Interpolate`, or the `MeshioInterpolateModeler`) — sample one mesh's field data onto another's geometry (nearest/barycentric, with extrapolation and conflict handling), and `ConservativeInterpolate` — its conservative sibling, redistributing by intersected cell measure so the summed quantity is preserved: the right choice for an extensive field (a mass, a heat load), the wrong one for an intensive one (a temperature)*

- **Regular grids and signed distance**

    * *`Grid` — a regular hexahedron lattice from nothing: the one *generator* in the library*

    * *`voxelize` — a lattice around a mesh, keeping the whole bounding box, only the cells a surface passes through, or only those inside it*

    * *`compute_sdf` and `DistanceToSurface` — signed distance to a surface, on a generated grid (uniform or adaptive octree) or attached to an existing model part; `CheckSurfaceWatertight` reports what is wrong with a skin in numbers rather than a bare flag*

    * *Set `"output"` to a registered variable name (`DISTANCE`) and the result comes back as real Kratos data — meshio++ names it `sdf:distance`, which no `Variable` can be, so without this it would be computed but unreachable*

- **Provenance**

    * *Every file written records how it was produced — the Kratos model part it came from, the target format and encoding, and the chain of operations applied on the way. `MeshioPlusPlusIO.GetProvenance()` reads it back, distinguishing a block meshio++ actually wrote from a leading comment left by something else. The `"provenance"` setting turns it off or makes it mandatory*

- **Kratos-native `mdpa` fidelity**

    * *Gapped and non-sequential node ids are read (they used to be rejected outright), which is what a deck left by a sub model part extraction or an entity removal actually looks like*

    * *`write_mdpa_ids` preserves the model part's own entity ids on write instead of renumbering to `1..n`*

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
> The application requires **meshio++ ABI 11** — **v10.17.0 or newer**, which includes **v10.20.0**. The pin is on `MESHIOPLUSPLUS_ABI_VERSION` rather than the release version: that counter moves only when the installed headers stop being compatible with an already-compiled consumer, so a release that cannot affect this application needs no rebuild (v10.20.0 is a pure version-number bump over v10.19.0 and holds the ABI at 11, as v10.0.0 did over v9.30.0 at ABI 6). It is a single ABI rather than a window on purpose — supporting several meant one `MESHIOPLUSPLUS_VERSION_AT_LEAST` guard per feature, and the operations exposed here would have needed one each. Getting it wrong is not silent — the C++ variants' `SOVERSION` is the ABI version and every translation unit carries a link-time sentinel naming it, so a skew is a link error rather than memory corruption. See meshio++'s [`doc/abi.md`](https://github.com/loumalouomega/meshioplusplus/blob/master/doc/abi.md) for the criterion.

> [!NOTE]
> Upgrading from any older meshio++ requires relinking this application once: the C++ variants' `SOVERSION` is the ABI version, so the needed library is `libmeshioplusplus_core_kratos.so.11`. Install the new meshio++ wholesale rather than part-upgrading — its headers reference a link-time sentinel an older library does not define, which fails closed rather than silently mixing.

> [!NOTE]
> Unlike the ABI 3 → 6 run, ABI 6 → 11 is one this application genuinely feels. ABI 5 and 6 came from format side-channel structs (`MedInfo`, `OpenFoamInfo`) it never names; 7 through 11 are all types it now passes by value — `RefineOptions` gained `mRecordHierarchy` (7), `RemeshOptions` grew twice (8, 9), `SmoothMethod` gained an explicit `: std::uint8_t` underlying type (10) and `MeshMetadata` gained the provenance fields (11, `sizeof` 256 → 288). What it picks up in exchange is the whole remeshing family (`remesh`, `remesh_volume`, `optimize_volume`, `decimate_volume`), `estimate_error` closing the adaptive loop against selective `refine`, `hessian`, `subdivide`, `agglomerate`, `undo_green`, `data_integrate`, `conservative_interpolate`, the GiD postprocess format in both directions, and provenance.

> [!TIP]
> Writing GiD needs a **zlib-enabled** meshio++ build (`gidpost` deflates unconditionally, so `MESHIOPLUSPLUS_WITH_GIDPOST` auto-disables without it), and the `hdf5` flavour additionally needs HDF5. Reading needs *nothing* for the ascii flavour, zlib for binary and HDF5 for hdf5 — so `gid` is readable in strictly more configurations than it is writable. `MeshioPlusPlusIO.IsFormatAvailable(Format.GID)` answers for the write side.

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

**Read (43):** `abaqus` `ansys` `ansysinp` `avsucd` `cgns` `dex` `dolfin` `ensight` `exodus` `flac3d` `flux` `freefem` `gid` `gmsh` `h5m` `hmf` `ip` `mdpa` `med` `medit` `mff` `mfm` `mphtxt` `nastran` `netgen` `obj` `off` `openfoam` `permas` `ply` `stl` `su2` `tecplot` `tetgen` `triangle` `ugrid` `unv` `vti` `vtk` `vtp` `vtu` `wkt` `xdmf`

**Write (46):** the same set (every readable format is writable since meshio++ v9.20.0 added the OpenFOAM polyMesh writer), plus `gmsh22`, `svg` and `tikz` (write-only).

`cgns`, `h5m`, `hmf`, `med` and the XDMF-HDF data path require an HDF5-enabled meshio++ build; `exodus` requires netCDF; `gid` requires zlib to write (see the tip above). A format compiled out is still resolved by extension and reports *why* it is unavailable rather than "unknown format".

`gid` has four on-disk flavours selected with `"gid_mode"` (`auto`/`ascii`/`binary`/`hdf5`/`ascii_zipped`): `ascii` writes a `<stem>.post.msh`/`<stem>.post.res` sibling pair, `binary` one `<stem>.post.bin`, `hdf5` one `<stem>.post.h5`. `auto` never resolves to `ascii_zipped` — no extension can express "zipped".

## ⚠️ Limitations:

- **Serial only.** meshio++ has no MPI, no distributed reader or writer and no communicator. The intended distributed workflow is `partition` with ghost layers, feeding an MPI assembly.
- `remesh` and `decimate` operate on triangle surfaces; `remesh_volume`, `optimize_volume` and `decimate_volume` on tetrahedra. Neither family accepts the other's cells — they are separate operations rather than modes of one, and say so by name when handed the wrong kind.
- `Matrix`-valued variables are not written, and a `Matrix`-valued `Properties` entry is skipped with a warning — meshio++ has no representation for it.
- `Properties` entries that are not numeric are left to the materials file: a `Begin Table` curve, and text values such as a constitutive law name (instantiating one needs Kratos's own registry, which meshio++ deliberately does not link). Everything numeric — scalars, integers and vector-valued variables alike — round-trips.
- The C++ `mdpa` reader throws by name on constructs it cannot represent (`Geometries`, `Mesh <id>`, `Constraints`, ...); `ReadOptions::mLenient` downgrades those to a warning and a skip.
- **GiD transient output buffers.** meshio++'s `write_gid_series` *pulls* every step through a callback while this IO is *pushed* one step per `PrintOutput`, so the steps are held in memory and the whole series is written in `CloseOutput()` — nothing is on disk before that. Unlike the XDMF writer, which streams. For a long run set `"time_series" : "file_series"` instead, which writes one `<stem>_<label>.post.msh` per step.
- Transient output is held open by the IO. Call `CloseOutput()` before deleting the output of a run: a live writer finalizes on destruction and would recreate the file, and since the series is opened in append mode the next run would then continue a series believed deleted. `MeshioOutputProcess` does this in `ExecuteFinalize`.
- The operations layer only writes a resulting array back onto the output model part when its name matches a registered `Variable<T>` with the right component count: Kratos stores non-historical/historical data keyed by `Variable` objects, not arbitrary strings, so an operation's own invented array names (`attach_quality`'s `"quality:scaled_jacobian"`, for instance) are computed but cannot be retrieved from the output model part — point the operation's own naming setting (`data_calc`'s `"output"`, `data_manage`'s renames, ...) at an existing variable name to get the result back.

## 🗎 Documentation:

The meshio++ library documentation, including a page per format and per operation, is available at [loumalouomega.github.io/meshioplusplus](https://loumalouomega.github.io/meshioplusplus).
