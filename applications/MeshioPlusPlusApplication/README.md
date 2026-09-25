# Meshio++ Application

<p align="center">
  <a href="https://github.com/loumalouomega/meshioplusplus"><img alt="meshio++" src="https://raw.githubusercontent.com/loumalouomega/meshioplusplus/master/doc/logo/logo-with-text.svg" width="100%"></a>
</p>

|            **Application**            |                                                                                                    **Description**                                                                                                    |                                              **Status**                                              |                                 **Authors**                                 |
|:-------------------------------------:|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|:----------------------------------------------------------------------------------------------------:|:---------------------------------------------------------------------------:|
| `MeshioPlusPlusApplication`           | The *Meshio++ Application* wraps the [meshio++](https://github.com/loumalouomega/meshioplusplus) library, bringing multi-format mesh input/output (78 readable, 66 writable formats) and a full mesh/data operations layer into *Kratos Multiphysics* | <img src="https://img.shields.io/badge/Status-%F0%9F%94%A7Maintained-blue"  width="300px"> | [*Vicente Mataix Ferrándiz*](mailto:vicente.mataix-ferrandiz@siemens.com) |

The application includes tests to check the proper functioning of the application.

## 😎 Features:

- **Multi-format mesh input/output through `MeshioPlusPlusIO`**

    * *78 readable and 66 writable formats, including the Kratos native `mdpa`, the GiD postprocess format, the ParaView `pvtu`/`pvtp`/`pvd` indexes, VTKHDF, the FEM interchange decks (Elmer, FEBio, Femap, libMesh, MFEM, Patran, Z88), LS-DYNA and Code_Aster decks, the solver results files of Abaqus, Ansys, Nastran, LS-DYNA, Marc, Radioss, FEBio and DOLFINx, CalculiX and MSC Nastran HDF5 results, point clouds and glTF*

    * *Format resolved explicitly or inferred from the file extension*

    * *Sub model parts mapped to meshio++ named regions, nesting included*

    * *Registered entity names (`SmallDisplacementElement3D4N`, ...) preserved across a round trip, instead of degrading to the generic cell-type name*

    * *`Properties` material data carried in both directions, so a `mdpa` round trip keeps the values and not just the ids*

    * *`"time_step"` (negative counts from the end) selects one step of a multi-step file for the formats meshio++ reads selectively — mdpa, med, exodus, gmsh, tecplot, ensight, cgns, openfoam, frd, unv, nastran_h5, vtkhdf and pvd — through a cheap native metadata reader where one exists; any other format ignores it and is always read whole*

    * *`"lenient"` downgrades an mdpa/med/vtkhdf construct the reader cannot represent from an error to a warning and a skip, instead of throwing*

    * *A partitioned file (`vtkhdf`, `pvtu`, `pvtp`, `pvd`) is read whole, one sub model part per piece; `"select_piece"` with `"piece"` (negative counts from the end) reads one piece instead, and `"ghosts" : "drop"` removes the halo cells a `pvtu`/`pvtp`/`pvd` carries as `vtkGhostType`*

    * *Solver results files (`abaqus_fil`, `ansys_rst`, `nastran_op2`, `lsdyna_d3plot`, `marc_t19`, `xplt`, `vtx`, ...) read one step per `"time_step"`, `GetTimeValues()` lists the steps, and `"read_field_data" : true` carries their arrays onto the registered variables of the same name*

    * *A format meshio++ cannot place by extension or file name — an Elmer mesh directory, an extensionless deck — is identified by its content*

    * *Companion field files: `"mfem_grid_functions"` and `"patran_result_files"` (`{"name", "path"}` lists) read the fields an MFEM `.gf` or Patran `.nod`/`.dis`/`.els` file holds; `"mfem_grid_functions_write"` writes the nodal data next to an MFEM mesh as `<stem>.<name>.gf`, ready for an MFEM solver; `"z88_results"` and `"z88_stubs"` control a Z88 deck's results and input stubs*

    * *Format-specific writers: `"gltf_settings"` (`container`, `up_axis`, `split_angle`, `color_by` with `cmap`/`range`/`component`, ...) for glTF, `"openfoam_label_bits"`/`"openfoam_scalar_bits"` with `"file_format" : "binary"` for a binary polyMesh, `"pcd_compressed"`/`"pcd_float64_points"` for PCL's `binary_compressed`, `"vtkhdf_gzip_level"`; every other format meshio++'s ascii/binary switch covers honours `"file_format"`*

    * *`"openfoam_region"` selects one region of an OpenFOAM case with no single `constant/polyMesh`, but `constant/<region>/polyMesh` per region; FLAC3D `ZGROUP`/`FGROUP` cell groups round-trip as named regions instead of always deferring to the Python fallback; a gmsh 4.1 export now allocates a physical tag for a named `Cell` region that has none of its own, instead of dropping it*

- **Transient output**

    * *XDMF temporal collections and VTKHDF `Steps` written in a single file, ParaView `.pvd` collections (one `.vtu` per step, index rewritten after every step), GiD multi-step `.post` series, file series for every other format*

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

    * *`subdivide` and `agglomerate` — polyhedral refinement and coarsening. Both emit polyhedral cells, which no Kratos `Element` can hold, so `"simplexify_result"` (on by default) decomposes the result into tetrahedra before it reaches the model part*

    * *`repair` — surface topology fix: rewinds triangles so neighbours agree, fans-fills boundary loops up to `"max_hole_edges"`, and duplicates edge-disconnected (bowtie) vertices. Reports both the input's and the output's defect counts, so what was fixed and what remains are both visible*

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

    * *`tensor_invariants` — von Mises, principal values, hydrostatic and deviatoric parts of a six- (`xx yy zz xy yz zx`) or nine-component tensor array; `"rename"` lands a result (`CAUCHY_STRESS_VECTOR_mises`) on a registered variable*

    * *`compute_normals` — unit point (and optionally cell) normals of a surface, angle/area/uniform weighted, with `"split"` duplicating the points where the surface creases by more than `"split_angle"`; `"output" : "NORMAL"` brings them back as Kratos data*

    * *`curvature` — per-vertex mean and Gaussian curvature of a surface (angle defect / cotangent Laplace-Beltrami), with `"total_angle_defect"` reporting the Gauss-Bonnet invariant — `2·π·χ` for a closed surface, whatever the tessellation — as a checkable oracle*

- **Geometry fitting and deformation**

    * *`Shrinkwrap` (`MeshioPlusPlusMeshOperations.Shrinkwrap`) — projects one mesh's points onto a target surface, optionally offset along the hit feature's pseudonormal; one projection, not an iteration, so a fit rather than a smoothing. Needs two independent meshes, like `Interpolate`*

    * *`sobolev_deform` — Sobolev-filters a raw nodal displacement field before moving the mesh's points by it, damping high-frequency noise a caller already computed (a `"length_scale"` of `0.0`, the default, applies it directly with no solve at all)*

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

Three formats need an optional dependency, all off by default: `-DMESHIOPLUSPLUS_WITH_ADIOS2=ON` for DOLFINx `vtx` (`.bp`), `-DMESHIOPLUSPLUS_WITH_TECIO=ON -DTECIO_ROOT=<dir>` for Tecplot `szplt`, and `-DMESHIOPLUSPLUS_WITH_BZIP2=ON` for bzip2-compressed libMesh (`.xda.bz2`). The configure line of this application reports which ones the install carries.

then point the Kratos configure at it:

```bash
-Dmeshioplusplus_DIR=<prefix>/lib/cmake/meshioplusplus
```

> [!IMPORTANT]
> The application requires **meshio++ ABI 16**, release **v16.13.0 or newer** — the release floor on top of the ABI pin is there because the newer formats and the `read_mfem`/`write_mfem`/`read_patran`/`read_z88`/`write_z88` overloads this application calls arrived in later ABI-16 releases (v16.2.0–v16.13.0), all as additions. The pin is on `MESHIOPLUSPLUS_ABI_VERSION` rather than the release version: that counter moves only when the installed headers stop being compatible with an already-compiled consumer, so a release that cannot affect this application needs no rebuild (v12.0.0 is a pure version-number bump over v11.6.0 and holds the ABI at 13, as v10.20.0 did over v10.19.0 at ABI 11). It is a single ABI rather than a window on purpose — supporting several meant one `MESHIOPLUSPLUS_VERSION_AT_LEAST` guard per feature, and the operations exposed here would have needed one each. Getting it wrong is not silent — the C++ variants' `SOVERSION` is the ABI version and every translation unit carries a link-time sentinel naming it, so a skew is a link error rather than memory corruption. See meshio++'s [`doc/abi.md`](https://github.com/loumalouomega/meshioplusplus/blob/master/doc/abi.md) for the criterion.

> [!NOTE]
> Upgrading from any older meshio++ requires relinking this application once: the C++ variants' `SOVERSION` is the ABI version, so the needed library is `libmeshioplusplus_core_kratos.so.16`. Install the new meshio++ wholesale rather than part-upgrading — its headers reference a link-time sentinel an older library does not define, which fails closed rather than silently mixing.

> [!NOTE]
> Unlike the ABI 3 → 6 run, ABI 6 → 11 is one this application genuinely feels. ABI 5 and 6 came from format side-channel structs (`MedInfo`, `OpenFoamInfo`) it never names; 7 through 11 are all types it now passes by value — `RefineOptions` gained `mRecordHierarchy` (7), `RemeshOptions` grew twice (8, 9), `SmoothMethod` gained an explicit `: std::uint8_t` underlying type (10) and `MeshMetadata` gained the provenance fields (11, `sizeof` 256 → 288). What it picks up in exchange is the whole remeshing family (`remesh`, `remesh_volume`, `optimize_volume`, `decimate_volume`), `estimate_error` closing the adaptive loop against selective `refine`, `hessian`, `subdivide`, `agglomerate`, `undo_green`, `data_integrate`, `conservative_interpolate`, the GiD postprocess format in both directions, and provenance.

> [!NOTE]
> ABI 11 → 13 touches this application far less. 12 (v10.35.0) is a Tier B ODR bump — no struct this application names changed size, but `NDArray::Size()`'s inline body used to report 0 for a rank-0 (scalar) array, so a 0-d `field_data` value was silently destroyed by every operation that clones a mesh; a plain relink against the new library fixes it, with no source change here. 13 (v11.4.0) grows `OpenFoamInfo` (96 → 128 bytes) with a new `mRegion` field for multi-region case selection — a struct this application does not name at this integration level, so the layout change is inert. What it picks up in exchange: `compute_curvature`, `repair`, `shrinkwrap` and `sobolev_deform`; the `vts`/`vtr`/`vtm` VTK XML structured/multiblock formats; `time_step`-selected reads for MED, CGNS, Tecplot, Gmsh, EnSight and OpenFOAM (each with a cheap native metadata reader); FLAC3D `ZGROUP`/`FGROUP` cell groups read and written natively as named regions instead of always deferring to the Python fallback; and Gmsh 4.1 export allocating a physical tag for an untagged named region instead of dropping it.

> [!NOTE]
> ABI 13 → 16 is a rebuild, not a relink. 14 (v14.0.0) grows `ReadOptions` 56 → 72 bytes with `mPiece`/`mPieceSet`, and 15 (v15.0.0) puts `mGhosts` into its tail padding (the size holds, but a consumer built against older headers leaves that byte indeterminate) — this application builds a `ReadOptions` for every read, so both reach it. 16 (v16.0.0) appends `CellType::Triangle7` before `CellType::Custom`, moving `Custom` 76 → 77 and growing the inline cell-type tables. What it picks up in exchange: the `vtkhdf`, `pvtu`/`pvtp`/`pvd`, `pcd`/`xyz`, `lsdyna`, `code_aster`, `mphbin` formats, `frd` and `nastran_h5` (read-only) and `gltf` (write-only); piece selection and ghost dropping; single-file VTKHDF and PVD time series; binary OpenFOAM; `compute_normals` and `tensor_invariants`. Upstream behaviour changes that reach this application's I/O: Nastran writes the real `CTETRA`/`CPYRA`/`CPENTA`/`CHEXA` keywords for quadratic solids and reads any bulk-data deck; COMSOL uses COMSOL's own quadratic node order and writes a version 4 Mesh; I-DEAS UNV swaps its group entity types and the hexahedron20/wedge15 mid-edge order; multi-zone Tecplot reads every zone; `.vtu`/`.vtp` carry `field_data`.

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

**Read (78):** `abaqus` `abaqus_fil` `ansys` `ansys_rst` `ansys_rst_cyclic` `ansysinp` `avsucd` `cgns` `code_aster` `dex` `dolfin` `elmer` `ensight` `exodus` `febio` `femap` `flac3d` `flux` `frd` `freefem` `gid` `gmsh` `h5m` `hmf` `ip` `libmesh` `lsdyna` `lsdyna_binout` `lsdyna_d3plot` `marc` `marc_t19` `mdpa` `med` `medit` `mfem` `mff` `mfm` `mphbin` `mphtxt` `nastran` `nastran_h5` `nastran_op2` `netgen` `obj` `off` `openfoam` `patran` `pcd` `permas` `ply` `pvd` `pvtp` `pvtu` `radioss` `radioss_anim` `radioss_th` `stl` `su2` `szplt` `tecplot` `tetgen` `triangle` `ugrid` `unv` `vti` `vtk` `vtkhdf` `vtm` `vtp` `vtr` `vts` `vtu` `vtx` `wkt` `xdmf` `xplt` `xyz` `z88`

**Write (66):** the same set minus the read-only ones — the solver results files `abaqus_fil`, `ansys_rst`, `ansys_rst_cyclic`, `frd`, `lsdyna_binout`, `lsdyna_d3plot`, `marc_t19`, `nastran_h5`, `nastran_op2`, `radioss_anim`, `radioss_th`, `szplt`, `vtx` and `xplt`, and the `marc` and `radioss` input decks — plus `gltf`, `gmsh22`, `svg` and `tikz` (write-only).

Resolution goes beyond the extension where meshio++ does: an existing `.mesh` starting with `MFEM ` is `mfem` (otherwise `medit`), an existing `.dat` holding a Marc deck is `marc` (otherwise `tecplot`), `.cdb` is `ansysinp` and `.plt` binary or ASCII `tecplot`; fixed file names win over any extension (`z88i1.txt` is `z88`, not `xyz`; `d3plot`, `binout`, `<run>A001`, `<run>T01`); and a path none of that places — an Elmer mesh *directory* — is sniffed by content. Writing Elmer names the format (`"format" : "elmer"`), since a directory that does not exist yet has no content to go by.

`pvtu`/`pvtp`/`pvd` write an index plus a sibling directory named after its stem holding the pieces (`case.pvtu` → `case/case_0000.vtu`). `pcd` and `xyz` are point clouds: only the nodes are written, and a read gives the points back. `gltf` (`.glb`, or `.gltf` with a `.bin` beside it) exports the surface of the mesh for the web — three.js, Blender, model viewers — with per-vertex normals, the nodal data as attributes and, with `color_by`, a baked colour map. `frd` and `nastran_h5` read every increment or result domain as a time step. `xyz` claims `.txt`, `.asc` and `.pts` too, and `nastran_h5` claims `.h5` (a `.post.h5` still resolves to GiD).

`vti`/`vts`/`vtr` (VTK XML ImageData/StructuredGrid/RectilinearGrid) all write only a regular lattice — the output of `Grid`, `voxelize` or `compute_sdf` — though `vtr`'s reader, unlike the other two, accepts a genuinely graded one. `vtm` (VTK XML MultiBlock) has no such requirement: it writes one `.vtu` piece per cell block and reads every piece back as a named `Cell` region.

`cgns`, `h5m`, `hmf`, `med`, `vtkhdf`, `nastran_h5` and the XDMF-HDF data path require an HDF5-enabled meshio++ build; `exodus` requires netCDF; `vtx` ADIOS2 and `szplt` TecIO (both absent from the supported lists of a build without them); `gid` requires zlib to write (see the tip above). A format compiled out is still resolved by extension and reports *why* it is unavailable rather than "unknown format".

`gid` has four on-disk flavours selected with `"gid_mode"` (`auto`/`ascii`/`binary`/`hdf5`/`ascii_zipped`): `ascii` writes a `<stem>.post.msh`/`<stem>.post.res` sibling pair, `binary` one `<stem>.post.bin`, `hdf5` one `<stem>.post.h5`. `auto` never resolves to `ascii_zipped` — no extension can express "zipped".

## ⚠️ Limitations:

- **Serial only.** meshio++ has no MPI, no distributed reader or writer and no communicator. The intended distributed workflow is `partition` with ghost layers, feeding an MPI assembly.
- meshio++'s `undo_green` is **not** exposed. It resolves a refinement's green closure through the colon-namespaced `refine:cell_id`/`refine:parent_id` arrays, and the write-back constraint below means those never reach a model part — the hierarchy is already gone by the time a refined mesh is a `ModelPart`, so any wrapper taking one would fail by construction.
- `remesh` and `decimate` operate on triangle surfaces; `remesh_volume`, `optimize_volume` and `decimate_volume` on tetrahedra. Neither family accepts the other's cells — they are separate operations rather than modes of one, and say so by name when handed the wrong kind.
- Results files carry their data under the solver's own names (`U`, `S11`, `DISP`, ...), which are rarely registered Kratos variables; `"read_field_data"` transfers only the arrays whose name *is* one, and skips the rest with a warning. `radioss_th` and `lsdyna_binout` hold field data on a bare point cloud, so they give nodes and no entities.
- Elmer's halo layer (`write_elmer(..., halo)`) is not exposed: it needs `partition:part` cell data, which no Kratos variable carries — the same reason multi-piece `pvtu` output is not.
- Higher-order Lagrange cells — MFEM meshes of order 3 and above, Z88's 16-node plate — have no Kratos geometry and cannot become entities.
- Kratos has no seven-node triangle, so a `triangle7` block (a Code_Aster `TRIA7`, a VTK/MED `TR7`) cannot become a model part entity.
- Only the halo of a *partitioned* file can be dropped: the application writes a single-piece `pvtu`/`pvtp` (Kratos has no `partition:part` variable to carve by), so multi-piece output is meshio++'s own `partition` → `write_pvtu_pieces_codec` route.
- `Matrix`-valued variables are not written, and a `Matrix`-valued `Properties` entry is skipped with a warning — meshio++ has no representation for it.
- `Properties` entries that are not numeric are left to the materials file: a `Begin Table` curve, and text values such as a constitutive law name (instantiating one needs Kratos's own registry, which meshio++ deliberately does not link). Everything numeric — scalars, integers and vector-valued variables alike — round-trips.
- The C++ `mdpa` reader throws by name on constructs it cannot represent (`Geometries`, `Mesh <id>`, `Constraints`, ...); the `"lenient"` setting downgrades those to a warning and a skip instead.
- `tessellate` (isoparametric subdivision of a higher-order cell) and the `pmsh`/`zarr`/`cae`/`usd` physics-ML dataset formats are Python-only in meshio++ itself — no C++ API exists for them, so neither is reachable from this application.
- **GiD transient output buffers.** meshio++'s `write_gid_series` *pulls* every step through a callback while this IO is *pushed* one step per `PrintOutput`, so the steps are held in memory and the whole series is written in `CloseOutput()` — nothing is on disk before that. Unlike the XDMF writer, which streams. For a long run set `"time_series" : "file_series"` instead, which writes one `<stem>_<label>.post.msh` per step.
- Transient output is held open by the IO. Call `CloseOutput()` before deleting the output of a run: a live writer finalizes on destruction and would recreate the file, and since the series is opened in append mode the next run would then continue a series believed deleted. `MeshioOutputProcess` does this in `ExecuteFinalize`.
- The operations layer only writes a resulting array back onto the output model part when its name matches a registered `Variable<T>` with the right component count: Kratos stores non-historical/historical data keyed by `Variable` objects, not arbitrary strings, so an operation's own invented array names (`attach_quality`'s `"quality:scaled_jacobian"`, for instance) are computed but cannot be retrieved from the output model part — point the operation's own naming setting (`data_calc`'s `"output"`, `data_manage`'s renames, ...) at an existing variable name to get the result back.

## 🗎 Documentation:

The meshio++ library documentation, including a page per format and per operation, is available at [loumalouomega.github.io/meshioplusplus](https://loumalouomega.github.io/meshioplusplus).
