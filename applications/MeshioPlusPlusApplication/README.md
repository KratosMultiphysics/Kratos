# Meshio++ Application

<p align="center">
  <a href="https://github.com/loumalouomega/meshioplusplus"><img alt="meshio++" src="https://raw.githubusercontent.com/loumalouomega/meshioplusplus/master/doc/logo/logo-with-text.svg" width="100%"></a>
</p>

|            **Application**            |                                                                                                    **Description**                                                                                                    |                                              **Status**                                              |                                 **Authors**                                 |
|:-------------------------------------:|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|:----------------------------------------------------------------------------------------------------:|:---------------------------------------------------------------------------:|
| `MeshioPlusPlusApplication`           | The *Meshio++ Application* wraps the [meshio++](https://github.com/loumalouomega/meshioplusplus) library, bringing multi-format mesh input/output (80 readable, 71 writable formats) and a full mesh/data operations layer into *Kratos Multiphysics* | <img src="https://img.shields.io/badge/Status-%F0%9F%94%A7Maintained-blue"  width="300px"> | [*Vicente Mataix Ferrándiz*](mailto:vicente.mataix-ferrandiz@siemens.com) |

The application includes tests to check the proper functioning of the application.

## 😎 Features:

- **Multi-format mesh input/output through `MeshioPlusPlusIO`**

    * *85 readable and 84 writable formats, including the Kratos native `mdpa`, the Kratos CAD JSON of isogeometric analysis, the GiD postprocess format, the ParaView `pvtu`/`pvtp`/`pvd` indexes, VTKHDF, the FEM interchange decks (Elmer, FEBio, Femap, libMesh, MFEM, Patran, Z88), LS-DYNA and Code_Aster decks, the solver results files of Abaqus, Ansys, Nastran, LS-DYNA, Marc, Radioss, FEBio and DOLFINx (all written since meshio++ v17.3.0), CalculiX and MSC Nastran HDF5 results, STEP and IGES CAD models, point clouds and glTF*

    * *Format resolved explicitly or inferred from the file extension*

    * *Sub model parts mapped to meshio++ named regions, nesting included*

    * *Registered entity names (`SmallDisplacementElement3D4N`, ...) preserved across a round trip, instead of degrading to the generic cell-type name*

    * *`Properties` material data carried in both directions, so a `mdpa` round trip keeps the values and not just the ids*

    * *`"time_step"` (negative counts from the end) selects one step of a multi-step file for the formats meshio++ reads selectively — mdpa, med, exodus, gmsh, tecplot, ensight, cgns, openfoam, frd, unv, nastran_h5, vtkhdf and pvd — through a cheap native metadata reader where one exists; any other format ignores it and is always read whole*

    * *`"lenient"` downgrades an mdpa/med/vtkhdf construct the reader cannot represent from an error to a warning and a skip, instead of throwing*
    * *`"mdpa_keep_extra_blocks"` reads an mdpa deck through meshio++'s `MdpaInfo` side channel (v16.27.0), so a deck with top-level tables, `Geometries`, `Mesh <id>` blocks, sub model part data, text `ModelPartData` or a `Constraints` block is read instead of refused; a model part cannot hold those blocks, so `MeshioPlusPlusIO.GetMdpaInfo()` reports what was found - including, since v16.29.0, how many geometry and constraint ids the sub model parts list under `SubModelPartGeometries`/`SubModelPartConstraints` - and a warning is logged*

    * *`"keep_brep"` reads a Kratos CAD JSON (`cad_json`), STEP or IGES file through meshio++'s `BrepInfo` side channel (v16.42.0), and `"mdpa_cad_json"` (`"auto"`, or a path) does the same for the `<stem>.cad.json` — or the `CadIoModeler` geometry file of a `ProjectParameters.json` — that goes with an mdpa deck. The model part receives the geometry sampled into cells (faces as surface cells, edges as lines; STEP/IGES are triangulated, see `"cad_deflection"` and `"cad_angular_deflection"`); a trimmed NURBS face is not a Kratos entity, so `MeshioPlusPlusIO.GetBrepInfo()` reports what the file holds (faces, trimmed faces, edges, vertices, patches), `WriteBrep(file)` writes the exact geometry back as Kratos CAD JSON (a lossless CAD JSON/STEP/IGES conversion), and it stays reachable through `MeshioPlusPlusIgaUtilities`*

    * *A partitioned file (`vtkhdf`, `pvtu`, `pvtp`, `pvd`) is read whole, one sub model part per piece; `"select_piece"` with `"piece"` (negative counts from the end) reads one piece instead, and `"ghosts" : "drop"` removes the halo cells a `pvtu`/`pvtp`/`pvd` carries as `vtkGhostType`*

    * *Transient output of the solver results formats (`xplt`, `lsdyna_d3plot`, `ansys_rst`, `marc_t19`, `nastran_h5`, `nastran_op2`, `abaqus_fil`, `lsdyna_binout`, `radioss_th`, `radioss_anim`, `frd`, `vtx`, `szplt` and the `openfoam` time directories) goes through meshio++'s own series writer in the default `"time_series" : "automatic"` mode: the mesh once and one state per `WriteModelPart` call, in one file (`radioss_anim` is one file per state by design, `openfoam` a time directory per step), closed by `CloseOutput()`. `"time_series" : "file_series"` still writes a file per step. `"d3plot_double"` and `"d3plot_max_file_bytes"` set the `lsdyna_d3plot` word size and the size at which the family is split*

    * *`"openfoam_region"` also selects the region on a write: `constant/<region>/polyMesh` of a multi-region case, with the fields under `<time>/<region>/`*

    * *Solver results files (`abaqus_fil`, `ansys_rst`, `nastran_op2`, `lsdyna_d3plot`, `marc_t19`, `xplt`, `vtx`, ...) read one step per `"time_step"`, `GetTimeValues()` lists the steps, and `"read_field_data" : true` carries their arrays onto the registered variables of the same name*

    * *A format meshio++ cannot place by extension or file name — an Elmer mesh directory, an extensionless deck — is identified by its content*

    * *Companion field files: `"mfem_grid_functions"` and `"patran_result_files"` (`{"name", "path"}` lists) read the fields an MFEM `.gf` or Patran `.nod`/`.dis`/`.els` file holds; `"mfem_grid_functions_write"` writes the nodal data next to an MFEM mesh as `<stem>.<name>.gf`, ready for an MFEM solver; `"z88_results"` and `"z88_stubs"` control a Z88 deck's results and input stubs*

    * *Format-specific writers: `"gltf_settings"` (`container`, `up_axis`, `split_angle`, `color_by` with `cmap`/`range`/`component`, ...) for glTF, `"openfoam_label_bits"`/`"openfoam_scalar_bits"` with `"file_format" : "binary"` for a binary polyMesh, `"pcd_compressed"`/`"pcd_float64_points"` for PCL's `binary_compressed`, `"vtkhdf_gzip_level"`, `"ensight_fortran"` for Fortran-record binary EnSight, `"radioss_stubs"` for a starter deck OpenRadioss accepts on its own (placeholder materials and properties); every other format meshio++'s ascii/binary switch covers honours `"file_format"`, and `"file_format" : "raw_appended"` writes a `.vtu` with one raw `<AppendedData>` section (no base64: the fastest VTU to write and read; any other format refuses it)*

    * *`"openfoam_region"` selects one region of an OpenFOAM case with no single `constant/polyMesh`, but `constant/<region>/polyMesh` per region; FLAC3D `ZGROUP`/`FGROUP` cell groups round-trip as named regions instead of always deferring to the Python fallback; a gmsh 4.1 export now allocates a physical tag for a named `Cell` region that has none of its own, instead of dropping it*

- **Transient output**

    * *XDMF temporal collections and VTKHDF `Steps` written in a single file, ParaView `.pvd` collections (one `.vtu` per step, index rewritten after every step), Femap neutral files (the mesh once, one output set per step), GiD multi-step `.post` series, file series for every other format*

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

    * *`check_quality` — tests thresholds against every cell's quality metrics (`"require" : "scaled_jacobian >= 0.2; aspect_ratio <= 5 @ 1%"`, `"max_inverted"`, `"max_degenerate"`) and reports `passed`, a `summary` and per-check violations and worst cell*

    * *`feature_edges` — the sharp, open, non-manifold and inconsistently wound edges of a surface (or of a volume's skin) as a line mesh (`"feature_angle"`, `"include_feature_edges"`, `"include_boundary_edges"`, `"include_non_manifold_edges"`, `"include_inconsistent_edges"`, `"region"`)*

    * *`segment` — splits a surface (or a volume's skin) into smooth patches by region growing bounded by its feature edges (`"feature_angle"`, `"normal_tolerance"`, `"curvature_tolerance"`, `"min_patch_cells"`, `"region"`); every patch comes back as a sub model part named `<"patch_prefix">:<k>`, and the report carries `number_of_patches`, `merged_patches` and the edge counts. `fit_primitives` — fits a plane, sphere, cylinder, cone or torus to each of those patches (`"fit_tolerance"`, `"fit_relative"`, `"min_inlier_fraction"`, `"fit_types"`, `"fit_max_iterations"`) and reports how many of each it found and how many it could not explain. meshio++'s `segment:patch` and `fit:*` arrays are colon-namespaced, so they never reach the model part; `fit_primitives` therefore runs the segmentation itself first, with the `segment` settings (meshio++ v16.39.0)*

    * *`tessellate` — the isoparametric subdivision of the quadratic cells (`triangle6`, `quad8`, `quad9`, `tetra10`, `hexahedron27`: Kratos' `Triangle2D6`, `Quadrilateral2D8/9`, `Tetrahedra3D10`, `Hexahedra3D27`) into linear triangles and tetrahedra, interpolating the nodal data (`"tessellate_levels"`, `"tessellate_curved"`, `"tessellate_fields"`, `"tessellate_record_stencil"`); every other cell passes through. The way to visualise or export a quadratic mesh to a format with linear cells only (meshio++ v16.42.0, in the C++ core since then)*

    * *`edit_regions` — union, intersection, difference, rename, retag and delete of the named regions (sub model parts) through an `"edits"` list; named Side regions also survive `refine`, `subdivide`, `convert_cells`, `decimate`, `repair` and `agglomerate` since meshio++ v16.27.0*

    * *`agglomerate` gains `"merge_coplanar_faces"`, `"coplanar_angle"` and `"min_sphericity"`, and reports `faces_merged` and `candidates_rejected`*

    * *`MeshioPlusPlusMeshOperations.HausdorffDistance` (sampled Hausdorff distance between two surfaces), `MatchPeriodicNodes` (the master node an affine transform maps each slave-region node onto, for periodic conditions) and `BlendSteps` (linear blend of two steps of one mesh). meshio++'s `resample_sequence` is driven by a file sequence and a pipeline rather than model parts, so it is not exposed; `BlendSteps` is its per-pair kernel*

- **Geometry fitting and deformation**

    * *`Shrinkwrap` (`MeshioPlusPlusMeshOperations.Shrinkwrap`) — projects one mesh's points onto a target surface, optionally offset along the hit feature's pseudonormal; one projection, not an iteration, so a fit rather than a smoothing. Needs two independent meshes, like `Interpolate`*

    * *`sobolev_deform` — Sobolev-filters a raw nodal displacement field before moving the mesh's points by it, damping high-frequency noise a caller already computed (a `"length_scale"` of `0.0`, the default, applies it directly with no solve at all)*

- **Skin to CAD and isogeometric models**

    * *`fit_cad` — turns a closed, consistently wound triangle skin into a CAD solid: it runs `segment` and `fit_primitives` itself (their settings apply), builds plane, sphere, cylinder, cone and torus faces from the accepted patches, plus B-spline, rolling-ball fillet and blend, chamfer, conic and vertex-blend and helical thread faces where `"cad_freeform"`, `"cad_fillets"` and `"cad_threads"` allow (the rest come back as planar facets, never silently wrong), sews and checks the solid, and writes it as AP214 STEP (`"cad_step_path"`) and/or Kratos CAD JSON (`"cad_json_path"`). The destination model part receives the solid re-meshed at `"cad_mesh_deflection"`; the report carries validity, the tolerance applied, the volume, the volume ratio against the skin, the maximum deviation and a count per face kind. Needs a meshio++ built with OpenCASCADE (`-DMESHIOPLUSPLUS_WITH_OCCT=ON`); without it the operation refuses by name (meshio++ v16.41.0–v17.2.0)*

    * *`MeshioPlusPlusIgaUtilities` — `Describe`, `Sample`, `Refine` and `Evaluate` for the NURBS models of Kratos CAD JSON, STEP and IGES files and of MFEM NURBS meshes (meshio++ v16.43.0). A NURBS model is not a `ModelPart`, so these take the model from a file: `Sample` fills a model part with it sampled into Lagrange cells (`"sample_subdivisions"`, `"trim"`, `"edges"`), `Refine` writes an exactly h- (`"subdivisions"`) and p- (`"elevate"`) refined model back as CAD JSON, or as an MFEM NURBS mesh, and `Evaluate` returns the points, first derivatives, normals and — for a trimmed face — the inside-the-trimming-loops flag at given parametric coordinates*

- **Regular grids and signed distance**

    * *`Grid` — a regular hexahedron lattice from nothing: the one *generator* in the library*

    * *`voxelize` — a lattice around a mesh, keeping the whole bounding box, only the cells a surface passes through, or only those inside it*

    * *`NeighborPairs` — the radius or k-nearest neighbour pairs of a model part's nodes (the cell-lattice search behind meshio++'s `proximity_graph`, optionally periodic), as node-id pairs*

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

Four features need an optional dependency, all off by default: `-DMESHIOPLUSPLUS_WITH_ADIOS2=ON` for DOLFINx `vtx` (`.bp`), `-DMESHIOPLUSPLUS_WITH_TECIO=ON -DTECIO_ROOT=<dir>` for Tecplot `szplt`, `-DMESHIOPLUSPLUS_WITH_BZIP2=ON` for bzip2-compressed libMesh (`.xda.bz2`), and `-DMESHIOPLUSPLUS_WITH_OCCT=ON` (OpenCASCADE 7.6 or newer, linked dynamically) for reading STEP and IGES and for the `fit_cad` operation. The configure line of this application reports which ones the install carries.

> [!WARNING]
> Link meshio++ against an ADIOS2 whose libraries use the *shared* libstdc++. Some distribution packages (seen with `adios2 2.12.1-6.1`) link libstdc++ statically and export its symbols, and in a Python process those interpose on Kratos' own (`std::locale` freed by the wrong copy: `free(): invalid size`). `readelf -d libadios2_cxx.so | grep libstdc++` must list it as `NEEDED`; otherwise build ADIOS2 from source (serial, C++ only, `-DBUILD_SHARED_LIBS=ON`) and point `ADIOS2_DIR` at it.

then point the Kratos configure at it:

```bash
-Dmeshioplusplus_DIR=<prefix>/lib/cmake/meshioplusplus
```

> [!IMPORTANT]
> The application requires **meshio++ ABI 30**, release **v17.3.0 or newer**. ABI 30 begins at v17.2.0, but the release floor is higher because the solver results writers this application now serves through the registry arrive with v17.3.0; the other entry points it calls were added along the way: `tessellate` in the C++ core, the `cad_json`/`iges` formats, STEP reading and the `BrepInfo` overloads of `read_mdpa` and the `NurbsInfo` ones of `read_mfem`/`write_mfem` (v16.42.0), `fit_cad` (v16.41.0, grown by v16.42.0–v17.2.0), `iga_evaluate`/`iga_refine`/`iga_sample` (v16.43.0), the `segment` and `fit_primitives` operations (v16.39.0), `resolve_write_format`, the Marc and Radioss writers, `FemapSeriesWriter` and Fortran EnSight (v16.17.0), raw appended VTU (v16.21.0), `MdpaInfo`'s full-deck round trip (v16.27.0), its sub model part geometry and constraint ids (v16.29.0), and `check_quality`, `feature_edges`, `hausdorff_distance`, `edit_regions`, `match_periodic_nodes` and `blend_steps` (v16.23.0–v16.25.0). The pin is on `MESHIOPLUSPLUS_ABI_VERSION` rather than the release version: that counter moves only when the installed headers stop being compatible with an already-compiled consumer, so a release that cannot affect this application needs no rebuild (v12.0.0 is a pure version-number bump over v11.6.0 and holds the ABI at 13, as v10.20.0 did over v10.19.0 at ABI 11). It is a single ABI rather than a window on purpose — supporting several meant one `MESHIOPLUSPLUS_VERSION_AT_LEAST` guard per feature, and the operations exposed here would have needed one each. Getting it wrong is not silent — the C++ variants' `SOVERSION` is the ABI version and every translation unit carries a link-time sentinel naming it, so a skew is a link error rather than memory corruption. See meshio++'s [`doc/abi.md`](https://github.com/loumalouomega/meshioplusplus/blob/master/doc/abi.md) for the criterion.

> [!NOTE]
> Upgrading from any older meshio++ requires relinking this application once: the C++ variants' `SOVERSION` is the ABI version, so the needed library is `libmeshioplusplus_core_kratos.so.30`. Install the new meshio++ wholesale rather than part-upgrading — its headers reference a link-time sentinel an older library does not define, which fails closed rather than silently mixing.

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

**Read (85, less any format compiled out — 84 without TecIO, 82 without TecIO and OpenCASCADE):** `abaqus` `abaqus_fil` `ansys` `ansys_rst` `ansys_rst_cyclic` `ansysinp` `avsucd` `cad_json` `cgns` `code_aster` `dex` `dolfin` `elmer` `ensight` `exodus` `febio` `femap` `flac3d` `flux` `frd` `freefem` `gid` `gltf` `gmsh` `h5m` `hmf` `iges` `ip` `libmesh` `lsdyna` `lsdyna_binout` `lsdyna_d3plot` `marc` `marc_t19` `mdpa` `med` `medit` `mfem` `mff` `mfm` `mphbin` `mphtxt` `nastran` `nastran_h5` `nastran_op2` `netgen` `neuroglancer` `obj` `off` `openfoam` `patran` `pcd` `permas` `ply` `pvd` `pvtp` `pvtu` `radioss` `radioss_anim` `radioss_th` `step` `stl` `su2` `svg` `szplt` `tecplot` `tetgen` `tikz` `triangle` `ugrid` `unv` `vti` `vtk` `vtkhdf` `vtm` `vtp` `vtr` `vts` `vtu` `vtx` `wkt` `xdmf` `xplt` `xyz` `z88`

**Write (84, less any format compiled out — 83 without TecIO):** the same set minus the read-only ones — `ansys_rst_cyclic` (an expansion, selected by name) and `iges` — plus `gmsh22` (write-only). The solver results files (`abaqus_fil`, `ansys_rst`, `lsdyna_binout`, `lsdyna_d3plot`, `marc_t19`, `nastran_h5`, `nastran_op2`, `radioss_anim`, `szplt`, `vtx`, `xplt`) are written since meshio++ v17.3.0, `frd` and `radioss_th` since v16.39.0. `gltf` is read back since v17.1.0 and `step` since v16.42.0 (OpenCASCADE only), `svg` and `tikz` since v16.39.0 (as a planar drawing, `z = 0`); `step` writes a closed surface — or a volume's skin — as a faceted AP214 solid (v16.39.0), `cad_json` writes `quad` cells as bilinear faces and `neuroglancer` triangles only, with no extension to resolve it by (name the format). The `marc` and `radioss` decks are written since meshio++ v16.17.0 (`.rad` is Radioss's; a Marc deck is written by naming `"format" : "marc"`).

Resolution goes beyond the extension where meshio++ does: an existing `.mesh` starting with `MFEM ` is `mfem` (otherwise `medit`), an existing `.dat` holding a Marc deck is `marc` (otherwise `tecplot`), `.cdb` is `ansysinp` and `.plt` binary or ASCII `tecplot`; fixed file names win over any extension (`z88i1.txt` is `z88`, not `xyz`; `d3plot`, `binout`, `<run>A001`, `<run>T01`); and a path none of that places — an Elmer mesh *directory* — is sniffed by content. Writing Elmer names the format (`"format" : "elmer"`), since a directory that does not exist yet has no content to go by. A *write* never looks at what the path already holds: `.dat` is Tecplot's and `.mesh` Medit's even over a Marc deck or an MFEM mesh — name `"marc"`/`"mfem"` to write those.

`pvtu`/`pvtp`/`pvd` write an index plus a sibling directory named after its stem holding the pieces (`case.pvtu` → `case/case_0000.vtu`). `pcd` and `xyz` are point clouds: only the nodes are written, and a read gives the points back. `gltf` (`.glb`, or `.gltf` with a `.bin` beside it) exports the surface of the mesh for the web — three.js, Blender, model viewers — with per-vertex normals, the nodal data as attributes and, with `color_by`, a baked colour map. `frd` and `nastran_h5` read every increment or result domain as a time step. `xyz` claims `.txt`, `.asc` and `.pts` too, and `nastran_h5` claims `.h5` (a `.post.h5` still resolves to GiD).

`vti`/`vts`/`vtr` (VTK XML ImageData/StructuredGrid/RectilinearGrid) all write only a regular lattice — the output of `Grid`, `voxelize` or `compute_sdf` — though `vtr`'s reader, unlike the other two, accepts a genuinely graded one. `vtm` (VTK XML MultiBlock) has no such requirement: it writes one `.vtu` piece per cell block and reads every piece back as a named `Cell` region.

`cgns`, `h5m`, `hmf`, `med`, `vtkhdf`, `nastran_h5` and the XDMF-HDF data path require an HDF5-enabled meshio++ build; `exodus` requires netCDF; `vtx` ADIOS2, `szplt` TecIO, and reading `step` and `iges` OpenCASCADE (all absent from the supported lists of a build without them); `gid` requires zlib to write (see the tip above). A format compiled out is still resolved by extension and reports *why* it is unavailable rather than "unknown format".

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
- The C++ `mdpa` reader throws by name on constructs it cannot represent (`Geometries`, `Mesh <id>`, `Constraints`, ...); the `"lenient"` setting downgrades those to a warning and a skip instead, and `"mdpa_keep_extra_blocks"` reads them into meshio++'s `MdpaInfo` (reported by `GetMdpaInfo()`, not transferred to the model part).
- The `pmsh`/`zarr`/`cae`/`usd` physics-ML dataset formats are Python-only in meshio++ itself — no C++ API exists for them, so none is reachable from this application.
- A NURBS model (the trimmed faces of `cad_json`, `step` and `iges`, the patches of an MFEM NURBS mesh) is not a `ModelPart` entity: a read gives its sampled cells, `GetBrepInfo()` a summary, and `MeshioPlusPlusIgaUtilities` the exact geometry from the file. A B-rep is not rebuilt from a model part (a `cad_json` write of one gives bilinear faces from `quad` cells; `WriteBrep()` writes what a read kept), and `fit_cad`'s `fit_cad:*` cell data are colon-namespaced, so they do not reach the destination.
- `step` and `iges` are read only by a meshio++ built with OpenCASCADE; otherwise they are listed as unavailable and say so. The VTK Bézier cell types (`VTK_BEZIER_*`) have no Kratos geometry either.
- **GiD transient output buffers.** meshio++'s `write_gid_series` *pulls* every step through a callback while this IO is *pushed* one step per `PrintOutput`, so the steps are held in memory and the whole series is written in `CloseOutput()` — nothing is on disk before that. Unlike the XDMF writer, which streams. For a long run set `"time_series" : "file_series"` instead, which writes one `<stem>_<label>.post.msh` per step.
- Transient output is held open by the IO. Call `CloseOutput()` before deleting the output of a run: a live writer finalizes on destruction and would recreate the file, and since the series is opened in append mode the next run would then continue a series believed deleted. `MeshioOutputProcess` does this in `ExecuteFinalize`.
- The operations layer only writes a resulting array back onto the output model part when its name matches a registered `Variable<T>` with the right component count: Kratos stores non-historical/historical data keyed by `Variable` objects, not arbitrary strings, so an operation's own invented array names (`attach_quality`'s `"quality:scaled_jacobian"`, for instance) are computed but cannot be retrieved from the output model part — point the operation's own naming setting (`data_calc`'s `"output"`, `data_manage`'s renames, ...) at an existing variable name to get the result back.

## 📝 Pending TODO:

- [ ] **Expose `fit_cad`'s per-cell data.** The `fit_cad:face`, `fit_cad:type` and `fit_cad:patch` arrays (and `segment:patch`) are colon-namespaced, so they never reach the destination model part; only the re-meshed solid and the report do. A sub model part per face (or per face kind) would carry them.
- [ ] **Neuroglancer by extension.** The format has no extension and holds triangles only, so it must be selected with `"format" : "neuroglancer"`; meshio++ itself registers no extension for it.
- [ ] **Cell types without a Kratos geometry.** Cubic and higher Lagrange cells (e.g. a degree-3 `MeshioPlusPlusIgaUtilities.Sample`), seven-node triangles and VTK Bézier cells are refused by name; they could be mapped to Kratos geometries, or tessellated automatically on read.
- [ ] **Decomposed OpenFOAM output.** `write_openfoam_decomposed` (`processorN/` cases) is not exposed: it needs the `partition:part` cell data, which no Kratos variable carries. An `"openfoam_decomposed"` setting driven by the model part's partition (or a sub model part per processor) would do.
- [ ] **Write the B-rep with a deck.** `WriteBrep()` writes the exact geometry of the last read as a CAD JSON file, but the sidecar is not written together with an `.mdpa` (meshio++'s `write_mdpa(..., BrepInfo)`).
- [ ] **Resume a transient solver-results series.** Like PVD and Femap, the series of the solver results formats starts afresh with each `MeshioPlusPlusIO` instance; only XDMF and VTKHDF continue a file left by a previous run.

## 🗎 Documentation:

The meshio++ library documentation, including a page per format and per operation, is available at [loumalouomega.github.io/meshioplusplus](https://loumalouomega.github.io/meshioplusplus).
