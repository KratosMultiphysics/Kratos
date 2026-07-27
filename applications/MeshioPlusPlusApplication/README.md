# Meshio++ Application

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

    * *`data_*` — expression calculator, conditioning, rename/drop/keep, point ⇄ cell averaging*

    * *`interpolate` — cross-mesh field transfer*

## 🛠️ Building:

meshio++ is consumed as a normal external dependency; **nothing is vendored into Kratos**. Build and install it with the C++ API and the `KRATOS` mesh backend:

```bash
cmake -S <meshioplusplus> -B build \
  -DMESHIOPLUSPLUS_INSTALL_CPP=ON \
  -DMESHIOPLUSPLUS_INSTALL_CPP_BACKENDS="KRATOS" \
  -DMESHIOPLUSPLUS_MESH_BACKEND=KRATOS \
  -DMESHIOPLUSPLUS_BUILD_PYTHON=OFF \
  -DBUILD_SHARED_LIBS=ON
cmake --build build && cmake --install build --prefix <prefix>
```

then point the Kratos configure at it:

```bash
-Dmeshioplusplus_DIR=<prefix>/lib/cmake/meshioplusplus
```

> [!IMPORTANT]
> The application requires **meshio++ 9.1 exactly**. The meshio++ C++ API makes no ABI promise — `Mesh`, `ModelPart` and `GeometricalEntity` are header-defined types whose layout changes with the headers — so the application must be rebuilt whenever meshio++ is.

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
- `Matrix`-valued variables are not written.
- **`Properties` values do not round-trip yet.** meshio++ 9.1.0 parses `Begin Properties` bodies into an `MdpaInfo` side channel, but its format registry — which this IO reads through — does not thread that through, so a read produces the properties *ids* referenced by the entities and no material data. The typed assignment is already implemented on the Kratos side (`ApplyMeshioProperty`) and starts working as soon as upstream carries the values on the mesh. Set material data from the materials file meanwhile.
- Reading a format that stores the properties id as a cell tag (the `gmsh:physical` convention, which the `mdpa` writer follows) produces an extra `gmsh_physical_<id>` sub model part alongside the real ones.
- The C++ `mdpa` reader throws by name on constructs it cannot represent (`Geometries`, `Mesh <id>`, `Constraints`, ...); `ReadOptions::mLenient` downgrades those to a warning and a skip.

## 🗎 Documentation:

The meshio++ library documentation, including a page per format and per operation, is available at [loumalouomega.github.io/meshioplusplus](https://loumalouomega.github.io/meshioplusplus).
