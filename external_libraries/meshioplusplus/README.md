# meshio++ (vendored)

C++ core of **meshio++**, a C++20 port of [meshio](https://github.com/nschloe/meshio)
supporting ~35 mesh file formats. Used by the Kratos core IO class
`MeshioPlusPlusIO` (`kratos/input_output/meshioplusplus_io.h`).

- **Upstream**: `https://github.com/loumalouomega/meshioplusplus` (local development tree: `/home/vicente/src/meshioplusplus`)
- **Vendored version**: `v6.2.1` (`7b77fe3dfeb1b618c0954aa16b10e0c443481651`)
- **License**: MIT (see `LICENSE`; © 2015-2021 Nico Schlömer et al., © 2026 Vicente Mataix Ferrándiz).
  The bundled `cpp/third_party/pugixml` is also MIT licensed (© Arseny Kapoulkine).

## What is vendored

| Path | Content |
|------|---------|
| `cpp/include/meshioplusplus/` | All public headers (mesh backends, format APIs, `kratos_bridge.hpp`, `detail/`) |
| `cpp/src/` | Compiled format readers/writers and the format registry |
| `cpp/third_party/pugixml/` | pugixml (XML parser, compiled into the static lib) |
| `LICENSE` | Upstream MIT license |

**Omitted from upstream**: the pure-Python package, pybind11/WASM bindings, tests,
benchmarks, the Eigen submodule (`MESHIOPLUSPLUS_HAS_EIGEN` is left off — the MED
writer uses its plain transpose fallback) and the upstream `CMakeLists.txt`
(replaced by a Kratos-specific one in this directory).

## Build configuration (Kratos-specific)

The library is built as a STATIC lib `meshioplusplus_core` and linked into
`KratosCore` (see `kratos/CMakeLists.txt`). Configuration highlights:

- `MESHIOPLUSPLUS_MESH_BACKEND_KRATOS`: the `Mesh` type is `KratosMesh`, which
  materializes a Kratos-style `meshioplusplus::ModelPart` (bridged to the real
  `Kratos::ModelPart` via the header-only `kratos_bridge.hpp`).
- `MESHIOPLUSPLUS_HAS_ZLIB`: always on — the Kratos build guarantees zlib.
- `MESHIOPLUSPLUS_HAS_HDF5` (CMake option `KRATOS_MESHIOPLUSPLUS_HDF5`, default `ON`):
  set when HDF5 is found; enables the `med`, `cgns`, `h5m`, `hmf` formats and the
  HDF data path of `xdmf`. Note the MedApplication also links HDF5 — using the
  same shared HDF5 library in one process is fine; mixing *different* HDF5
  versions triggers the usual HDF5 runtime version warning.
- `MESHIOPLUSPLUS_HAS_NETCDF` (CMake option `KRATOS_MESHIOPLUSPLUS_NETCDF`, default `ON`):
  set when netCDF is found; enables the `exodus` format.
- Parallel backend: OpenMP when available (matching Kratos), sequential otherwise.

Formats compiled out due to a missing dependency still resolve by extension and
raise a descriptive error naming the missing dependency
(`meshioplusplus::registry_compiled_out`).

## How to resync with upstream

```bash
UPSTREAM=/path/to/meshioplusplus
DST=external_libraries/meshioplusplus
rm -rf $DST/cpp
cp -r $UPSTREAM/cpp/include $DST/cpp/include
cp -r $UPSTREAM/cpp/src     $DST/cpp/src
mkdir -p $DST/cpp/third_party/pugixml
cp $UPSTREAM/cpp/third_party/pugixml/{pugixml.hpp,pugixml.cpp,pugiconfig.hpp} $DST/cpp/third_party/pugixml/
cp $UPSTREAM/LICENSE $DST/LICENSE
# then update the "Vendored version" above (prefer a tagged release when one exists):
git -C $UPSTREAM describe --tags
git -C $UPSTREAM rev-parse HEAD
```

Do **not** modify the vendored sources directly — fix issues upstream and resync.
