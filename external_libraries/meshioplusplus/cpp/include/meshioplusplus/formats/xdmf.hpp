//  ██████   ██████ ██████████  █████████  █████   █████ █████    ███████
// ░░██████ ██████ ░░███░░░░░█ ███░░░░░███░░███   ░░███ ░░███   ███░░░░░███      ███         ███
//  ░███░█████░███  ░███  █ ░ ░███    ░░░  ░███    ░███  ░███  ███     ░░███    ░███        ░███
//  ░███░░███ ░███  ░██████   ░░█████████  ░███████████  ░███ ░███      ░███ ███████████ ███████████
//  ░███ ░░░  ░███  ░███░░█    ░░░░░░░░███ ░███░░░░░███  ░███ ░███      ░███░░░░░███░░░ ░░░░░███░░░
//  ░███      ░███  ░███ ░   █ ███    ░███ ░███    ░███  ░███ ░░███     ███     ░███        ░███
//  █████     █████ ██████████░░█████████  █████   █████ █████ ░░░███████░      ░░░         ░░░
// ░░░░░     ░░░░░ ░░░░░░░░░░  ░░░░░░░░░  ░░░░░   ░░░░░ ░░░░░    ░░░░░░░
//
//
//  License:         MIT License
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//
#pragma once

/**
 * @file xdmf.hpp
 * @brief XDMF3 (.xdmf/.xmf) C++ reader/writer — "light data" XML with
 *        XML/Binary/HDF "heavy data" DataItem payloads.
 *
 * An XDMF file is `<Xdmf Version="3.x"><Domain><Grid><Topology
 * Type|TopologyType=".." NumberOfElements=".."><DataItem
 * DataType|NumberType="Int|UInt|Float" Precision="1|2|4|8"
 * Dimensions=".." Format="XML|Binary|HDF">...</DataItem></Topology>
 * <Geometry Type|GeometryType="X|XY|XYZ"><DataItem .../></Geometry>
 * <Attribute Name=".." AttributeType="Scalar|Vector|Tensor|Tensor6|Matrix"
 * Center="Node|Cell|Grid"><DataItem .../></Attribute></Grid></Domain>
 * </Xdmf>`. Both XDMF2 and XDMF3 exist in the wild (dispatched by the major
 * version digit in the root `Version` attribute), but **the C++ core only
 * implements version 3** — any `Version="2.x"` file throws ReadError and
 * falls back to Python. XDMF3 accepts either `Type=` or `TopologyType=`/
 * `GeometryType=` but errors if both are given on the same element.
 *
 * `Format="XML"` DataItem text is whitespace-separated inline numbers;
 * `Format="Binary"` text is a raw-binary sibling file path; `Format="HDF"`
 * text is `"<file.h5>:/path/to/dataset"` (resolved relative to the `.xdmf`
 * file). The HDF path is handled by the C++ core **only** when built with
 * `MESHIOPLUSPLUS_HAS_HDF5` and `compression` is `None`/`"gzip"`; otherwise
 * it throws and the Python `h5py` fallback takes over — this includes the
 * always-Python case of a non-HDF5 build (the `#ifdef`-guarded HDF code
 * compiles to an empty/throwing path).
 *
 * `Mixed` topology encodes a flat array of `(xdmf_type_index, node0, ...)`
 * tuples concatenated across all cells (shared type-index table with the
 * per-type `TopologyType` names, e.g. `0x6`=tetra, `0x9`=hexahedron,
 * `0x26`=tetra10). A `line`/`Polyline` entry in a Mixed array carries an
 * extra "point count" field that **must equal exactly 2** — anything else
 * throws ReadError. The C++ type table is a strict subset of the Python
 * one: it covers through `hexahedron27` but omits the higher-order
 * `hexahedron64`..`hexahedron1331` types, and does not implement
 * `Reference="XML"`/XPath DataItem references or the XDMF2-only
 * `Information`-based `field_data` — all of these throw and fall back to
 * Python. Points are restricted to dimension <=3 on write.
 *
 * Temporal XDMF (`TimeSeriesWriter`/`TimeSeriesReader`) is unrelated to this
 * header and remains pure Python regardless of the C++ core.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Write a mesh as an XDMF3 file.
 *
 * Emits `<Topology>` (single-type or `Mixed`), `<Geometry Type="XYZ">` (or
 * `X`/`XY` per point dimension), `<Attribute>` elements for point_data/
 * cell_data, with DataItem payloads stored per `data_format`.
 *
 * @param rPath filesystem path to write (companion `.h5`/`.bin` sibling
 *        files are written alongside it for `"HDF"`/`"Binary"`)
 * @param rMesh the mesh to write
 * @param rDataFormat one of `"XML"` (inline text), `"Binary"` (external raw
 *        sibling files), or `"HDF"` (companion `.h5` file, requires an
 *        HDF5-enabled build)
 * @param gzip_level gzip compression level for `"HDF"` DataItems; `-1`
 *        (default) means uncompressed. Ignored for `"XML"`/`"Binary"`.
 * @throws WriteError if the mesh mixes cell types that cannot share one
 *         `Topology` block, points exceed dimension 3, `data_format="HDF"`
 *         is requested on a build without HDF5 support, or `data_format` is
 *         otherwise unrecognized
 * @note point_data/cell_data map generically to `<Attribute Center="Node"|
 *       "Cell">` elements, keyed by the raw attribute name.
 */
void write_xdmf(const std::string& rPath, const Mesh& rMesh, const std::string& rDataFormat,
                int gzip_level = -1);

/**
 * @brief Read an XDMF3 file's first `<Grid>`.
 *
 * Parses `<Topology>` (resolving `Mixed` via the numeric type-index table),
 * `<Geometry>`, and `<Attribute>` elements, decoding each `<DataItem>`
 * according to its `Format` (`XML` inline, `Binary` external file, or `HDF`
 * companion dataset when built with HDF5 support).
 *
 * @param rPath filesystem path to read
 * @return the read Mesh
 * @throws ReadError if the file is XDMF2 (`Version="2.x"`), uses a
 *         `Reference` DataItem attribute, an XDMF2 `Information` field-data
 *         block, a Mixed `Polyline` entry with a point count other than 2, a
 *         cell type outside the C++ type table (e.g. `hexahedron64`+), or a
 *         `Format="HDF"` DataItem on a build without HDF5 support — the shim
 *         then falls back to the Python/`h5py` reader.
 * @note `<Attribute>` elements map generically to `point_data`/`cell_data`.
 */
Mesh read_xdmf(const std::string& rPath);

}  // namespace meshioplusplus
