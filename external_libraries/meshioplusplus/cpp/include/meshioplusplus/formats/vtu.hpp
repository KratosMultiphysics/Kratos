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
 * @file vtu.hpp
 * @brief VTK XML UnstructuredGrid (.vtu) C++ reader/writer, ascii and inline
 *        binary (uncompressed or zlib), single `<Piece>` only.
 *
 * A `.vtu` file is `<VTKFile type="UnstructuredGrid" ...><UnstructuredGrid>
 * <Piece ...><Points>/<Cells>(connectivity,offsets,types[,faces,
 * faceoffsets])/<PointData>/<CellData></Piece><FieldData></UnstructuredGrid>
 * [<AppendedData>...]</VTKFile>`. Cell reconstruction from
 * connectivity+offsets+types shares the exact helper used by the VTK 5.1
 * reader (`detail/vtk_cells.hpp`) — zero-copy-friendly when the connectivity
 * block is contiguous and node order is identity.
 *
 * **Binary encoding** matches VTK's own convention exactly so files
 * round-trip byte-for-byte with other VTK tools: uncompressed is
 * `base64(header[header_type: total_nbytes] + raw_bytes)`; zlib-compressed
 * is `base64(header[nblocks, blocksize=32768, last_block_size,
 * csize_0..csize_{n-1}])` followed by a **separate** base64 blob of
 * `concat(compressed_block_0..n-1)` — header fields use the file's declared
 * `header_type` dtype throughout. The C++ writer always declares
 * `byte_order="LittleEndian"` (unlike the Python writer, which records the
 * host's native order).
 *
 * The C++ core explicitly does **not** implement several paths, each of
 * which throws ReadError/WriteError to force the Python fallback:
 * lzma compression; `<AppendedData>` (raw/appended binary, including the
 * regex-based manual XML-repair the Python reader falls back to when
 * appended raw bytes break XML parsing); polyhedron cells (both read and
 * write — they also cannot mix with other cell types, a Python-side
 * ValueError); multiple `<Piece>` elements (the Python reader concatenates
 * them, C++ requires exactly one); and any non-default `header_type` other
 * than `UInt32`.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Write a mesh as a `.vtu` file.
 *
 * Emits `<Points>`, `<Cells>` (connectivity/offsets/types, using the full
 * VTK cell-type set including Lagrange high-order cells), `<PointData>`/
 * `<CellData>` (generic key mapping; `cell_sets` round-trip as extra data
 * arrays), and `<FieldData>`. 2D points are silently padded to 3D. Byte
 * order is always declared `LittleEndian`.
 *
 * @param rPath filesystem path to write
 * @param rMesh the mesh to write
 * @param binary true for base64-encoded binary DataArrays, false for inline
 *        ascii text
 * @param zlib true additionally zlib-compresses binary DataArrays in
 *        32768-byte blocks (ignored when `binary` is false); lzma is not
 *        supported by the C++ writer
 * @throws WriteError if the mesh contains polyhedron cells, an unknown cell
 *         type, or the output path cannot be opened
 * @note cell_data key handled specially: `cell_sets` become extra
 *       `<CellData>` arrays (VTU has no native set concept).
 */
void write_vtu(const std::string& rPath, const Mesh& rMesh, bool binary, bool zlib);

/**
 * @brief Read a `.vtu` file.
 *
 * Parses the single `<Piece>`'s `<Points>`/`<Cells>`/`<PointData>`/
 * `<CellData>` and the grid's `<FieldData>`, decoding ascii, uncompressed
 * binary, or zlib-compressed binary `DataArray` payloads per the encoding
 * described above.
 *
 * @param rPath filesystem path to read
 * @return the read Mesh
 * @throws ReadError if the file uses lzma compression, an `<AppendedData>`
 *         section, more than one `<Piece>`, polyhedron cells, or a
 *         non-`UInt32` `header_type` — the shim then falls back to the
 *         Python reader, which supports all of these.
 * @note `<FieldData>` -> `mesh.field_data`; `<PointData>`/`<CellData>` map
 *       generically to `point_data`/`cell_data`.
 */
Mesh read_vtu(const std::string& rPath);

}  // namespace meshioplusplus
