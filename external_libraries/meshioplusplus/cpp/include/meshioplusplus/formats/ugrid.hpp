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
 * @file ugrid.hpp
 * @brief AFLR UGRID (.ugrid) C++ reader/writer, ascii and every binary
 *        flavour.
 *
 * The byte layout is entirely determined by the file's **penultimate**
 * filename suffix, e.g. `foo.lb8.ugrid` -> flavour `lb8`: no suffix = ascii
 * (native types); `b8l`/`b8`/`b4` = C-layout, big-endian, `{8,8,4}`-byte
 * floats and `{8,4,4}`-byte ints respectively; `lb8l`/`lb8`/`lb4` = the same
 * but little-endian; `r8`/`r4`/`lr8`/`lr4` = Fortran-unformatted-record
 * variants (big/little-endian, 8/4-byte floats, always 4-byte ints) wrapping
 * the whole body in exactly 2 Fortran records (record 1 = the 7-integer
 * header, record 2 = everything else) — each record framed by a leading and
 * trailing integer byte-count that is written but **not validated on
 * re-read**. All binary byte-swapping is done host-relative via
 * `detail/byteswap.hpp` intrinsics, never a per-byte loop.
 *
 * Body layout (fixed order): header of 7 ints (`num_points, num_triangle,
 * num_quad, num_tetra, num_pyramid, num_wedge, num_hexahedron`), then
 * points, triangle connectivity (1-based on disk, decremented on read),
 * quad connectivity, triangle boundary tags, quad boundary tags, tetra
 * connectivity, pyramid connectivity (**permuted** `[1,0,3,4,2]` on read /
 * `[1,0,4,2,3]` on write — the one place in this format where a node-order
 * mistake would silently invert cell volumes rather than error, hence a
 * dedicated signed-volume regression test), wedge connectivity, hexahedron
 * connectivity. Volume cell types (tetra/pyramid/wedge/hexahedron) get
 * zero-filled boundary tags synthesized for uniformity with the surface
 * tags, since UGRID has no native per-volume-element tag concept.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Write a mesh to a UGRID file, flavour taken from `path`'s
 *        penultimate suffix.
 *
 * Enforces **at most one** cell block per known type (triangle, quad,
 * tetra, pyramid, wedge, hexahedron) — throws otherwise; unknown types are
 * skipped with a warning. Boundary tags come from the first integer-typed
 * `cell_data` array found (warns if more than one candidate exists),
 * defaulting to all-`1` if none is present; volume elements always get
 * all-zero tags. Pyramid connectivity is permuted `[1,0,4,2,3]` before
 * writing.
 *
 * @param path filesystem path to write; its penultimate suffix (e.g. `lb8`
 *        in `out.lb8.ugrid`, or none for ascii) selects the on-disk flavour
 * @param mesh the mesh to write
 * @throws WriteError if the mesh has more than one cell block of the same
 *         known type, or the output file cannot be opened
 * @note cell_data key consumed: the first integer-typed array (used as
 *       `"ugrid:ref"` boundary tags on the surface blocks).
 */
void write_ugrid(const std::string& rPath, const Mesh& rMesh);

/**
 * @brief Read a UGRID file, flavour taken from `path`'s penultimate suffix.
 *
 * Reads the 7-integer header then points/connectivity/tags in the fixed
 * body order described above, decoding Fortran-record framing when the
 * flavour requires it and byte-swapping when the flavour's endianness
 * differs from host order (single-instruction intrinsics, never per-byte).
 * Pyramid connectivity is permuted `[1,0,3,4,2]` after reading; 1-based
 * on-disk connectivity is decremented to meshio++'s 0-based convention.
 *
 * @param path filesystem path to read
 * @return the read Mesh, with cell blocks in the fixed type order
 *         (triangle, quad, tetra, pyramid, wedge, hexahedron) for whichever
 *         counts are non-zero in the header
 * @throws ReadError on a truncated file or malformed Fortran-record framing
 * @note cell_data key produced: `"ugrid:ref"`, one array per written cell
 *       block (real boundary tags for triangle/quad, all-zero for the
 *       volume types).
 */
Mesh read_ugrid(const std::string& rPath);

}  // namespace meshioplusplus
