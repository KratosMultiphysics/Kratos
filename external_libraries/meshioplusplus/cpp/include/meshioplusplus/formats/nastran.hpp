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
 * @file nastran.hpp
 * @brief MSC/NX Nastran bulk-data (.bdf/.fem/.nas) C++ writer + reader.
 *
 * Nastran bulk data is fixed-width card text (`GRID`, `CTRIA3`, `CTETRA`,
 * `CHEXA`, …) between a `"BEGIN BULK"` line and `"ENDDATA"`, in
 * small-field (10 x 8-char fields), large-field (8 + 4x16 + 8 chars, `*`
 * continuation), or free (comma-separated) layout. This C++ implementation
 * only emits/parses the `fixed-large`/`fixed-small` point/cell-format
 * combination.
 *
 * **The reader is sentinel-gated**: it only accepts files whose first `$`
 * comment line is the exact literal string the C++ writer itself emits
 * (`"meshioplusplus-cpp-nastran"`). Any real-world Nastran file — including
 * this project's own reference `.fem` fixtures — lacks that sentinel and
 * is therefore always parsed by the more permissive Python reader instead;
 * this is the single most consequential interop rule for this format (see
 * doc/formats/nastran.md). The shim likewise only attempts the C++ writer
 * for the exact `fixed-large`/`fixed-small` combination and only when no
 * `nastran:ref` data is present.
 *
 * Cell-type map includes `CTETRA`/`CPYRAM`/`CPENTA`/`CHEXA` auto-upgraded
 * to their 10/13/15/20-node quadratic meshio++ counterparts whenever a card
 * lists more node ids than the linear element's base count (a heuristic,
 * not a version flag). Node-order permutations are applied for
 * `triangle6`/`CTRAX6`/`CTRIAX6` (to-VTK `[0,2,4,1,3,5]`, to-Nastran
 * `[0,3,1,4,2,5]`), `hexahedron20`
 * (`[0,1,2,3,4,5,6,7,8,9,10,11,16,17,18,19,12,13,14,15]`), and `wedge15`
 * (`[0,1,2,3,4,5,6,7,8,12,13,14,9,10,11]`).
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Write a Mesh to a Nastran bulk-data file (fixed-large/fixed-small
 *        layout only).
 *
 * Emits `GRID*` large-field point cards (16-character floats found by
 * searching increasing precision, 0 through 11, for the shortest string
 * that round-trips exactly via `strtod`, then trimming trailing mantissa
 * zeros — not guaranteed byte-identical to the Python writer's
 * `np.format_float_scientific(precision=11)` + `e`->`E` approach, but
 * targeting the same 16-char field) and fixed-small element cards, plus
 * the `"meshioplusplus-cpp-nastran"` sentinel comment as the first `$`
 * line (required for this writer's own output to be read back by
 * #read_nastran). 2D points are force-promoted to 3D with a warning.
 *
 * @param rPath filesystem path to the .bdf/.fem/.nas file to create/overwrite
 * @param rMesh the mesh to write
 * @throws WriteError on an unsupported cell type
 */
void write_nastran(const std::string& rPath, const Mesh& rMesh);

/**
 * @brief Read a Nastran bulk-data file into a Mesh — only accepts files
 *        carrying this writer's sentinel comment.
 *
 * Parses `GRID`/`GRID*` and element cards between `"BEGIN BULK"` and
 * `"ENDDATA"`, decoding Nastran's compressed-exponent float notation
 * (e.g. `1.5+1`) and re-merging large-field continuation lines. Applies
 * the inverse node-order permutation for `triangle6`, `hexahedron20`, and
 * `wedge15`. `CBAR`/`CBEAM`/`CBUSH`/`CBUSH1D`/`CGAP` cards only keep their
 * first 2 node ids (a 3rd orientation/grid-id field is discarded).
 *
 * @param rPath filesystem path to the .bdf/.fem/.nas file to read
 * @return the read Mesh
 * @throws ReadError if the file's first `$` comment line is not exactly
 *         `"meshioplusplus-cpp-nastran"` (routes real-world Nastran files
 *         to the Python fallback), or on a malformed card
 * @note point_data/cell_data keys are not produced by this reader (unlike
 *       the Python reference, which populates `"nastran:ref"` from the
 *       optional GRID/element reference field)
 */
Mesh read_nastran(const std::string& rPath);

}  // namespace meshioplusplus
