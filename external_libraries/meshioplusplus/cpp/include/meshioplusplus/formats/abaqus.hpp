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
 * @file abaqus.hpp
 * @brief Abaqus input-deck (.inp) C++ reader/writer.
 *
 * The Abaqus format is a keyword-driven ASCII deck: comment lines start
 * `**`; keyword lines start `*` and are matched on
 * `line.partition(",")[0].strip().replace("*","").upper()`. This
 * implementation covers `*NODE` (comma-separated `id, x, y, [z]` rows) and
 * `*ELEMENT, TYPE=<abaqus_type>[, ELSET=<name>]` (comma-separated integer
 * rows, flattened and split into fixed-width records of
 * `node_count(TYPE) + 1`). Any other recognized keyword —`*NSET`, `*ELSET`,
 * `*INCLUDE` — is refused outright by the C++ reader (see @ref read_abaqus),
 * deferring the whole file to the Python fallback, since `point_sets`/
 * `cell_sets` (built from those keywords) are not carried by the Mesh
 * conversion layer.
 *
 * Cell types go through the Abaqus <-> meshio++ element-name table (trusses,
 * beams, shells, solids -> `line`/`line3`/`triangle`/`triangle6`/`quad`/
 * `quad8`/`quad9`/`tetra`/`tetra10`/`hexahedron`/`hexahedron20`/`wedge`/
 * `wedge15`, plus the asymmetric `C3D4H` -> `"tetra4"` entry); see
 * doc/formats/abaqus.md for the full table and its "known table quirk" note.
 * The reverse (meshio++ -> Abaqus) map is lossy: several Abaqus names
 * collapse onto one meshio++ type, so the writer always emits whichever name
 * is last in the internal table for that type.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Write `mesh` as an Abaqus .inp file (`*NODE`/`*ELEMENT` only).
 *
 * Emits one `*NODE` block (1-based ids matching row position) followed by
 * one `*ELEMENT, TYPE=<abaqus_type>` block per cell block, translating each
 * meshio++ cell type through the Abaqus element-name table.
 *
 * @param rPath filesystem path to write
 * @param rMesh the mesh to write
 * @throws WriteError if a cell block's type has no Abaqus element-name
 *         mapping (`"Abaqus writer: unsupported cell type ..."`)
 * @note the shim only attempts this C++ path when `float_fmt == ".16e"`,
 *       `translate_cell_names == True`, and the mesh has no `point_sets`/
 *       `cell_sets` — anything else falls back to the Python writer, which
 *       also supports `translate_cell_names=False` (verbatim type strings).
 */
void write_abaqus(const std::string& rPath, const Mesh& rMesh);

/**
 * @brief Read an Abaqus .inp file (`*NODE`/`*ELEMENT` only).
 *
 * Parses `*NODE` rows into points (keyed by the file's own, possibly
 * non-contiguous, node ids) and `*ELEMENT, TYPE=...` rows into cell blocks,
 * looking up each Abaqus type name first upper-cased then, if that misses,
 * case-sensitively as written (a leniency the plain-dict Python reader does
 * not have).
 *
 * @param rPath filesystem path to read
 * @return the read Mesh (no point_data/cell_data/field_data — this reader
 *         never populates them)
 * @throws ReadError if the file can't be opened, an `*ELEMENT` card has no
 *         `TYPE=`, the type isn't in the lookup table, the node count for a
 *         type is unknown, a data row has the wrong stride, a referenced
 *         node id is unknown, or the file uses `*NSET`/`*ELSET`/`*INCLUDE`
 *         (always deferred to the Python fallback, which supports them,
 *         including `GENERATE` ranges and recursive `*ELSET` references)
 */
Mesh read_abaqus(const std::string& rPath);

}  // namespace meshioplusplus
