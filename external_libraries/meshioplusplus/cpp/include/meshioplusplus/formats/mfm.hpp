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
 * @file mfm.hpp
 * @brief Modulef Formatted Mesh (.mfm) ASCII C++ reader/writer.
 *
 * MFM (used by FEconv, a simplified NOPO/Modulef mesh) holds a **single
 * non-hybrid element type** per file: an 8-integer header
 * `nel nnod nver dim lnn lnv lne lnf`, then a flat whitespace-token stream
 * (no line-structure requirement) laid out as connectivity `mm`
 * (`nel × lnv`, 1-based, element-major), a face reference array `nrc`
 * (`nel × lnf`, present only if `dim == 3`), an edge reference array `nra`
 * (`nel × lne`, present only if `dim >= 2`), a vertex reference array `nrv`
 * (`nel × lnv`, always present), vertex coordinates `z` (`nver × dim`,
 * vertex-major), and a per-element subdomain array `nsd` (`nel` values ->
 * `cell_data["mfm:ref"]`). `nrc`/`nra`/`nrv` are read-and-discarded (no
 * meshio++-side representation) and always written back as zeros.
 *
 * The element type is recovered from `(lnv, lne, lnf)` plus `lnn == lnv`
 * (`line`, `triangle`, `quad`, `tetra`, `hexahedron`, `wedge` — linear only;
 * `lnn != lnv` or `nnod != nver` would imply a second-order element MFM
 * cannot store without losing curvature, so those are rejected outright).
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Write a Mesh to an MFM (.mfm) file.
 *
 * Requires every cell in the mesh to share exactly one linear type
 * (`line`, `triangle`, `quad`, `tetra`, `hexahedron`, or `wedge`); a
 * mixed-type mesh raises `WriteError` since MFM is fundamentally
 * single-type ("non-hybrid"). Emits the 8-int header, 1-based
 * element-major connectivity, all-zero `nrc`/`nra`/`nrv` placeholders,
 * vertex-major coordinates formatted with `float_fmt`, and the per-element
 * subdomain array from `cell_data["mfm:ref"]` (defaulting to all-ones if
 * absent).
 *
 * @param rPath filesystem path to the .mfm file to create/overwrite
 * @param rMesh the mesh to write (must be single-cell-type, linear)
 * @param rFloatFmt coordinate format string (e.g. `".16e"`)
 * @throws WriteError if the mesh has more than one cell type, a
 *         higher-order cell type, or is otherwise unsupported
 * @note reads `cell_data["mfm:ref"]` if present
 */
void write_mfm(const std::string& rPath, const Mesh& rMesh, const std::string& rFloatFmt);

/**
 * @brief Read an MFM (.mfm) file into a Mesh.
 *
 * Parses the 8-int header to recover `nel`/`nver`/`dim` and the element
 * type from `(lnv, lne, lnf)`/`lnn`, then reads connectivity (shifted from
 * 1-based to 0-based), skips `nrc` (if `dim == 3`)/`nra` (if `dim >= 2`)/
 * `nrv` unconditionally-present sections without storing them, reads
 * vertex coordinates, and reads the per-element subdomain array into
 * `cell_data["mfm:ref"]`.
 *
 * @param rPath filesystem path to the .mfm file to read
 * @return the read Mesh, single cell block, with `cell_data["mfm:ref"]`
 *         populated from the file's `nsd` array
 * @throws ReadError if `lnn != lnv` (would imply a second-order element),
 *         `nnod != nver`, or the `(lnv, lne, lnf)` triple doesn't match a
 *         known linear type
 */
Mesh read_mfm(const std::string& rPath);

}  // namespace meshioplusplus
