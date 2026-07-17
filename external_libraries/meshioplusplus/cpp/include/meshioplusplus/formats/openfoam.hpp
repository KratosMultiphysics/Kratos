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
 * @file openfoam.hpp
 * @brief OpenFOAM polyMesh (read-only) C++ reader.
 *
 * A polyMesh is a directory of sibling "FoamFile"-headered files (`points`,
 * `faces`, `owner`, `neighbour`, `boundary`), each ASCII or binary
 * (little-endian only; `label=32/64`, `scalar=32/64` per the file's `arch`
 * header string). `points` (`vectorField`) and `owner`/`neighbour`
 * (`labelList`) are flat contiguous buffers read directly; `faces`
 * (`faceList`) is non-contiguous (each face is its own length-prefixed
 * `labelList`) and is read via a two-pass CSR gather bounded in peak
 * memory. `boundary` is a `patch_name -> {type, nFaces, startFace}` table
 * parsed with a brace-matching regex.
 *
 * Cells are reconstructed from the owner/neighbour/face topology: each
 * face is oriented outward from its owning cell (reversed if the cell is
 * that face's neighbour), then classified by `(n_faces, n_points)` into
 * `tetra` (4,4), `pyramid` (5,5), `wedge` (5,6), `hexahedron` (6,8) — each
 * with a dedicated orientation-fixing builder that flips node order if a
 * scalar triple product comes out negative — or, for any other signature,
 * a general `polyhedron<N>` **ragged** cell block (ragged data crosses the
 * C++/Python boundary as a copied list of face-node arrays, never
 * zero-copy). Boundary faces become `triangle`/`quad`/`polygon<N>` blocks,
 * one per patch/size combination.
 *
 * This reader is **read-only** — there is no OpenFOAM writer at all, in
 * C++ or Python. Only mesh topology is read; OpenFOAM field files (`U`,
 * `p`, `T`, …) under a case's time directories are never read by this
 * module, so no `point_data`/`field_data` is ever produced.
 */

// System includes
#include <cstdint>
#include <map>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Side-channel struct carrying OpenFOAM boundary-patch tag data
 *        that the zero-copy Mesh conversion layer cannot carry (`Python
 *        mesh.cell_tags` is a custom Mesh attribute, not `cell_data`). The
 *        binding layer `setattr`s this onto the returned Python `Mesh`.
 */
struct OpenFoamInfo {
    // MED-style negative family id -> {patch name}.
    /**
     * `family_id -> [patch_name]`, mirroring Python `mesh.cell_tags`. Each
     * boundary patch gets a distinct negative "MED-style family id"
     * `-(patch_index+1)` (assigned once per patch and reused across
     * whichever face-size cell blocks that patch's faces fall into); the
     * matching `cell_data["cell_tags"]` array on the returned Mesh holds
     * `0` for every volume-cell block and the patch's family id for its
     * boundary-face blocks. This lets a subsequent MED write bridge patch
     * names through the same family mechanism used for Gmsh physical
     * groups (see doc/formats/med.md).
     */
    std::map<std::int64_t, std::vector<std::string>> mCellTags;
};

// `path` may be a `.foam` marker file, a case directory, or a polyMesh
// directory (resolved like the Python reader).
/**
 * @brief Read an OpenFOAM polyMesh into a Mesh.
 *
 * `path` may be a `.foam` marker file (looks for
 * `<parent>/constant/polyMesh`), a directory literally named `polyMesh`
 * (used as-is), or any other directory (checked for `constant/polyMesh`
 * then `polyMesh` as subdirectories) — resolved identically to the Python
 * reader's `_resolve_polymesh`. Reconstructs volume cells
 * (tetra/pyramid/wedge/hexahedron/general polyhedron) and boundary faces
 * (triangle/quad/polygon) from the `points`/`faces`/`owner`/`neighbour`/
 * `boundary` files, auto-detecting ASCII vs binary and label/scalar width
 * per file. Degenerate volume cells that match a named type's
 * `(n_faces, n_points)` signature but whose topology doesn't cleanly
 * resolve are silently skipped (logged as a warning count) rather than
 * demoted to a general polyhedron.
 *
 * @param rPath a `.foam` file, case directory, or polyMesh directory
 * @param rInfo output side-channel struct populated with boundary-patch
 *        family ids and names (see #OpenFoamInfo)
 * @return the read Mesh: points, volume + boundary cell blocks,
 *         `cell_data["cell_tags"]` (0 for volume blocks, a per-patch
 *         negative id for boundary blocks), `mesh.point_tags` always set
 *         to `{}` (OpenFOAM has no point-tag concept; present only for
 *         interface symmetry with the MED-derived tag convention) — no
 *         point_data or field_data
 * @throws ReadError / std::filesystem-related errors if no polyMesh
 *         directory can be resolved, or on a malformed/unsupported file;
 *         callers (the Python shim) catch this and retry with the
 *         pure-Python reader
 */
Mesh read_openfoam(const std::string& rPath, OpenFoamInfo& rInfo);

}  // namespace meshioplusplus
