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
 * @file med.hpp
 * @brief MED/Salome (.med) HDF5-backed C++ reader/writer.
 *
 * MED (Salome/Code-Aster) is the most structurally involved format
 * meshio++ supports. On disk it is HDF5 with groups written in order:
 * `INFOS_GENERALES` (MAJ/MIN/REL version triple), `ENS_MAA/<mesh_name>`
 * (mesh-level attrs `DIM`/`ESP`/`UNT`/`UNI`/`NOM`/`DES`) holding a single
 * time-step group with `NOE` (nodes: `COO` Fortran-order-flattened
 * coordinates, optional `FAM` per-point family id) and `MAI` (one group per
 * cell block keyed by MED type, e.g. `HE8`, each with a Fortran-order
 * 1-based `NOD` connectivity and optional `FAM`), plus `FAS/<mesh_name>`
 * (family/group definitions: `FAMILLE_ZERO`, `NOEUD/FAM_<id>_.../GRO/NOM`,
 * `ELEME/...` with the same layout) and `CHA/<field_name>` (fields).
 * Coordinate/connectivity arrays are stored **Fortran-ordered**
 * (column-major); this C++ implementation flattens/unflattens explicitly
 * since C++ has no native Fortran-order array type (`med.cpp`'s transpose
 * uses Eigen when `MESHIOPLUSPLUS_HAS_EIGEN`, else a hand-written
 * transpose).
 *
 * **What the C++ path handles** (matching the Python output byte-for-byte):
 * points, point/cell tags, families with `GRO` group names, mesh-level
 * metadata (`mesh_name`/`description`/`unit_time`/`unit_coords`/
 * `point_tag_groups`/`cell_tag_groups`, all carried via #MedInfo), the
 * fixed node-orientation permutations for linear 3D types (`tetra`,
 * `pyramid`, `wedge`, `hexahedron` — see `_med_node_perm` in
 * doc/formats/med.md), and ragged `POG`/`POG2` polygon cell blocks (CSR
 * `NOD`+`INN` offset arrays, copied — not zero-copy — across the C++/Python
 * boundary since ragged data cannot be viewed in place). The `MAI` cell
 * blocks are iterated in HDF5 **creation order** (matching h5py's
 * `track_order`) since block order must align with `cell_data`/`cell_sets`.
 *
 * **What always falls back to Python** (the C++ functions `throw` and the
 * `meshioplusplus.med` shim catches and retries with the pure-Python/h5py
 * implementation): any file/mesh carrying `CHA` **fields** (MED-4.1
 * bitmask attributes, `field_data["med:field_units"]`/`["med:step_meta"]`,
 * and multi-timestep field-name grouping are Python-only), the
 * `gmsh:physical`→family **bridging** performed on write, non-default
 * **profiles** / `ELGA` support, and **multi-mesh** files
 * (`read_med_multi`/`write_med_multi`, which have no C++ equivalent at
 * all). Quadratic 3D types (`tetra10`, `hexahedron20`, `pyramid13`,
 * `wedge15`) share the linear types' orientation convention but have no
 * implemented corners+midpoints permutation yet — they round-trip
 * unconverted (a warning is logged the first time one is seen).
 */

#ifdef MESHIOPLUSPLUS_HAS_HDF5

// System includes
#include <cstdint>
#include <map>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Side-channel struct carrying MED mesh-level metadata that the
 *        zero-copy Mesh conversion layer cannot carry (custom attributes
 *        on the Python `meshio.Mesh`, not `point_data`/`cell_data`/
 *        `field_data`). The binding layer `setattr`s these fields onto the
 *        Python `Mesh` on read, and reads them back off it on write.
 */
struct MedInfo {
    /**
     * Per-point family/tag membership: `set_id -> [subset_name, ...]`,
     * corresponding to Python `mesh.point_tags`. Populated on read from
     * `FAS/NOEUD` family group names; consumed on write to build the
     * `NOEUD` family definitions (and each point's `FAM` id, which is
     * itself carried in `point_data["point_tags"]`, not here).
     */
    std::map<std::int64_t, std::vector<std::string>> mPointTags;
    /**
     * Per-cell-block family/tag membership: `set_id -> [subset_name, ...]`,
     * corresponding to Python `mesh.cell_tags`. Populated on read from
     * `FAS/ELEME` family group names (including the synthetic families
     * created for OpenFOAM boundary patches or Gmsh physical groups when
     * bridged elsewhere); consumed on write analogously to `point_tags`.
     */
    std::map<std::int64_t, std::vector<std::string>> mCellTags;
    /** `field_data["med:nom"]` — one component-name list per field, in
     * field-iteration order (point_data fields first, then cell_data
     * fields); each field's `NOM` attribute is the 16-char-padded
     * concatenation of its entry here. */
    std::vector<std::vector<std::string>> mMedNom;  // field_data["med:nom"]

    // Mesh-level metadata attributes (custom attributes on the Python Mesh).
    /** Python `mesh.mesh_name` — the `ENS_MAA` group name and the mesh's
     * `NOM` attribute value; defaults to `"mesh"` when absent. */
    std::string mMeshName = "mesh";
    /** Python `mesh.description` — the `DES` attribute; defaults to
     * `"Mesh created with meshio++"` on write when unset. */
    std::string mDescription;
    /** Python `mesh.unit_time` — the `UNT` attribute (physical unit of the
     * time axis, e.g. `"s"`). */
    std::string mUnitTime;
    /** Python `mesh.unit_coords` — the `UNI` attribute (physical unit of
     * the coordinate axes, e.g. `"m"`). Values round-trip through
     * `latin-1` and are stripped of surrounding whitespace/NUL padding on
     * read. */
    std::string mUnitCoords;
    // set_id -> family link name (e.g. "FAM_2_Side").
    /** `set_id -> "FAM_<id>..."` short family link name, mirroring Python
     * `mesh.point_tag_groups`; always present (possibly empty) after any
     * Python `read()` regardless of whether the source file had a `FAS`
     * section. */
    std::map<std::int64_t, std::string> mPointTagGroups;
    /** `set_id -> "FAM_<id>..."` short family link name for cell-block
     * families, mirroring Python `mesh.cell_tag_groups`; same defaulting
     * behavior as `point_tag_groups`. */
    std::map<std::int64_t, std::string> mCellTagGroups;
};

/**
 * @brief Read a MED (.med) HDF5 file into a Mesh, handling the
 *        mesh-representation subset described in the file-level docs.
 *
 * Reads points (un-transposing the Fortran-ordered `COO` dataset), the
 * `MAI` cell blocks in HDF5 creation order, per-point/per-cell `FAM` tag
 * arrays (exposed as `point_data["point_tags"]`/`cell_data["cell_tags"]`),
 * family/group names from `FAS` (searched first under the mesh's own
 * time-step group, then at the top level), mesh-level metadata, the fixed
 * node-orientation permutation for linear 3D types, and ragged
 * `POG`/`POG2` polygon blocks (materialized as a copied `list`-like ragged
 * `CellBlock` since they cannot be represented as a rectangular NDArray
 * without loss).
 *
 * @param rPath filesystem path to the .med file to read
 * @param rInfo output side-channel struct populated with tags, families,
 *        and mesh-level metadata (see #MedInfo)
 * @return the read Mesh (points, cells, point_data["point_tags"],
 *         cell_data["cell_tags"], arbitrary named point/cell data from
 *         `CHA` fields except those excluded below)
 * @throws ReadError — on any `CHA` field, non-default profile, `ELGA`
 *         support, or multi-mesh file; on malformed/unsupported HDF5
 *         layout. Callers (the Python shim) catch this and retry with the
 *         pure-Python/h5py reader.
 */
Mesh read_med(const std::string& rPath, MedInfo& rInfo);

/**
 * @brief Write a Mesh to a MED (.med) HDF5 file, handling the
 *        mesh-representation subset described in the file-level docs.
 *
 * Writes `INFOS_GENERALES` (parsed from `med_version`, falling back to
 * `4, 1, 0` if unparsable), `ENS_MAA/<rInfo.mMeshName>` with points
 * (Fortran-order-flattened) and one `MAI/<MED type>` group per cell block
 * (rejecting up front with `WriteError` if two blocks share a MED type,
 * since MED cannot represent that), `FAS` family definitions built from
 * `rInfo.mPointTags`/`rInfo.mCellTags` (a family with no groups omits `GRO`
 * entirely), and the fixed node-orientation permutation applied to linear
 * 3D cell types before writing `NOD`. Ragged `polygon`/`polygon2` blocks
 * are written as `POG`/`POG2` CSR data. Family names longer than 80 bytes
 * after `latin-1` encoding raise `WriteError` rather than truncating.
 *
 * @param rPath filesystem path to the .med file to create/overwrite
 * @param rMesh the mesh to write
 * @param rInfo side-channel struct supplying tags, families, and mesh-level
 *        metadata (see #MedInfo); read from the Python Mesh's custom
 *        attributes by the binding layer
 * @param rMedVersion the `MAJ.MIN.REL` triple written to
 *        `INFOS_GENERALES` (default `"4.1.0"`)
 * @throws WriteError — if the mesh carries `CHA`-worthy fields (any
 *         point_data/cell_data beyond `point_tags`/`cell_tags` that this
 *         path doesn't handle), `gmsh:physical` bridging is needed, two
 *         cell blocks share one MED type, or a family name exceeds 80
 *         bytes. Callers (the Python shim) catch this and retry with the
 *         pure-Python/h5py writer.
 * @note point_data/cell_data keys produced/consumed: `"point_tags"`,
 *       `"cell_tags"`.
 */
void write_med(const std::string& rPath, const Mesh& rMesh, const MedInfo& rInfo,
               const std::string& rMedVersion = "4.1.0");

}  // namespace meshioplusplus

#endif  // MESHIOPLUSPLUS_HAS_HDF5
