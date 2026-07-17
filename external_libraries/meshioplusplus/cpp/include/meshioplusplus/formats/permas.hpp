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
 * @file permas.hpp
 * @brief PERMAS (.post/.dato) plain-text C++ reader/writer.
 *
 * PERMAS files are `$`-delimited keyword sections in plain text (`!` starts
 * a comment). `$COOR...` introduces the node block (`<gid> <x> <y> <z...>`
 * rows, building a `gid -> running index` map); `$ELEMENT TYPE=<permas
 * type>` introduces an element block, where PERMAS uses a **trailing `!`
 * as a line-continuation marker** — node ids accumulate across lines until
 * one does *not* end in `!`, at which point the accumulated ids become one
 * completed cell (a standalone `!` separator line between blocks must
 * yield **no** cell, not an empty one — getting this wrong previously
 * caused an out-of-bounds crash during development, now fixed and tested).
 * `$NSET`/`$ESET` blocks (including `GENERATE`, using exclusive-stop
 * `np.arange`-style semantics) are parsed but never attached to the
 * returned Mesh — a currently-dead read path, kept only for parity with
 * the Python reference. All other keywords are silently ignored.
 *
 * PERMAS <-> meshio++ type map (several PERMAS names collapse onto one
 * meshio++ type on read; the C++ writer hardcodes the same canonical
 * PERMAS name per meshio++ type that Python's dict-insertion-order
 * "last wins" reverse map produces): `PLOT1`->`vertex`, beam/rod names
 * (`FSCPIPE2` + 10 others)->`line`, `PLOTL3`->`line3`, `TRIMS3` + 6
 * others->`triangle`, `TRIMS6`->`triangle6`, `SHELL4` + 5 others->`quad`,
 * `QUAMS8`->`quad8`, `QUAMS9`->`quad9`, `HEXFO8`->`hexahedron`,
 * `HEXE20`->`hexahedron20`, `HEXE27`->`hexahedron27`, `TET4`->`tetra`,
 * `TET10`->`tetra10`, `PYRA5`->`pyramid`, `PENTA6`->`wedge`,
 * `PENTA15`->`wedge15`.
 *
 * PERMAS produces no `point_data`/`cell_data`/`field_data` at all.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Write a Mesh to a PERMAS (.post/.dato) file.
 *
 * Emits `!PERMAS DataFile Version 18.0`, a `!written by meshio++ (C++
 * core)` credit line, `$ENTER COMPONENT NAME=DFLT_COMP`, `$STRUCTURE`,
 * `$COOR` with sequential 1-based node indices (the original PERMAS gid,
 * if any, is not tracked and not reproduced), then one `$ELEMENT
 * TYPE=<permas type>` block per cell type with a continuously-incrementing
 * element id across all blocks and 1-based connectivity. **Write-only**
 * node-order permutations are applied for second-order types (no inverse
 * exists on read — see the file-level quirk): `triangle6`
 * `[0,3,1,4,2,5]`, `tetra10` `[0,4,1,5,2,6,7,8,9,3]`, `quad9`
 * `[0,4,1,7,8,5,3,6,2]`, `wedge15`
 * `[0,6,1,7,2,8,9,10,11,3,12,4,13,5,14]`. Ends with `$END STRUCTURE` /
 * `$EXIT COMPONENT` / `$FIN`.
 *
 * @param rPath filesystem path to the .post/.dato file to create/overwrite
 * @param rMesh the mesh to write
 * @throws WriteError on an unsupported cell type
 */
void write_permas(const std::string& rPath, const Mesh& rMesh);

/**
 * @brief Read a PERMAS (.post/.dato) file into a Mesh.
 *
 * Parses `$COOR` node rows (building a gid -> index map) and `$ELEMENT
 * TYPE=...` blocks, resolving node ids through that map and handling the
 * trailing-`!` line-continuation convention. `$NSET`/`$ESET` blocks are
 * parsed (including `GENERATE` expansion) but discarded — they are never
 * attached to the returned Mesh. No node-order permutation is applied on
 * read (the write-side second-order reorders have no read-side inverse:
 * see the **asymmetric quadratic round-trip** quirk — a file written by
 * this writer and read back by this reader does not restore the original
 * node order for `triangle6`/`tetra10`/`quad9`/`wedge15` without external
 * correction).
 *
 * @param rPath filesystem path to the .post/.dato file to read
 * @return the read Mesh (no point_data/cell_data/field_data — PERMAS
 *         carries none)
 * @throws ReadError on a malformed file (e.g. an unrecognized element
 *         type, or a continuation line with no terminating non-`!` line)
 */
Mesh read_permas(const std::string& rPath);

}  // namespace meshioplusplus
