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
 * @file mff.hpp
 * @brief Modulef Formatted Field (.mff) C++ reader/writer.
 *
 * An MFF file is the field companion to a Modulef mesh (.mfm): an integer
 * value count followed by a flat list of double-precision floats, with no
 * geometry or component/location metadata. Read here as a geometry-less Mesh
 * (no cells, `points` with zero columns) carrying `point_data["mff:field"]`;
 * written from the first `point_data` array (or first non-`unv:pid`
 * `cell_data` array). Standalone, only field values round-trip.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/** @brief Read a Modulef Formatted Field (.mff) into a geometry-less Mesh. */
Mesh read_mff(const std::string& rPath);

/** @brief Write a mesh's first field as a Modulef Formatted Field (.mff). */
void write_mff(const std::string& rPath, const Mesh& rMesh);

}  // namespace meshioplusplus
