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
 * @file ip.hpp
 * @brief ANSYS Fluent interpolation file (.ip) C++ reader/writer.
 *
 * An IP file stores one or more fields over a set of points: version, spatial
 * dimension, point count, component count, the component names, then a section
 * of all values for each coordinate (x, y, z) and a section of all values for
 * each field component -- in version 3 each section is wrapped in `(`/`)`.
 * Only text files (versions 2 and 3) are supported. Read here as a
 * geometry-less Mesh (no cells) with `points` from the coordinate sections and
 * one `point_data` entry per field; written as a version-3 file.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/** @brief Read an ANSYS Fluent interpolation file (.ip) into a Mesh. */
Mesh read_ip(const std::string& rPath);

/** @brief Write a mesh's nodal fields as a version-3 interpolation file (.ip). */
void write_ip(const std::string& rPath, const Mesh& rMesh);

}  // namespace meshioplusplus
