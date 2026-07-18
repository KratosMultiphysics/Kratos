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
 * @file svg.hpp
 * @brief SVG (Scalable Vector Graphics) mesh writer (write-only).
 *
 * Draws the mesh's `line`/`triangle`/`quad` cells as `<path>` elements in a
 * single `<svg>` document — a visualization format with no reader. For 2D or
 * flat-3D input (all z ~ 0) the classic flat path is used, byte-identical to
 * previous releases. A genuinely non-flat 3D mesh is *rendered*: if it
 * contains supported volume cells its boundary skin is extracted first
 * (`extract_skin`, linearized), then the surface faces are projected with an
 * orthographic camera (azimuth/elevation/roll, default the classic CAD
 * isometric view) and painted back-to-front (painter's algorithm on face
 * centroid depth — see `detail/projection.hpp`). The y-axis is flipped
 * (`max_y + min_y - y`) to convert the mesh/math convention (y-up) to SVG's
 * screen convention (y-down). Any non-drawable cell type is silently
 * skipped. No point_data/cell_data/field_data is emitted.
 */

// System includes
#include <optional>
#include <string>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Write a mesh's `line`/`triangle`/`quad` cells as an SVG document.
 *
 * @param rPath        filesystem path to write
 * @param rMesh        the mesh to write (only line/triangle/quad contribute)
 * @param rFloatFmt    printf-style float format for coordinates without the
 *                     leading '%' (e.g. `".3f"`)
 * @param rStrokeWidth explicit stroke width; `std::nullopt` auto-computes it as
 *                     1% of the on-canvas width
 * @param rImageWidth  output width in user units; `std::nullopt` keeps the
 *                     mesh's own width (no scaling)
 * @param rFill        cell fill colour
 * @param rStroke      edge stroke colour
 * @param azimuth      camera azimuth in degrees (3D input only)
 * @param elevation    camera elevation in degrees (3D input only); the
 *                     default pair (45, atan(1/sqrt(2))) is the classic CAD
 *                     isometric view
 * @param roll         in-screen camera roll in degrees (3D input only)
 * @throws WriteError on an unopenable output path
 */
void write_svg(const std::string& rPath, const Mesh& rMesh, const std::string& rFloatFmt = ".3f",
               const std::optional<std::string>& rStrokeWidth = std::nullopt,
               const std::optional<double>& rImageWidth = 100.0,
               const std::string& rFill = "#c8c5bd", const std::string& rStroke = "#000080",
               double azimuth = 45.0, double elevation = 35.264389682754654, double roll = 0.0);

}  // namespace meshioplusplus
