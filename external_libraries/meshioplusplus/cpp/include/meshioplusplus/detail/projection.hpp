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
 * @file projection.hpp
 * @brief Orthographic camera projection + painter's-algorithm face ordering
 * for the SVG/TikZ 3D rendering path.
 *
 * The camera is parameterized by azimuth/elevation/roll in degrees: the view
 * direction (from the scene toward the camera) is
 * `w = (cos el * cos az, cos el * sin az, sin el)`, the screen-right axis is
 * `u = normalize(cross((0,0,1), w))` (falling back to the y axis as the up
 * reference when looking straight along z), the screen-up axis is
 * `v = cross(w, u)`, and `roll` rotates `(u, v)` about `w`. The defaults
 * `azimuth = 45`, `elevation = 35.264389682754654` (`atan(1/sqrt(2))`) give
 * `w` proportional to `(1,1,1)` — the classic CAD isometric view. Projection
 * is orthographic (`x' = p.u`, `y' = p.v`), and `depth = p.w` orders faces
 * back-to-front (painter's algorithm; larger depth = closer to the camera).
 *
 * KEEP IN SYNC: `src/meshioplusplus/_projection.py` is the Python twin used
 * by the pure-Python fallback writers — the arithmetic (down to expression
 * order, which fixes the floating-point rounding) must stay identical so
 * TikZ output is byte-identical across the two implementations.
 */

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {
namespace detail {

/** @brief Orthonormal camera basis: screen-right `u`, screen-up `v`, and the
 * view direction `w` (scene toward camera). */
struct CameraBasis {
    std::array<double, 3> mU;
    std::array<double, 3> mV;
    std::array<double, 3> mW;
};

/**
 * @brief Build the orthographic camera basis for azimuth/elevation/roll.
 * @param azimuth Camera azimuth in degrees (rotation about z, from +x).
 * @param elevation Camera elevation in degrees above the xy plane.
 * @param roll In-screen rotation in degrees about the view direction.
 * @return The `{u, v, w}` basis described in the file header.
 */
inline CameraBasis camera_basis(double azimuth, double elevation, double roll) {
    constexpr double deg = 3.141592653589793 / 180.0;
    const double az = azimuth * deg;
    const double el = elevation * deg;
    const double caz = std::cos(az);
    const double saz = std::sin(az);
    const double cel = std::cos(el);
    const double sel = std::sin(el);
    const double w0 = cel * caz;
    const double w1 = cel * saz;
    const double w2 = sel;

    double u0, u1, u2;
    if (std::fabs(w2) > 1.0 - 1.0e-12) {
        // Looking straight along z: use y as the up reference.
        // u = normalize(cross((0,1,0), w))
        const double n = std::sqrt(w2 * w2 + w0 * w0);
        u0 = w2 / n;
        u1 = 0.0;
        u2 = -w0 / n;
    } else {
        // u = normalize(cross((0,0,1), w))
        const double n = std::sqrt(w1 * w1 + w0 * w0);
        u0 = -w1 / n;
        u1 = w0 / n;
        u2 = 0.0;
    }
    // v = cross(w, u)
    double v0 = w1 * u2 - w2 * u1;
    double v1 = w2 * u0 - w0 * u2;
    double v2 = w0 * u1 - w1 * u0;

    if (roll != 0.0) {
        const double cr = std::cos(roll * deg);
        const double sr = std::sin(roll * deg);
        const double ru0 = cr * u0 + sr * v0;
        const double ru1 = cr * u1 + sr * v1;
        const double ru2 = cr * u2 + sr * v2;
        const double rv0 = -sr * u0 + cr * v0;
        const double rv1 = -sr * u1 + cr * v1;
        const double rv2 = -sr * u2 + cr * v2;
        u0 = ru0;
        u1 = ru1;
        u2 = ru2;
        v0 = rv0;
        v1 = rv1;
        v2 = rv2;
    }
    return CameraBasis{{u0, u1, u2}, {v0, v1, v2}, {w0, w1, w2}};
}

/** @brief One drawable face of a projected surface: 2-4 corner point ids. */
struct ProjectedFace {
    std::array<std::int64_t, 4> mNodes;
    std::uint8_t mNumNodes;  // 2 (line), 3 (triangle), or 4 (quad)
    bool mIsLine;
    double mDepth;  // view-space depth of the face centroid
};

/** @brief A surface mesh projected to screen space, faces sorted
 * back-to-front. */
struct ProjectedSurface {
    std::vector<double> mX;
    std::vector<double> mY;
    std::vector<ProjectedFace> mFaces;
};

/**
 * @brief Project a surface mesh's points with an orthographic camera and
 * gather its drawable faces sorted back-to-front (painter's algorithm).
 *
 * Drawable blocks are `line`, `triangle`, `quad`, and the corner-linearized
 * higher-order surface types (`triangle6`, `quad8`, `quad9` — corner nodes
 * only); anything else is skipped. Faces keep enumeration order within
 * equal depths (`std::stable_sort`, matching numpy's stable argsort).
 *
 * @param rMesh The surface mesh to project (already skin-extracted when it
 *              came from a volume mesh).
 * @param azimuth Camera azimuth in degrees.
 * @param elevation Camera elevation in degrees.
 * @param roll In-screen rotation in degrees.
 * @return Projected screen coordinates per point and the sorted face list.
 */
inline ProjectedSurface project_surface(const Mesh& rMesh, double azimuth, double elevation,
                                        double roll) {
    const CameraBasis cam = camera_basis(azimuth, elevation, roll);
    const NDArray& points = rMesh.Points();
    const std::size_t num_points = rMesh.NumPoints();
    const std::size_t dim = rMesh.PointDim();

    ProjectedSurface out;
    out.mX.resize(num_points);
    out.mY.resize(num_points);
    std::vector<double> depth(num_points);
    for (std::size_t i = 0; i < num_points; ++i) {
        const double px = (0 < dim) ? read_double(points, i * dim + 0) : 0.0;
        const double py = (1 < dim) ? read_double(points, i * dim + 1) : 0.0;
        const double pz = (2 < dim) ? read_double(points, i * dim + 2) : 0.0;
        out.mX[i] = px * cam.mU[0] + py * cam.mU[1] + pz * cam.mU[2];
        out.mY[i] = px * cam.mV[0] + py * cam.mV[1] + pz * cam.mV[2];
        depth[i] = px * cam.mW[0] + py * cam.mW[1] + pz * cam.mW[2];
    }

    for (const auto cb : rMesh.CellRange()) {
        const std::string& type = cb.Type();
        std::uint8_t n_corner = 0;
        bool is_line = false;
        if (type == "line") {
            n_corner = 2;
            is_line = true;
        } else if (type == "triangle" || type == "triangle6") {
            n_corner = 3;
        } else if (type == "quad" || type == "quad8" || type == "quad9") {
            n_corner = 4;
        } else {
            continue;
        }
        const NDArray& conn = cb.Conn();
        const std::size_t ncols = cols(conn);
        const std::size_t n = cb.NumCells();
        for (std::size_t r = 0; r < n; ++r) {
            ProjectedFace face;
            face.mNodes = {-1, -1, -1, -1};
            face.mNumNodes = n_corner;
            face.mIsLine = is_line;
            double d = 0.0;
            for (std::uint8_t k = 0; k < n_corner; ++k) {
                const std::int64_t p = read_int(conn, r * ncols + k);
                face.mNodes[k] = p;
                d += depth[static_cast<std::size_t>(p)];
            }
            face.mDepth = d / static_cast<double>(n_corner);
            out.mFaces.push_back(face);
        }
    }

    std::stable_sort(
        out.mFaces.begin(), out.mFaces.end(),
        [](const ProjectedFace& rA, const ProjectedFace& rB) { return rA.mDepth < rB.mDepth; });
    return out;
}

}  // namespace detail
}  // namespace meshioplusplus
