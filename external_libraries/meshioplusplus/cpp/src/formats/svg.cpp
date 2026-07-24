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

// System includes
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/detail/projection.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/svg.hpp"
#include "meshioplusplus/skin.hpp"

namespace meshioplusplus {

namespace {

// Format a single double with a printf-style spec (spec without leading '%',
// e.g. ".3f"), mirroring the Python reference's `format(x, float_fmt)`.
std::string svg_fmt_num(double value, const std::string& rSpec) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), ("%" + rSpec).c_str(), value);
    return buf;
}

// 3D rendering path: project the (already skin-extracted) surface mesh with
// the orthographic camera and emit its faces back-to-front. Replicates the
// flat path's bbox -> y-flip -> image_width scaling -> <path> emission, but
// over the sorted projected faces.
void svg_proj_write(const std::string& rPath, const Mesh& rDrawMesh, const std::string& rFloatFmt,
                    const std::optional<std::string>& rStrokeWidth,
                    const std::optional<double>& rImageWidth, const std::string& rFill,
                    const std::string& rStroke, double azimuth, double elevation, double roll) {
    detail::ProjectedSurface ps = detail::project_surface(rDrawMesh, azimuth, elevation, roll);
    std::vector<double>& x = ps.mX;
    std::vector<double>& y = ps.mY;
    const std::size_t num_points = x.size();

    double min_x = 0.0, max_x = 0.0, min_y = 0.0, max_y = 0.0;
    if (num_points > 0) {
        min_x = max_x = x[0];
        min_y = max_y = y[0];
        for (std::size_t i = 1; i < num_points; ++i) {
            min_x = std::min(min_x, x[i]);
            max_x = std::max(max_x, x[i]);
            min_y = std::min(min_y, y[i]);
            max_y = std::max(max_y, y[i]);
        }
    }

    // Flip y (projected math convention y-up -> SVG screen convention y-down).
    for (std::size_t i = 0; i < num_points; ++i)
        y[i] = max_y + min_y - y[i];

    double width = max_x - min_x;
    double height = max_y - min_y;

    if (rImageWidth.has_value() && width != 0.0) {
        const double scaling_factor = *rImageWidth / width;
        min_x *= scaling_factor;
        min_y *= scaling_factor;
        width *= scaling_factor;
        height *= scaling_factor;
        for (std::size_t i = 0; i < num_points; ++i) {
            x[i] *= scaling_factor;
            y[i] *= scaling_factor;
        }
    }

    std::string stroke_width;
    if (rStrokeWidth.has_value()) {
        stroke_width = *rStrokeWidth;
    } else {
        char buf[64];
        std::snprintf(buf, sizeof(buf), "%g", width / 100.0);
        stroke_width = buf;
    }

    std::ofstream os(rPath, std::ios::binary);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);

    os << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" viewBox=\""
       << svg_fmt_num(min_x, rFloatFmt) << ' ' << svg_fmt_num(min_y, rFloatFmt) << ' '
       << svg_fmt_num(width, rFloatFmt) << ' ' << svg_fmt_num(height, rFloatFmt) << "\">";

    os << "<style>path {fill: " << rFill << "; stroke: " << rStroke
       << "; stroke-width: " << stroke_width << "; stroke-linejoin:bevel}</style>";

    for (const detail::ProjectedFace& face : ps.mFaces) {
        std::string d;
        for (std::uint8_t k = 0; k < face.mNumNodes; ++k) {
            const std::int64_t p = face.mNodes[k];
            d += (k == 0) ? "M " : "L ";
            d += svg_fmt_num(x[static_cast<std::size_t>(p)], rFloatFmt);
            d += ' ';
            d += svg_fmt_num(y[static_cast<std::size_t>(p)], rFloatFmt);
        }
        if (!face.mIsLine)
            d += "Z";
        os << "<path d=\"" << d << "\" />";
    }

    os << "</svg>";
}

}  // namespace

void write_svg(const std::string& rPath, const Mesh& rMesh, const std::string& rFloatFmt,
               const std::optional<std::string>& rStrokeWidth,
               const std::optional<double>& rImageWidth, const std::string& rFill,
               const std::string& rStroke, double azimuth, double elevation, double roll) {
    const NDArray& points = rMesh.Points();
    const std::size_t num_points = rMesh.NumPoints();
    const std::size_t dim = rMesh.PointDim();

    // A genuinely non-flat 3D mesh takes the projected-rendering path (skin
    // extraction for volume cells + orthographic camera); a flat one (every
    // z ~ 0) keeps the classic 2D path below, byte-identical to before.
    if (dim == 3) {
        for (std::size_t i = 0; i < num_points; ++i) {
            if (std::fabs(detail::read_double(points, i * dim + 2)) > 1.0e-14) {
                if (has_skinnable_cells(rMesh)) {
                    svg_proj_write(rPath, extract_skin(rMesh, /*linearize=*/true), rFloatFmt,
                                   rStrokeWidth, rImageWidth, rFill, rStroke, azimuth, elevation,
                                   roll);
                } else {
                    svg_proj_write(rPath, rMesh, rFloatFmt, rStrokeWidth, rImageWidth, rFill,
                                   rStroke, azimuth, elevation, roll);
                }
                return;
            }
        }
    }

    // Copy the first two coordinate columns.
    std::vector<double> x(num_points), y(num_points);
    for (std::size_t i = 0; i < num_points; ++i) {
        x[i] = (0 < dim) ? detail::read_double(points, i * dim + 0) : 0.0;
        y[i] = (1 < dim) ? detail::read_double(points, i * dim + 1) : 0.0;
    }

    double min_x = 0.0, max_x = 0.0, min_y = 0.0, max_y = 0.0;
    if (num_points > 0) {
        min_x = max_x = x[0];
        min_y = max_y = y[0];
        for (std::size_t i = 1; i < num_points; ++i) {
            min_x = std::min(min_x, x[i]);
            max_x = std::max(max_x, x[i]);
            min_y = std::min(min_y, y[i]);
            max_y = std::max(max_y, y[i]);
        }
    }

    // Flip y (mesh math convention y-up -> SVG screen convention y-down).
    for (std::size_t i = 0; i < num_points; ++i)
        y[i] = max_y + min_y - y[i];

    double width = max_x - min_x;
    double height = max_y - min_y;

    if (rImageWidth.has_value() && width != 0.0) {
        const double scaling_factor = *rImageWidth / width;
        min_x *= scaling_factor;
        min_y *= scaling_factor;
        width *= scaling_factor;
        height *= scaling_factor;
        for (std::size_t i = 0; i < num_points; ++i) {
            x[i] *= scaling_factor;
            y[i] *= scaling_factor;
        }
    }

    std::string stroke_width;
    if (rStrokeWidth.has_value()) {
        stroke_width = *rStrokeWidth;
    } else {
        char buf[64];
        std::snprintf(buf, sizeof(buf), "%g", width / 100.0);
        stroke_width = buf;
    }

    std::ofstream os(rPath, std::ios::binary);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);

    // viewBox: "min_x min_y width height", each float_fmt-formatted.
    os << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" viewBox=\""
       << svg_fmt_num(min_x, rFloatFmt) << ' ' << svg_fmt_num(min_y, rFloatFmt) << ' '
       << svg_fmt_num(width, rFloatFmt) << ' ' << svg_fmt_num(height, rFloatFmt) << "\">";

    // Use path (not polygon): svgo rewrites polygons to paths but drops style.
    os << "<style>path {fill: " << rFill << "; stroke: " << rStroke
       << "; stroke-width: " << stroke_width << "; stroke-linejoin:bevel}</style>";

    for (const auto cb : rMesh.CellRange()) {
        const std::string& type = cb.Type();
        if (type != "line" && type != "triangle" && type != "quad")
            continue;

        const NDArray& conn = cb.Conn();
        const std::size_t ncols = detail::cols(conn);
        const std::size_t n = cb.NumCells();
        for (std::size_t r = 0; r < n; ++r) {
            std::string d;
            for (std::size_t k = 0; k < ncols; ++k) {
                const std::int64_t p = detail::read_int(conn, r * ncols + k);
                // "M x y" for the first vertex, "L x y" for the rest — no
                // separating space before the command letter (matches the
                // Python reference's concatenated format strings).
                d += (k == 0) ? "M " : "L ";
                d += svg_fmt_num(x[static_cast<std::size_t>(p)], rFloatFmt);
                d += ' ';
                d += svg_fmt_num(y[static_cast<std::size_t>(p)], rFloatFmt);
            }
            // triangle/quad are closed; line stays open.
            if (type != "line")
                d += "Z";
            os << "<path d=\"" << d << "\" />";
        }
    }

    os << "</svg>";
}

}  // namespace meshioplusplus
