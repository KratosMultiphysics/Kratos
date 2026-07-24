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
#include "meshioplusplus/formats/tikz.hpp"
#include "meshioplusplus/skin.hpp"

namespace meshioplusplus {

namespace {

// Format a single double with a printf-style spec (spec without leading '%',
// e.g. ".6f").
std::string tikz_fmt_num(double value, const std::string& rSpec) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), ("%" + rSpec).c_str(), value);
    return buf;
}

// 3D rendering path: project the (already skin-extracted) surface mesh with
// the orthographic camera and emit its faces back-to-front as \draw commands.
// Mirrors the flat path's emission; no y-flip (TikZ is y-up).
void tikz_proj_write(const std::string& rPath, const Mesh& rDrawMesh, const std::string& rFloatFmt,
                     bool Standalone, const std::optional<std::string>& rLineWidth,
                     const std::string& rFill, const std::string& rDraw,
                     const std::optional<double>& rScale, double azimuth, double elevation,
                     double roll) {
    const detail::ProjectedSurface ps =
        detail::project_surface(rDrawMesh, azimuth, elevation, roll);

    std::string fill_style = "fill=" + rFill + ", draw=" + rDraw;
    std::string line_style = "draw=" + rDraw;
    if (rLineWidth.has_value()) {
        fill_style += ", line width=" + *rLineWidth;
        line_style += ", line width=" + *rLineWidth;
    }

    std::vector<std::string> lines;
    for (const detail::ProjectedFace& face : ps.mFaces) {
        std::string path;
        for (std::uint8_t k = 0; k < face.mNumNodes; ++k) {
            if (k)
                path += " -- ";
            const std::size_t p = static_cast<std::size_t>(face.mNodes[k]);
            path += "(" + tikz_fmt_num(ps.mX[p], rFloatFmt) + "," +
                    tikz_fmt_num(ps.mY[p], rFloatFmt) + ")";
        }
        if (face.mIsLine)
            lines.push_back("  \\draw[" + line_style + "] " + path + ";");
        else
            lines.push_back("  \\draw[" + fill_style + "] " + path + " -- cycle;");
    }

    std::string pic_opts;
    if (rScale.has_value()) {
        char buf[64];
        std::snprintf(buf, sizeof(buf), "scale=%g", *rScale);
        pic_opts = buf;
    }
    if (rLineWidth.has_value()) {
        if (!pic_opts.empty())
            pic_opts += ", ";
        pic_opts += "line width=" + *rLineWidth;
    }
    const std::string pic_opt_str = pic_opts.empty() ? "" : ("[" + pic_opts + "]");

    std::vector<std::string> out;
    if (Standalone) {
        out.push_back("\\documentclass{standalone}");
        out.push_back("\\usepackage{tikz}");
        out.push_back("\\begin{document}");
    }
    out.push_back("\\begin{tikzpicture}" + pic_opt_str);
    for (const auto& l : lines)
        out.push_back(l);
    out.push_back("\\end{tikzpicture}");
    if (Standalone)
        out.push_back("\\end{document}");

    std::ofstream os(rPath, std::ios::binary);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);

    for (std::size_t i = 0; i < out.size(); ++i)
        os << out[i] << '\n';
}

}  // namespace

void write_tikz(const std::string& rPath, const Mesh& rMesh, const std::string& rFloatFmt,
                bool Standalone, const std::optional<std::string>& rLineWidth,
                const std::string& rFill, const std::string& rDraw,
                const std::optional<double>& rScale, double azimuth, double elevation,
                double roll) {
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
                    tikz_proj_write(rPath, extract_skin(rMesh, /*linearize=*/true), rFloatFmt,
                                    Standalone, rLineWidth, rFill, rDraw, rScale, azimuth,
                                    elevation, roll);
                } else {
                    tikz_proj_write(rPath, rMesh, rFloatFmt, Standalone, rLineWidth, rFill, rDraw,
                                    rScale, azimuth, elevation, roll);
                }
                return;
            }
        }
    }

    // TikZ/PGF uses the math convention (y-up), so — unlike SVG — no y-flip.
    auto coord = [&](std::int64_t p) {
        const std::size_t idx = static_cast<std::size_t>(p);
        const double px = (0 < dim) ? detail::read_double(points, idx * dim + 0) : 0.0;
        const double py = (1 < dim) ? detail::read_double(points, idx * dim + 1) : 0.0;
        return "(" + tikz_fmt_num(px, rFloatFmt) + "," + tikz_fmt_num(py, rFloatFmt) + ")";
    };

    // Per-path style option lists.
    std::string fill_style = "fill=" + rFill + ", draw=" + rDraw;
    std::string line_style = "draw=" + rDraw;
    if (rLineWidth.has_value()) {
        fill_style += ", line width=" + *rLineWidth;
        line_style += ", line width=" + *rLineWidth;
    }

    std::vector<std::string> lines;
    for (const auto cb : rMesh.CellRange()) {
        const std::string& type = cb.Type();
        if (type != "line" && type != "triangle" && type != "quad")
            continue;

        const NDArray& conn = cb.Conn();
        const std::size_t ncols = detail::cols(conn);
        const std::size_t n = cb.NumCells();
        for (std::size_t r = 0; r < n; ++r) {
            std::string path;
            for (std::size_t k = 0; k < ncols; ++k) {
                if (k)
                    path += " -- ";
                path += coord(detail::read_int(conn, r * ncols + k));
            }
            if (type == "line")
                lines.push_back("  \\draw[" + line_style + "] " + path + ";");
            else
                lines.push_back("  \\draw[" + fill_style + "] " + path + " -- cycle;");
        }
    }

    // tikzpicture options (scale / line width) — emitted only when set.
    std::string pic_opts;
    if (rScale.has_value()) {
        char buf[64];
        std::snprintf(buf, sizeof(buf), "scale=%g", *rScale);
        pic_opts = buf;
    }
    if (rLineWidth.has_value()) {
        if (!pic_opts.empty())
            pic_opts += ", ";
        pic_opts += "line width=" + *rLineWidth;
    }
    const std::string pic_opt_str = pic_opts.empty() ? "" : ("[" + pic_opts + "]");

    std::vector<std::string> out;
    if (Standalone) {
        out.push_back("\\documentclass{standalone}");
        out.push_back("\\usepackage{tikz}");
        out.push_back("\\begin{document}");
    }
    out.push_back("\\begin{tikzpicture}" + pic_opt_str);
    for (const auto& l : lines)
        out.push_back(l);
    out.push_back("\\end{tikzpicture}");
    if (Standalone)
        out.push_back("\\end{document}");

    std::ofstream os(rPath, std::ios::binary);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);

    for (std::size_t i = 0; i < out.size(); ++i)
        os << out[i] << '\n';
}

}  // namespace meshioplusplus
