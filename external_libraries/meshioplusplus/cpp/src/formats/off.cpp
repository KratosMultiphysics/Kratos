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
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/obj_off.hpp"

namespace meshioplusplus {

namespace {

std::string off_strip(const std::string& rS) {
    std::size_t b = 0, e = rS.size();
    while (b < e && std::isspace(static_cast<unsigned char>(rS[b])))
        ++b;
    while (e > b && std::isspace(static_cast<unsigned char>(rS[e - 1])))
        --e;
    return rS.substr(b, e - b);
}

}  // namespace

Mesh read_off(const std::string& rPath) {
    std::ifstream in(rPath);
    if (!in)
        throw ReadError("Could not open file: " + rPath);

    std::string line;
    if (!std::getline(in, line) || off_strip(line) != "OFF")
        throw ReadError("Expected the first line to be 'OFF'");

    // Skip comments / blank lines to the counts line.
    std::string counts;
    while (std::getline(in, line)) {
        std::string s = off_strip(line);
        if (!s.empty() && s[0] != '#') {
            counts = s;
            break;
        }
    }
    std::istringstream cs(counts);
    long long num_verts = 0, num_faces = 0, num_edges = 0;
    cs >> num_verts >> num_faces >> num_edges;

    Mesh mesh;
    NDArray pts(DType::Float64, {static_cast<std::size_t>(num_verts), 3});
    double* pp = pts.As<double>();
    for (long long i = 0; i < num_verts * 3; ++i) {
        if (!(in >> pp[i]))
            throw ReadError("OFF: not enough vertex coordinates");
    }
    mesh.AssignPoints(std::move(pts));

    NDArray cells(DType::Int64, {static_cast<std::size_t>(num_faces), 3});
    std::int64_t* cp = cells.As<std::int64_t>();
    for (long long f = 0; f < num_faces; ++f) {
        long long n;
        if (!(in >> n))
            throw ReadError("OFF: not enough faces");
        if (n != 3)
            throw ReadError("OFF: can only read triangular faces");
        in >> cp[f * 3 + 0] >> cp[f * 3 + 1] >> cp[f * 3 + 2];
    }
    mesh.AddCellBlock("triangle", std::move(cells));
    return mesh;
}

void write_off(const std::string& rPath, const Mesh& rMesh) {
    std::ofstream os(rPath, std::ios::binary);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);

    const NDArray& points = rMesh.Points();
    const std::size_t num_points = rMesh.NumPoints();
    const std::size_t dim = rMesh.PointDim();

    // Gather triangles (OFF supports triangles only).
    std::vector<std::int64_t> tri;
    std::size_t ntri = 0;
    for (const auto cb : rMesh.CellRange()) {
        if (cb.Type() != "triangle")
            continue;
        const NDArray& conn = cb.Conn();
        for (std::size_t r = 0; r < cb.NumCells(); ++r) {
            for (int k = 0; k < 3; ++k)
                tri.push_back(detail::read_int(conn, r * 3 + k));
            ++ntri;
        }
    }

    os << "OFF\n# Created by meshio++ (C++ core)\n\n";
    os << num_points << ' ' << ntri << " 0\n\n";

    char buf[96];
    for (std::size_t r = 0; r < num_points; ++r) {
        double x = (0 < dim) ? detail::read_double(points, r * dim + 0) : 0.0;
        double y = (1 < dim) ? detail::read_double(points, r * dim + 1) : 0.0;
        double z = (2 < dim) ? detail::read_double(points, r * dim + 2) : 0.0;
        std::snprintf(buf, sizeof(buf), "%.17g %.17g %.17g\n", x, y, z);
        os << buf;
    }
    for (std::size_t t = 0; t < ntri; ++t)
        os << "3 " << tri[t * 3] << ' ' << tri[t * 3 + 1] << ' ' << tri[t * 3 + 2] << '\n';
}

}  // namespace meshioplusplus
