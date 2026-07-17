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
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

// Project includes
#include "meshioplusplus/formats/stl.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus {

namespace {

// First-occurrence de-duplication of 3-component rows. Returns per-row unique
// index; appends unique rows (raw bytes) to `rOutPoints`.
std::vector<std::int64_t> dedup(const unsigned char* pRows, std::size_t nrows, std::size_t isz,
                                std::vector<unsigned char>& rOutPoints) {
    std::unordered_map<std::string, std::int64_t> seen;
    seen.reserve(nrows);
    std::vector<std::int64_t> idx(nrows);
    const std::size_t rowbytes = 3 * isz;
    for (std::size_t i = 0; i < nrows; ++i) {
        std::string key(reinterpret_cast<const char*>(pRows) + i * rowbytes, rowbytes);
        auto it = seen.find(key);
        if (it == seen.end()) {
            std::int64_t id = static_cast<std::int64_t>(seen.size());
            seen.emplace(std::move(key), id);
            rOutPoints.insert(rOutPoints.end(), pRows + i * rowbytes, pRows + (i + 1) * rowbytes);
            idx[i] = id;
        } else {
            idx[i] = it->second;
        }
    }
    return idx;
}

Mesh build_mesh(std::vector<unsigned char>& rVertBytes, DType dt,
                std::vector<unsigned char>* pNormalBytes) {
    std::size_t isz = dtype_size(dt);
    std::size_t nverts = rVertBytes.size() / (3 * isz);

    std::vector<unsigned char> point_bytes;
    std::vector<std::int64_t> idx = dedup(rVertBytes.data(), nverts, isz, point_bytes);

    Mesh mesh;
    std::size_t num_unique = point_bytes.size() / (3 * isz);
    NDArray pts(dt, {num_unique, 3});
    if (!point_bytes.empty())
        std::memcpy(pts.Data(), point_bytes.data(), point_bytes.size());
    mesh.AssignPoints(std::move(pts));

    std::size_t ntri = nverts / 3;
    // An empty STL has no cells (match the Python reader, which returns no
    // cell blocks rather than an empty triangle block).
    if (ntri == 0)
        return mesh;

    NDArray cells(DType::Int64, {ntri, 3});
    std::int64_t* cp = cells.As<std::int64_t>();
    for (std::size_t i = 0; i < ntri * 3; ++i)
        cp[i] = idx[i];
    mesh.AddCellBlock("triangle", std::move(cells));

    if (pNormalBytes && !pNormalBytes->empty()) {
        NDArray nrm(DType::Float64, {ntri, 3});
        std::memcpy(nrm.Data(), pNormalBytes->data(), pNormalBytes->size());
        mesh.AppendCellData("facet_normals", std::move(nrm));
    }
    return mesh;
}

bool starts_with(const std::string& rS, const char* pP) {
    return rS.rfind(pP, 0) == 0;
}

bool is_comment_line(const std::string& rS) {
    return starts_with(rS, "solid") || starts_with(rS, "outer loop") ||
           starts_with(rS, "endloop") || starts_with(rS, "endfacet") || starts_with(rS, "endsolid");
}

std::string lstrip(const std::string& rS) {
    std::size_t b = 0;
    while (b < rS.size() && std::isspace(static_cast<unsigned char>(rS[b])))
        ++b;
    return rS.substr(b);
}

Mesh read_ascii(std::ifstream& rIn) {
    // Collect the last 3 numbers of every non-comment line; rows 0,4,8,... are
    // facet normals, the rest are vertices.
    std::vector<double> data;
    std::string line;
    while (std::getline(rIn, line)) {
        std::string s = lstrip(line);
        if (s.empty() || is_comment_line(s))
            continue;
        std::istringstream iss(s);
        std::vector<std::string> tok;
        std::string t;
        while (iss >> t)
            tok.push_back(t);
        if (tok.size() < 3)
            continue;
        for (std::size_t j = tok.size() - 3; j < tok.size(); ++j)
            data.push_back(std::strtod(tok[j].c_str(), nullptr));
    }
    std::size_t nrows = data.size() / 3;
    if (nrows % 4 != 0)
        throw ReadError("Malformed ascii STL");

    std::vector<unsigned char> verts, normals;
    for (std::size_t r = 0; r < nrows; ++r) {
        const double* row = data.data() + r * 3;
        std::vector<unsigned char>& dst = (r % 4 == 0) ? normals : verts;
        dst.insert(dst.end(), reinterpret_cast<const unsigned char*>(row),
                   reinterpret_cast<const unsigned char*>(row) + 3 * sizeof(double));
    }
    return build_mesh(verts, DType::Float64, &normals);
}

Mesh read_binary(std::ifstream& rIn, std::uint32_t num_tri) {
    std::vector<unsigned char> verts;
    verts.reserve(num_tri * 9 * sizeof(float));
    unsigned char tri[50];
    for (std::uint32_t i = 0; i < num_tri; ++i) {
        rIn.read(reinterpret_cast<char*>(tri), 50);
        if (rIn.gcount() != 50)
            throw ReadError("Truncated binary STL");
        // bytes [12, 48) are the 9 float32 vertex coords (host is little-endian).
        verts.insert(verts.end(), tri + 12, tri + 48);
    }
    return build_mesh(verts, DType::Float32, nullptr);
}

}  // namespace

Mesh read_stl(const std::string& rPath) {
    std::ifstream in(rPath, std::ios::binary);
    if (!in)
        throw ReadError("Could not open file: " + rPath);
    in.seekg(0, std::ios::end);
    std::streamoff filesize = in.tellg();
    in.seekg(0, std::ios::beg);

    if (filesize < 80)
        return read_ascii(in);

    char header[80];
    in.read(header, 80);
    std::uint32_t num_tri = 0;
    in.read(reinterpret_cast<char*>(&num_tri), 4);  // little-endian host
    if (static_cast<std::streamoff>(84 + std::uint64_t(num_tri) * 50) == filesize)
        return read_binary(in, num_tri);

    // Fall back to ascii: rewind, skip the first line.
    in.clear();
    in.seekg(0, std::ios::beg);
    std::string first;
    std::getline(in, first);
    return read_ascii(in);
}

namespace {

void gather_triangles(const Mesh& rMesh, std::vector<std::array<double, 9>>& rTris,
                      std::vector<std::array<double, 3>>& rNormals) {
    const bool have_normals = rMesh.HasCellData("facet_normals");
    const std::size_t normal_blocks = have_normals ? rMesh.CellDataNumBlocks("facet_normals") : 0;
    const NDArray& points = rMesh.Points();
    std::size_t dim = points.Shape().size() >= 2 ? points.Shape()[1] : 3;

    std::size_t block = 0;
    for (const auto cb : rMesh.CellRange()) {
        if (cb.Type() != "triangle") {
            ++block;
            continue;
        }
        std::size_t nc = cb.NumCells();
        const NDArray& conn = cb.Conn();
        const NDArray* nrm = nullptr;
        if (have_normals && block < normal_blocks)
            nrm = &rMesh.CellData("facet_normals", block);
        for (std::size_t r = 0; r < nc; ++r) {
            std::array<double, 9> tri{};
            double v[3][3];
            for (int k = 0; k < 3; ++k) {
                std::int64_t pi = detail::read_int(conn, r * 3 + k);
                for (int c = 0; c < 3; ++c)
                    v[k][c] =
                        (std::size_t(c) < dim) ? detail::read_double(points, pi * dim + c) : 0.0;
                tri[k * 3 + 0] = v[k][0];
                tri[k * 3 + 1] = v[k][1];
                tri[k * 3 + 2] = v[k][2];
            }
            rTris.push_back(tri);

            std::array<double, 3> n{};
            if (nrm) {
                for (int c = 0; c < 3; ++c)
                    n[c] = detail::read_double(*nrm, r * 3 + c);
            } else {
                double a[3] = {v[1][0] - v[0][0], v[1][1] - v[0][1], v[1][2] - v[0][2]};
                double b[3] = {v[2][0] - v[0][0], v[2][1] - v[0][1], v[2][2] - v[0][2]};
                n[0] = a[1] * b[2] - a[2] * b[1];
                n[1] = a[2] * b[0] - a[0] * b[2];
                n[2] = a[0] * b[1] - a[1] * b[0];
                double len = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
                if (len > 0) {
                    n[0] /= len;
                    n[1] /= len;
                    n[2] /= len;
                }
            }
            rNormals.push_back(n);
        }
        ++block;
    }
}

}  // namespace

void write_stl(const std::string& rPath, const Mesh& rMesh, bool binary) {
    std::vector<std::array<double, 9>> tris;
    std::vector<std::array<double, 3>> normals;
    gather_triangles(rMesh, tris, normals);

    std::ofstream os(rPath, std::ios::binary);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);

    if (binary) {
        char header[80];
        std::memset(header, 'X', 80);
        const char* msg = "meshio++ (C++ core) binary STL";
        std::memcpy(header, msg, std::strlen(msg));
        os.write(header, 80);
        std::uint32_t n = static_cast<std::uint32_t>(tris.size());
        os.write(reinterpret_cast<const char*>(&n), 4);
        for (std::size_t i = 0; i < tris.size(); ++i) {
            float buf[12];
            for (int c = 0; c < 3; ++c)
                buf[c] = static_cast<float>(normals[i][c]);
            for (int c = 0; c < 9; ++c)
                buf[3 + c] = static_cast<float>(tris[i][c]);
            os.write(reinterpret_cast<const char*>(buf), 48);
            std::uint16_t attr = 0;
            os.write(reinterpret_cast<const char*>(&attr), 2);
        }
    } else {
        auto wr3 = [&](const char* prefix, const double* p) {
            char line[160];
            std::snprintf(line, sizeof(line), "%s %.17g %.17g %.17g\n", prefix, p[0], p[1], p[2]);
            os << line;
        };
        os << "solid\n";
        for (std::size_t i = 0; i < tris.size(); ++i) {
            wr3("facet normal", normals[i].data());
            os << " outer loop\n";
            wr3("  vertex", &tris[i][0]);
            wr3("  vertex", &tris[i][3]);
            wr3("  vertex", &tris[i][6]);
            os << " endloop\nendfacet\n";
        }
        os << "endsolid\n";
    }
}

}  // namespace meshioplusplus
