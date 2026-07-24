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

struct FaceBlock {
    std::size_t mSize = 0;
    std::vector<std::int64_t> mIdx;  // flat, 0-based
    std::vector<std::int64_t> mGids;
    std::size_t mCount = 0;
};

std::string cell_type_for(std::size_t n) {
    if (n == 3)
        return "triangle";
    if (n == 4)
        return "quad";
    return "polygon";
}

NDArray make_point_data(const std::vector<std::vector<double>>& rRows) {
    std::size_t n = rRows.size();
    std::size_t nc = n ? rRows[0].size() : 0;
    NDArray a(DType::Float64, {n, nc});
    double* p = a.As<double>();
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < nc; ++j)
            p[i * nc + j] = rRows[i][j];
    return a;
}

}  // namespace

Mesh read_obj(const std::string& rPath) {
    std::ifstream in(rPath);
    if (!in)
        throw ReadError("Could not open file: " + rPath);

    std::vector<std::array<double, 3>> points;
    std::vector<std::vector<double>> vn, vt;
    std::vector<FaceBlock> blocks;
    std::int64_t group_id = -1;

    std::string line;
    while (std::getline(in, line)) {
        // strip
        std::size_t b = 0, e = line.size();
        while (b < e && std::isspace(static_cast<unsigned char>(line[b])))
            ++b;
        while (e > b && std::isspace(static_cast<unsigned char>(line[e - 1])))
            --e;
        if (b == e || line[b] == '#')
            continue;

        std::istringstream iss(line.substr(b, e - b));
        std::string tag;
        iss >> tag;
        if (tag == "v") {
            std::array<double, 3> p{0, 0, 0};
            iss >> p[0] >> p[1] >> p[2];
            points.push_back(p);
        } else if (tag == "vn") {
            std::vector<double> row;
            double x;
            while (iss >> x)
                row.push_back(x);
            vn.push_back(row);
        } else if (tag == "vt") {
            std::vector<double> row;
            double x;
            while (iss >> x)
                row.push_back(x);
            vt.push_back(row);
        } else if (tag == "f") {
            std::vector<std::int64_t> dat;
            std::string item;
            while (iss >> item) {
                std::size_t slash = item.find('/');
                std::string num = (slash == std::string::npos) ? item : item.substr(0, slash);
                dat.push_back(static_cast<std::int64_t>(std::stoll(num)) - 1);
            }
            std::size_t sz = dat.size();
            if (blocks.empty() || (blocks.back().mCount > 0 && blocks.back().mSize != sz)) {
                FaceBlock fb;
                fb.mSize = sz;
                blocks.push_back(std::move(fb));
            }
            FaceBlock& cur = blocks.back();
            if (cur.mCount == 0)
                cur.mSize = sz;
            cur.mIdx.insert(cur.mIdx.end(), dat.begin(), dat.end());
            cur.mGids.push_back(group_id);
            ++cur.mCount;
        } else if (tag == "g") {
            FaceBlock fb;
            blocks.push_back(std::move(fb));
            ++group_id;
        }
        // 's' and others: ignored.
    }

    // Drop empty blocks (e.g. from trailing 'g').
    std::vector<FaceBlock> nonempty;
    for (auto& fb : blocks)
        if (fb.mCount > 0)
            nonempty.push_back(std::move(fb));

    Mesh mesh;
    std::size_t np = points.size();
    NDArray pts(DType::Float64, {np, 3});
    double* pp = pts.As<double>();
    for (std::size_t i = 0; i < np; ++i)
        for (int c = 0; c < 3; ++c)
            pp[i * 3 + c] = points[i][c];
    mesh.AssignPoints(std::move(pts));

    if (!vt.empty())
        mesh.AddPointData("obj:vt", make_point_data(vt));
    if (!vn.empty())
        mesh.AddPointData("obj:vn", make_point_data(vn));

    if (!nonempty.empty()) {
        std::vector<NDArray> gid_blocks;
        for (auto& fb : nonempty) {
            NDArray data(DType::Int64, {fb.mCount, fb.mSize});
            std::int64_t* dp = data.As<std::int64_t>();
            for (std::size_t i = 0; i < fb.mIdx.size(); ++i)
                dp[i] = fb.mIdx[i];
            mesh.AddCellBlock(cell_type_for(fb.mSize), std::move(data));

            NDArray g(DType::Int64, {fb.mCount});
            std::int64_t* gp = g.As<std::int64_t>();
            for (std::size_t i = 0; i < fb.mCount; ++i)
                gp[i] = fb.mGids[i];
            gid_blocks.push_back(std::move(g));
        }
        mesh.AddCellData("obj:group_ids", std::move(gid_blocks));
    }
    return mesh;
}

void write_obj(const std::string& rPath, const Mesh& rMesh) {
    for (const auto cb : rMesh.CellRange())
        if (cb.Type() != "triangle" && cb.Type() != "quad" && cb.Type() != "polygon")
            throw WriteError(
                "Wavefront .obj files can only contain triangle, quad, "
                "or polygon cells.");

    std::ofstream os(rPath, std::ios::binary);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);

    const NDArray& points = rMesh.Points();
    const std::size_t num_points = rMesh.NumPoints();
    const std::size_t dim = rMesh.PointDim();

    os << "# Created by meshio++ (C++ core)\n";
    char buf[96];
    for (std::size_t r = 0; r < num_points; ++r) {
        double x = (0 < dim) ? detail::read_double(points, r * dim + 0) : 0.0;
        double y = (1 < dim) ? detail::read_double(points, r * dim + 1) : 0.0;
        double z = (2 < dim) ? detail::read_double(points, r * dim + 2) : 0.0;
        std::snprintf(buf, sizeof(buf), "v %.17g %.17g %.17g\n", x, y, z);
        os << buf;
    }

    auto write_pd = [&](const char* key, const char* tag) {
        if (!rMesh.HasPointData(key))
            return;
        const NDArray& d = rMesh.PointData(key);
        std::size_t nc = d.Shape().size() >= 2 ? d.Shape()[1] : 1;
        for (std::size_t r = 0; r < (d.Shape().empty() ? 0 : d.Shape()[0]); ++r) {
            os << tag;
            for (std::size_t c = 0; c < nc; ++c) {
                std::snprintf(buf, sizeof(buf), " %.17g", detail::read_double(d, r * nc + c));
                os << buf;
            }
            os << '\n';
        }
    };
    write_pd("obj:vn", "vn");
    write_pd("obj:vt", "vt");

    for (const auto cb : rMesh.CellRange()) {
        const NDArray& conn = cb.Conn();
        std::size_t k = conn.Shape().size() >= 2 ? conn.Shape()[1] : 1;
        for (std::size_t r = 0; r < cb.NumCells(); ++r) {
            os << 'f';
            for (std::size_t j = 0; j < k; ++j)
                os << ' ' << (detail::read_int(conn, r * k + j) + 1);
            os << '\n';
        }
    }
}

}  // namespace meshioplusplus
