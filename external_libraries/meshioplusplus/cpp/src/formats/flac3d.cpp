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
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/flac3d.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {

namespace {

// meshio type -> simplified FLAC3D base type (zone = 3D, face = 2D), or "".
std::string zone_key(const std::string& rT) {
    static const std::unordered_map<std::string, std::string> m = {{"tetra", "tetra"},
                                                                   {"tetra10", "tetra"},
                                                                   {"pyramid", "pyramid"},
                                                                   {"pyramid13", "pyramid"},
                                                                   {"wedge", "wedge"},
                                                                   {"wedge12", "wedge"},
                                                                   {"wedge15", "wedge"},
                                                                   {"wedge18", "wedge"},
                                                                   {"hexahedron", "hexahedron"},
                                                                   {"hexahedron20", "hexahedron"},
                                                                   {"hexahedron24", "hexahedron"},
                                                                   {"hexahedron27", "hexahedron"}};
    auto it = m.find(rT);
    return it == m.end() ? std::string() : it->second;
}
std::string face_key(const std::string& rT) {
    static const std::unordered_map<std::string, std::string> m = {
        {"triangle", "triangle"}, {"triangle6", "triangle"}, {"triangle7", "triangle"},
        {"quad", "quad"},         {"quad8", "quad"},         {"quad9", "quad"}};
    auto it = m.find(rT);
    return it == m.end() ? std::string() : it->second;
}

const std::unordered_map<int, std::string>& numnodes_type(int dim) {
    static const std::unordered_map<int, std::string> z = {
        {4, "tetra"}, {5, "pyramid"}, {6, "wedge"}, {8, "hexahedron"}};
    static const std::unordered_map<int, std::string> fc = {{3, "triangle"}, {4, "quad"}};
    return dim == 3 ? z : fc;
}

const char* flac3d_type(const std::string& rKey) {
    if (rKey == "triangle")
        return "T3";
    if (rKey == "quad")
        return "Q4";
    if (rKey == "tetra")
        return "T4";
    if (rKey == "pyramid")
        return "P5";
    if (rKey == "wedge")
        return "W6";
    return "B8";  // hexahedron
}

const std::vector<int>& f2m_order(const std::string& rKey) {
    static const std::unordered_map<std::string, std::vector<int>> m = {
        {"triangle", {0, 1, 2}},       {"quad", {0, 1, 2, 3}},
        {"tetra", {0, 1, 2, 3}},       {"pyramid", {0, 1, 4, 2, 3}},
        {"wedge", {0, 1, 3, 2, 4, 5}}, {"hexahedron", {0, 1, 4, 2, 3, 6, 7, 5}}};
    return m.at(rKey);
}
const std::vector<int>& m2f_order(const std::string& rKey) {
    static const std::unordered_map<std::string, std::vector<int>> m = {
        {"triangle", {0, 1, 2}},       {"quad", {0, 1, 2, 3}},
        {"tetra", {0, 1, 2, 3}},       {"pyramid", {0, 1, 3, 4, 2}},
        {"wedge", {0, 1, 3, 2, 4, 5}}, {"hexahedron", {0, 1, 3, 4, 2, 7, 5, 6}}};
    return m.at(rKey);
}
const std::vector<int>& m2f_order2(const std::string& rKey) {
    static const std::unordered_map<std::string, std::vector<int>> m = {
        {"tetra", {0, 2, 1, 3}},
        {"pyramid", {0, 3, 1, 4, 2}},
        {"wedge", {0, 2, 3, 1, 5, 4}},
        {"hexahedron", {0, 3, 1, 4, 2, 5, 7, 6}}};
    return m.at(rKey);
}

// little-endian binary scalar I/O (host assumed little-endian)
std::uint32_t ru32(std::istream& rIn) {
    std::uint32_t v;
    rIn.read(reinterpret_cast<char*>(&v), 4);
    if (rIn.gcount() != 4)
        throw ReadError("FLAC3D: unexpected end of file");
    return v;
}
double rf64(std::istream& rIn) {
    double v;
    rIn.read(reinterpret_cast<char*>(&v), 8);
    if (rIn.gcount() != 8)
        throw ReadError("FLAC3D: unexpected end of file");
    return v;
}
void wu32(std::ostream& rOs, std::uint32_t v) {
    rOs.write(reinterpret_cast<const char*>(&v), 4);
}
void wf64(std::ostream& rOs, double v) {
    rOs.write(reinterpret_cast<const char*>(&v), 8);
}

// Accumulating raw cell block: meshio node order will be applied later.
struct Flac3dRawBlock {
    std::string mType;                             // meshio type
    std::vector<std::vector<std::int64_t>> mRows;  // 0-based point indices
};

void add_cell(std::vector<Flac3dRawBlock>& rBlocks, const std::string& rType,
              std::vector<std::int64_t>&& cell) {
    if (rBlocks.empty() || rBlocks.back().mType != rType)
        rBlocks.push_back(Flac3dRawBlock{rType, {}});
    rBlocks.back().mRows.push_back(std::move(cell));
}

std::vector<std::string> flac3d_split_ws(const std::string& rS) {
    std::vector<std::string> out;
    std::istringstream iss(rS);
    std::string t;
    while (iss >> t)
        out.push_back(t);
    return out;
}

}  // namespace

Mesh read_flac3d(const std::string& rPath) {
    // Sniff binary (a null byte in the first 8 bytes).
    bool binary = false;
    {
        std::ifstream sniff(rPath, std::ios::binary);
        if (!sniff)
            throw ReadError("Could not open file: " + rPath);
        char block[8] = {0};
        sniff.read(block, 8);
        std::streamsize got = sniff.gcount();
        for (std::streamsize i = 0; i < got; ++i)
            if (block[i] == '\0') {
                binary = true;
                break;
            }
    }

    std::vector<double> points;                                // flat xyz
    std::unordered_map<std::int64_t, std::int64_t> point_ids;  // file id -> index
    std::vector<Flac3dRawBlock> z_blocks, f_blocks;
    std::vector<std::int64_t> z_ids, f_ids;

    if (binary) {
        std::ifstream in(rPath, std::ios::binary);
        char hdr[8];
        in.read(hdr, 8);  // unknown header
        std::uint32_t num_nodes = ru32(in);
        points.reserve(num_nodes * 3);
        for (std::uint32_t i = 0; i < num_nodes; ++i) {
            std::uint32_t pid = ru32(in);
            double x = rf64(in), y = rf64(in), z = rf64(in);
            point_ids[pid] = static_cast<std::int64_t>(i);
            points.push_back(x);
            points.push_back(y);
            points.push_back(z);
        }
        for (int fi = 0; fi < 2; ++fi) {
            int dim = (fi == 0) ? 3 : 2;
            std::vector<Flac3dRawBlock>& blocks = (fi == 0) ? z_blocks : f_blocks;
            std::vector<std::int64_t>& ids = (fi == 0) ? z_ids : f_ids;
            std::uint32_t num_cells = ru32(in);
            const auto& tmap = numnodes_type(dim);
            for (std::uint32_t k = 0; k < num_cells; ++k) {
                std::uint32_t cid = ru32(in);
                std::uint32_t nv = ru32(in);
                std::vector<std::int64_t> cell(nv);
                for (std::uint32_t j = 0; j < nv; ++j)
                    cell[j] = point_ids.at(ru32(in));
                if (nv == 7)
                    cell.push_back(cell.back());
                auto it = tmap.find(static_cast<int>(cell.size()));
                if (it == tmap.end())
                    throw ReadError("FLAC3D: bad cell node count");
                ids.push_back(cid);
                add_cell(blocks, it->second, std::move(cell));
            }
            std::uint32_t num_groups = ru32(in);
            if (num_groups > 0)
                throw ReadError("FLAC3D: cell groups handled by Python fallback");
        }
    } else {
        std::ifstream in(rPath, std::ios::binary);
        std::string line;
        while (std::getline(in, line)) {
            std::vector<std::string> s = flac3d_split_ws(line);
            if (s.empty())
                continue;
            if (s[0] == "G") {
                std::int64_t pid = std::strtoll(s[1].c_str(), nullptr, 10);
                point_ids[pid] = static_cast<std::int64_t>(points.size() / 3);
                for (std::size_t j = 2; j < s.size(); ++j)
                    points.push_back(std::strtod(s[j].c_str(), nullptr));
            } else if (s[0] == "Z" || s[0] == "F") {
                int dim = (s[0] == "Z") ? 3 : 2;
                std::int64_t cid = std::strtoll(s[2].c_str(), nullptr, 10);
                bool is_b7 = (s[1] == "B7");
                std::vector<std::int64_t> cell;
                for (std::size_t j = 3; j < s.size(); ++j)
                    cell.push_back(point_ids.at(std::strtoll(s[j].c_str(), nullptr, 10)));
                if (is_b7)
                    cell.push_back(cell.back());
                const auto& tmap = numnodes_type(dim);
                auto it = tmap.find(static_cast<int>(cell.size()));
                if (it == tmap.end())
                    throw ReadError("FLAC3D: bad cell node count");
                if (dim == 3) {
                    z_ids.push_back(cid);
                    add_cell(z_blocks, it->second, std::move(cell));
                } else {
                    f_ids.push_back(cid);
                    add_cell(f_blocks, it->second, std::move(cell));
                }
            } else if (s[0] == "ZGROUP" || s[0] == "FGROUP") {
                throw ReadError("FLAC3D: cell groups handled by Python fallback");
            }
            // other lines (comments starting with '*') are ignored
        }
    }

    // Assemble: faces first, then zones (matching the Python reader).
    Mesh mesh;
    const std::int64_t npoints = static_cast<std::int64_t>(points.size() / 3);
    NDArray pts(DType::Float64, {static_cast<std::size_t>(npoints), 3});
    std::memcpy(pts.Data(), points.data(), points.size() * sizeof(double));
    mesh.AssignPoints(std::move(pts));

    std::vector<std::size_t> block_sizes;
    auto emit = [&](std::vector<Flac3dRawBlock>& blocks) {
        for (auto& b : blocks) {
            const std::vector<int>& ord =
                f2m_order(zone_key(b.mType).empty() ? face_key(b.mType) : zone_key(b.mType));
            std::size_t n = b.mRows.size();
            std::size_t k = ord.size();
            NDArray data(DType::Int64, {n, k});
            std::int64_t* dp = data.As<std::int64_t>();
            for (std::size_t r = 0; r < n; ++r)
                for (std::size_t j = 0; j < k; ++j)
                    dp[r * k + j] = b.mRows[r][ord[j]];
            mesh.AddCellBlock(b.mType, std::move(data));
            block_sizes.push_back(n);
        }
    };
    emit(f_blocks);
    emit(z_blocks);

    // Global cell ids -> cell_data["cell_ids"], split per block.
    if (mesh.NumCellBlocks() != 0) {
        std::int64_t z_offset = static_cast<std::int64_t>(f_ids.size());
        std::vector<std::int64_t> all_ids;
        all_ids.reserve(f_ids.size() + z_ids.size());
        for (auto v : f_ids)
            all_ids.push_back(v);
        for (auto v : z_ids)
            all_ids.push_back(v + z_offset);

        std::vector<NDArray> id_blocks;
        std::size_t off = 0;
        for (std::size_t sz : block_sizes) {
            NDArray a(DType::Int64, {sz});
            for (std::size_t r = 0; r < sz; ++r)
                a.As<std::int64_t>()[r] = all_ids[off + r];
            off += sz;
            id_blocks.push_back(std::move(a));
        }
        mesh.AddCellData("cell_ids", std::move(id_blocks));
    }

    return mesh;
}

namespace {

// Reorder one zone cell to FLAC3D order, choosing the right-handed permutation
// via the scalar triple product of the first four ordered corners.
std::vector<std::int64_t> zone_cell_flac3d(const NDArray& rPoints, const NDArray& rData,
                                           std::size_t row, const std::string& rKey) {
    const std::vector<int>& o1 = m2f_order(rKey);
    const std::vector<int>& o2 = m2f_order2(rKey);
    const std::size_t ncols = detail::cols(rData);

    auto node = [&](int local) -> std::int64_t {
        return detail::read_int(rData, row * ncols + local);
    };
    auto coord = [&](std::int64_t p, int c) -> double {
        return detail::read_double(rPoints, static_cast<std::size_t>(p) * 3 + c);
    };

    // first four corners in FLAC3D order
    std::int64_t c0 = node(o1[0]), c1 = node(o1[1]), c2 = node(o1[2]), c3 = node(o1[3]);
    double a[3], b[3], c[3];
    for (int i = 0; i < 3; ++i) {
        a[i] = coord(c1, i) - coord(c0, i);
        b[i] = coord(c2, i) - coord(c0, i);
        c[i] = coord(c3, i) - coord(c0, i);
    }
    double cross0 = b[1] * c[2] - b[2] * c[1];
    double cross1 = b[2] * c[0] - b[0] * c[2];
    double cross2 = b[0] * c[1] - b[1] * c[0];
    double det = a[0] * cross0 + a[1] * cross1 + a[2] * cross2;

    const std::vector<int>& ord = (det > 0) ? o1 : o2;
    std::vector<std::int64_t> out(ord.size());
    for (std::size_t j = 0; j < ord.size(); ++j)
        out[j] = node(ord[j]);
    return out;
}

}  // namespace

void write_flac3d(const std::string& rPath, const Mesh& rMesh, const std::string& rFloatFmt,
                  bool binary) {
    // Split blocks by FLAC3D category.
    std::vector<std::size_t> zone_idx, face_idx;
    for (std::size_t i = 0; i < rMesh.NumCellBlocks(); ++i) {
        if (!zone_key(rMesh.Cells(i).Type()).empty())
            zone_idx.push_back(i);
        else if (!face_key(rMesh.Cells(i).Type()).empty())
            face_idx.push_back(i);
    }

    std::ofstream f(rPath, std::ios::binary);
    if (!f)
        throw WriteError("Could not open file for writing: " + rPath);

    const std::size_t npts = rMesh.NumPoints();
    const std::size_t pdim = rMesh.PointDim();
    const NDArray& points = rMesh.Points();

    if (binary) {
        wu32(f, 1375135718u);
        wu32(f, 3u);
        // points
        wu32(f, static_cast<std::uint32_t>(npts));
        for (std::size_t i = 0; i < npts; ++i) {
            wu32(f, static_cast<std::uint32_t>(i + 1));
            for (int c = 0; c < 3; ++c)
                wf64(f, c < static_cast<int>(pdim) ? detail::read_double(points, i * pdim + c)
                                                   : 0.0);
        }
        std::uint32_t gid = 0;
        // zones
        std::uint32_t nz = 0;
        for (auto i : zone_idx)
            nz += static_cast<std::uint32_t>(rMesh.Cells(i).NumCells());
        wu32(f, nz);
        for (auto i : zone_idx) {
            const auto cb = rMesh.Cells(i);
            const NDArray& conn = cb.Conn();
            std::string key = zone_key(cb.Type());
            std::size_t n = cb.NumCells();
            // Right-handed reorder per row is independent -> compute in
            // parallel, then stream sequentially.
            std::vector<std::vector<std::int64_t>> zcells(n);
            parallel_for(n, [&](std::size_t r) {
                zcells[r] = zone_cell_flac3d(points, conn, r, key);
            });
            for (std::size_t r = 0; r < n; ++r) {
                const auto& cell = zcells[r];
                wu32(f, ++gid);
                wu32(f, static_cast<std::uint32_t>(cell.size()));
                for (auto v : cell)
                    wu32(f, static_cast<std::uint32_t>(v + 1));
            }
        }
        wu32(f, 0u);  // zone groups
        // faces
        std::uint32_t nf = 0;
        for (auto i : face_idx)
            nf += static_cast<std::uint32_t>(rMesh.Cells(i).NumCells());
        wu32(f, nf);
        for (auto i : face_idx) {
            const auto cb = rMesh.Cells(i);
            const NDArray& conn = cb.Conn();
            std::string key = face_key(cb.Type());
            const std::vector<int>& ord = m2f_order(key);
            std::size_t n = cb.NumCells();
            std::size_t ncols = detail::cols(conn);
            for (std::size_t r = 0; r < n; ++r) {
                wu32(f, ++gid);
                wu32(f, static_cast<std::uint32_t>(ord.size()));
                for (int local : ord)
                    wu32(f,
                         static_cast<std::uint32_t>(detail::read_int(conn, r * ncols + local) + 1));
            }
        }
        wu32(f, 0u);  // face groups
        return;
    }

    // ASCII
    f << "* FLAC3D grid produced by meshio++ (C++ core)\n";
    f << "* GRIDPOINTS\n";
    char buf[64];
    for (std::size_t i = 0; i < npts; ++i) {
        f << "G\t" << (i + 1) << "\t";
        for (int c = 0; c < 3; ++c) {
            double v =
                c < static_cast<int>(pdim) ? detail::read_double(points, i * pdim + c) : 0.0;
            std::snprintf(buf, sizeof(buf), ("%" + rFloatFmt).c_str(), v);
            f << buf << (c == 2 ? '\n' : '\t');
        }
    }

    std::int64_t gid = 0;
    f << "* ZONES\n";
    for (auto i : zone_idx) {
        const auto cb = rMesh.Cells(i);
        const NDArray& conn = cb.Conn();
        std::string key = zone_key(cb.Type());
        const char* abbr = flac3d_type(key);
        std::size_t n = cb.NumCells();
        // Right-handed reorder per row is independent -> compute in parallel,
        // then stream sequentially.
        std::vector<std::vector<std::int64_t>> zcells(n);
        parallel_for(n,
                     [&](std::size_t r) { zcells[r] = zone_cell_flac3d(points, conn, r, key); });
        for (std::size_t r = 0; r < n; ++r) {
            f << "Z " << abbr << " " << (++gid);
            for (auto v : zcells[r])
                f << " " << (v + 1);
            f << "\n";
        }
    }
    f << "* ZONE GROUPS\n";

    f << "* FACES\n";
    for (auto i : face_idx) {
        const auto cb = rMesh.Cells(i);
        const NDArray& conn = cb.Conn();
        std::string key = face_key(cb.Type());
        const char* abbr = flac3d_type(key);
        const std::vector<int>& ord = m2f_order(key);
        std::size_t n = cb.NumCells();
        std::size_t ncols = detail::cols(conn);
        for (std::size_t r = 0; r < n; ++r) {
            f << "F " << abbr << " " << (++gid);
            for (int local : ord)
                f << " " << (detail::read_int(conn, r * ncols + local) + 1);
            f << "\n";
        }
    }
    f << "* FACE GROUPS\n";
}

}  // namespace meshioplusplus
