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
#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/su2.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {

namespace {

int su2_numnodes(int t) {
    switch (t) {
        case 3:
            return 2;  // line
        case 5:
            return 3;  // triangle
        case 9:
            return 4;  // quad
        case 10:
            return 4;  // tetra
        case 12:
            return 8;  // hexahedron
        case 13:
            return 6;  // wedge
        case 14:
            return 5;  // pyramid
        default:
            return 0;
    }
}
std::string su2_to_meshio(int t) {
    switch (t) {
        case 3:
            return "line";
        case 5:
            return "triangle";
        case 9:
            return "quad";
        case 10:
            return "tetra";
        case 12:
            return "hexahedron";
        case 13:
            return "wedge";
        case 14:
            return "pyramid";
        default:
            return "";
    }
}
int meshio_to_su2(const std::string& rT) {
    if (rT == "line")
        return 3;
    if (rT == "triangle")
        return 5;
    if (rT == "quad")
        return 9;
    if (rT == "tetra")
        return 10;
    if (rT == "hexahedron")
        return 12;
    if (rT == "wedge")
        return 13;
    if (rT == "pyramid")
        return 14;
    return -1;
}

std::string su2_strip(const std::string& rS) {
    std::size_t b = rS.find_first_not_of(" \t\r\n");
    if (b == std::string::npos)
        return "";
    std::size_t e = rS.find_last_not_of(" \t\r\n");
    return rS.substr(b, e - b + 1);
}
std::vector<std::string> su2_tokens(const std::string& rS) {
    std::vector<std::string> out;
    std::istringstream iss(rS);
    std::string t;
    while (iss >> t)
        out.push_back(t);
    return out;
}

struct Blk {
    std::string mType;
    int mN = 0;
    std::vector<std::int64_t> mConn;
    std::vector<std::int32_t> mTag;
    std::size_t mCount = 0;
};

// Parse `count` element lines (each "vtk_type n0 n1 ... [extra]") into type-
// grouped blocks (sorted by vtk type code, matching numpy.unique), all with
// the given tag.
void read_elem_block(const std::vector<std::string>& rLines, std::size_t& rLi, std::size_t count,
                     std::int32_t tag, std::vector<Blk>& rOut) {
    std::vector<std::pair<int, std::vector<std::int64_t>>> elems;
    std::set<int> types;
    for (std::size_t e = 0; e < count; ++e) {
        auto t = su2_tokens(rLines.at(rLi++));
        int vt = std::stoi(t[0]);
        int nn = su2_numnodes(vt);
        if (nn == 0)
            throw ReadError("SU2: unsupported element type " + t[0]);
        std::vector<std::int64_t> nodes(nn);
        for (int j = 0; j < nn; ++j)
            nodes[j] = std::strtoll(t[1 + j].c_str(), nullptr, 10);
        elems.emplace_back(vt, std::move(nodes));
        types.insert(vt);
    }
    for (int vt : types) {  // std::set is sorted
        Blk b;
        b.mType = su2_to_meshio(vt);
        b.mN = su2_numnodes(vt);
        for (auto& e : elems) {
            if (e.first != vt)
                continue;
            b.mConn.insert(b.mConn.end(), e.second.begin(), e.second.end());
            b.mTag.push_back(tag);
            ++b.mCount;
        }
        rOut.push_back(std::move(b));
    }
}

}  // namespace

Mesh read_su2(const std::string& rPath) {
    std::ifstream in(rPath);
    if (!in)
        throw ReadError("Could not open file: " + rPath);
    std::vector<std::string> lines;
    std::string l;
    while (std::getline(in, l))
        lines.push_back(l);

    int dim = 0;
    Mesh mesh;
    std::vector<Blk> blocks;
    std::int32_t next_tag_id = 0;

    std::size_t li = 0;
    while (li < lines.size()) {
        std::string line = su2_strip(lines[li]);
        if (line.empty() || line[0] == '%') {
            ++li;
            continue;
        }
        std::size_t eq = line.find('=');
        if (eq == std::string::npos) {
            ++li;
            continue;
        }
        std::string name = su2_strip(line.substr(0, eq));
        std::string rest = su2_strip(line.substr(eq + 1));
        ++li;

        if (name == "NDIME") {
            dim = std::stoi(rest);
            if (dim != 2 && dim != 3)
                throw ReadError("SU2: invalid NDIME");
        } else if (name == "NPOIN") {
            std::size_t npoin = static_cast<std::size_t>(std::stoll(su2_tokens(rest)[0]));
            NDArray pts(DType::Float64, {npoin, static_cast<std::size_t>(dim)});
            double* pp = pts.As<double>();
            for (std::size_t i = 0; i < npoin; ++i) {
                auto t = su2_tokens(lines.at(li++));
                for (int c = 0; c < dim; ++c)
                    pp[i * dim + c] = std::strtod(t[c].c_str(), nullptr);
            }
            mesh.AssignPoints(std::move(pts));
        } else if (name == "NELEM") {
            std::size_t ne = static_cast<std::size_t>(std::stoll(rest));
            read_elem_block(lines, li, ne, 0, blocks);
        } else if (name == "NMARK") {
            // handled implicitly via MARKER_TAG/MARKER_ELEMS
        } else if (name == "MARKER_TAG") {
            try {
                std::size_t pos;
                int v = std::stoi(rest, &pos);
                if (pos == rest.size())
                    next_tag_id = v;
                else {
                    ++next_tag_id;
                }
            } catch (...) {
                ++next_tag_id;
            }
        } else if (name == "MARKER_ELEMS") {
            std::size_t ne = static_cast<std::size_t>(std::stoll(rest));
            read_elem_block(lines, li, ne, next_tag_id, blocks);
        }
    }

    // Merge boundary blocks of the same type (lines in 2D; tris/quads in 3D).
    std::vector<std::string> btypes = (dim == 2) ? std::vector<std::string>{"line"}
                                                 : std::vector<std::string>{"triangle", "quad"};
    for (const auto& bt : btypes) {
        int first = -1;
        for (std::size_t i = 0; i < blocks.size(); ++i) {
            if (blocks[i].mType != bt)
                continue;
            if (first < 0) {
                first = static_cast<int>(i);
                continue;
            }
            Blk& dst = blocks[first];
            Blk& src = blocks[i];
            dst.mConn.insert(dst.mConn.end(), src.mConn.begin(), src.mConn.end());
            dst.mTag.insert(dst.mTag.end(), src.mTag.begin(), src.mTag.end());
            dst.mCount += src.mCount;
            src.mCount = 0;  // mark for removal
            src.mConn.clear();
        }
    }

    std::vector<NDArray> tags;
    for (auto& b : blocks) {
        if (b.mCount == 0)
            continue;  // merged-away or empty
        NDArray data(DType::Int64, {b.mCount, static_cast<std::size_t>(b.mN)});
        std::memcpy(data.Data(), b.mConn.data(), b.mConn.size() * sizeof(std::int64_t));
        mesh.AddCellBlock(b.mType, std::move(data));
        NDArray tg(DType::Int32, {b.mCount});
        std::memcpy(tg.Data(), b.mTag.data(), b.mTag.size() * sizeof(std::int32_t));
        tags.push_back(std::move(tg));
    }
    mesh.AddCellData("su2:tag", std::move(tags));
    return mesh;
}

void write_su2(const std::string& rPath, const Mesh& rMesh) {
    std::ofstream os(rPath);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);

    const NDArray& points = rMesh.Points();
    const std::size_t dim = rMesh.PointDim();
    const std::size_t npoin = rMesh.NumPoints();

    os << "NDIME= " << dim << "\n";
    os << "NPOIN= " << npoin << "\n";
    {
        // Format point rows in parallel (snprintf per row, bytes unchanged),
        // then stream sequentially.
        std::vector<std::string> rows(npoin);
        parallel_for(npoin, [&](std::size_t i) {
            char buf[64];
            std::string& row = rows[i];
            for (std::size_t c = 0; c < dim; ++c) {
                std::snprintf(buf, sizeof(buf), "%.16e",
                              detail::read_double(points, i * dim + c));
                row += buf;
                row += (c + 1 == dim ? '\n' : ' ');
            }
        });
        for (const auto& row : rows)
            os << row;
    }

    std::vector<std::string> vtypes =
        (dim == 2) ? std::vector<std::string>{"triangle", "quad"}
                   : std::vector<std::string>{"tetra", "hexahedron", "wedge", "pyramid"};
    std::vector<std::string> btypes = (dim == 2) ? std::vector<std::string>{"line"}
                                                 : std::vector<std::string>{"triangle", "quad"};
    auto in = [](const std::vector<std::string>& v, const std::string& t) {
        return std::find(v.begin(), v.end(), t) != v.end();
    };

    // Volume cells.
    std::size_t nelem = 0;
    for (const auto cb : rMesh.CellRange())
        if (in(vtypes, cb.Type()))
            nelem += cb.NumCells();
    os << "NELEM= " << nelem << "\n";
    for (const auto cb : rMesh.CellRange()) {
        if (!in(vtypes, cb.Type()))
            continue;
        int st = meshio_to_su2(cb.Type());
        const NDArray& conn = cb.Conn();
        std::size_t k = conn.Shape().size() >= 2 ? conn.Shape()[1] : 1;
        for (std::size_t r = 0; r < cb.NumCells(); ++r) {
            os << st;
            for (std::size_t j = 0; j < k; ++j)
                os << " " << detail::read_int(conn, r * k + j);
            os << "\n";
        }
    }

    // Boundary markers from su2:tag (first int cell_data).
    std::string tag_key;
    for (const auto& name : rMesh.CellDataNames()) {
        if (rMesh.CellDataNumBlocks(name) == 0)
            continue;
        DType t = rMesh.CellData(name, 0).Dtype();
        if (t == DType::Int8 || t == DType::Int16 || t == DType::Int32 || t == DType::Int64 ||
            t == DType::UInt8 || t == DType::UInt16 || t == DType::UInt32 || t == DType::UInt64) {
            tag_key = name;
            break;
        }
    }

    // Collect unique tags (with total counts) over boundary cell blocks.
    std::map<std::int64_t, std::size_t> tag_counts;
    for (std::size_t bi = 0; bi < rMesh.NumCellBlocks(); ++bi) {
        const auto cb = rMesh.Cells(bi);
        if (!in(btypes, cb.Type()))
            continue;
        for (std::size_t r = 0; r < cb.NumCells(); ++r) {
            std::int64_t tg = 1;
            if (!tag_key.empty())
                tg = detail::read_int(rMesh.CellData(tag_key, bi), r);
            ++tag_counts[tg];
        }
    }

    os << "NMARK= " << tag_counts.size() << "\n";
    for (const auto& tc : tag_counts) {
        std::int64_t tag = tc.first;
        os << "MARKER_TAG= " << tag << "\n";
        os << "MARKER_ELEMS= " << tc.second << "\n";
        for (std::size_t bi = 0; bi < rMesh.NumCellBlocks(); ++bi) {
            const auto cb = rMesh.Cells(bi);
            if (!in(btypes, cb.Type()))
                continue;
            int st = meshio_to_su2(cb.Type());
            const NDArray& conn = cb.Conn();
            std::size_t k = conn.Shape().size() >= 2 ? conn.Shape()[1] : 1;
            for (std::size_t r = 0; r < cb.NumCells(); ++r) {
                std::int64_t tg = 1;
                if (!tag_key.empty())
                    tg = detail::read_int(rMesh.CellData(tag_key, bi), r);
                if (tg != tag)
                    continue;
                os << st;
                for (std::size_t j = 0; j < k; ++j)
                    os << " " << detail::read_int(conn, r * k + j);
                os << "\n";
            }
        }
    }
}

}  // namespace meshioplusplus
