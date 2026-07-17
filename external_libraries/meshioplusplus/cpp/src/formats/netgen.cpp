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
#include <iterator>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/netgen.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/types.hpp"

namespace meshioplusplus {

namespace {

// netgen cell node count -> meshio type, per topological dimension.
const std::unordered_map<int, std::string>& netgen_type(int dim) {
    static const std::unordered_map<int, std::string> d0 = {{1, "vertex"}};
    static const std::unordered_map<int, std::string> d1 = {{2, "line"}};
    static const std::unordered_map<int, std::string> d2 = {
        {3, "triangle"}, {6, "triangle6"}, {4, "quad"}, {8, "quad8"}};
    static const std::unordered_map<int, std::string> d3 = {
        {4, "tetra"},    {5, "pyramid"},    {6, "wedge"},    {8, "hexahedron"},
        {10, "tetra10"}, {13, "pyramid13"}, {15, "wedge15"}, {20, "hexahedron20"}};
    switch (dim) {
        case 0:
            return d0;
        case 1:
            return d1;
        case 2:
            return d2;
        default:
            return d3;
    }
}

// netgen -> meshio node permutation: meshio[i] = netgen[pmap[i]].
const std::unordered_map<std::string, std::vector<int>>& n2m_pmap() {
    static const std::unordered_map<std::string, std::vector<int>> m = {
        {"vertex", {0}},
        {"line", {0, 1}},
        {"triangle", {0, 1, 2}},
        {"triangle6", {0, 1, 2, 5, 3, 4}},
        {"quad", {0, 1, 2, 3}},
        {"quad8", {0, 1, 2, 3, 4, 7, 5, 6}},
        {"tetra", {0, 2, 1, 3}},
        {"tetra10", {0, 2, 1, 3, 5, 7, 4, 6, 9, 8}},
        {"pyramid", {0, 3, 2, 1, 4}},
        {"pyramid13", {0, 3, 2, 1, 4, 7, 6, 8, 5, 9, 12, 11, 10}},
        {"wedge", {0, 2, 1, 3, 5, 4}},
        {"wedge15", {0, 2, 1, 3, 5, 4, 7, 8, 6, 13, 14, 12, 9, 11, 10}},
        {"hexahedron", {0, 3, 2, 1, 4, 7, 6, 5}},
        {"hexahedron20", {0, 3, 2, 1, 4, 7, 6, 5, 10, 9, 11, 8, 16, 19, 18, 17, 14, 13, 15, 12}},
    };
    return m;
}

// meshio -> netgen node permutation (inverse of n2m_pmap).
const std::unordered_map<std::string, std::vector<int>>& m2n_pmap() {
    static const std::unordered_map<std::string, std::vector<int>> m = [] {
        std::unordered_map<std::string, std::vector<int>> out;
        for (const auto& kv : n2m_pmap()) {
            const auto& p = kv.second;
            std::vector<int> inv(p.size());
            for (std::size_t i = 0; i < p.size(); ++i)
                inv[p[i]] = static_cast<int>(i);
            out.emplace(kv.first, std::move(inv));
        }
        return out;
    }();
    return m;
}

int topo_dim(const std::string& rType) {
    auto it = topological_dimension().find(rType);
    return it == topological_dimension().end() ? -1 : it->second;
}

std::vector<std::string> split_ws(const std::string& rS) {
    std::vector<std::string> out;
    std::istringstream iss(rS);
    std::string tok;
    while (iss >> tok)
        out.push_back(tok);
    return out;
}

std::string strip(const std::string& rS) {
    std::size_t a = 0, b = rS.size();
    while (a < b && std::isspace(static_cast<unsigned char>(rS[a])))
        ++a;
    while (b > a && std::isspace(static_cast<unsigned char>(rS[b - 1])))
        --b;
    return rS.substr(a, b - a);
}

// Cursor over the file's lines, with comment/blank handling like the Python
// reader's _fast_forward_over_blank_lines.
struct LineCursor {
    std::vector<std::string> mLines;
    std::size_t mPos = 0;

    explicit LineCursor(std::istream& rIn) {
        std::string line;
        while (std::getline(rIn, line))
            mLines.push_back(line);
    }

    bool Eof() const { return mPos >= mLines.size(); }

    // Next non-blank, non-comment line (stripped). Sets is_eof when exhausted.
    std::string NextReal(bool& rIsEof) {
        while (mPos < mLines.size()) {
            std::string s = strip(mLines[mPos++]);
            if (!s.empty() && s[0] != '#') {
                rIsEof = false;
                return s;
            }
        }
        rIsEof = true;
        return "";
    }

    // Next line raw (stripped), used for count lines that directly follow a
    // keyword; skips any stray blank/comment lines defensively.
    std::string NextCount() {
        bool eof = false;
        return NextReal(eof);
    }
};

struct RawBlock {
    std::string mType;
    std::vector<std::vector<std::int64_t>> mRows;  // meshio node order, 0-based
    std::vector<std::int64_t> mIndex;
};

void read_cells(LineCursor& rC, const std::string& rSection, std::vector<RawBlock>& rBlocks) {
    int dim, pi0, i_index, fixed_nump = -1;
    if (rSection == "pointelements") {
        dim = 0;
        pi0 = 0;
        i_index = 1;
        fixed_nump = 1;
    } else if (rSection.rfind("edgesegments", 0) == 0) {
        dim = 1;
        pi0 = 2;
        i_index = 0;
        fixed_nump = 2;
    } else if (rSection.rfind("surfaceelements", 0) == 0) {
        dim = 2;
        pi0 = 5;
        i_index = 1;
    } else if (rSection == "volumeelements") {
        dim = 3;
        pi0 = 2;
        i_index = 0;
    } else {
        throw ReadError("Netgen: unknown cell section '" + rSection + "'");
    }

    std::int64_t num_cells = std::strtoll(rC.NextCount().c_str(), nullptr, 10);
    const auto& tmap = netgen_type(dim);

    for (std::int64_t k = 0; k < num_cells; ++k) {
        bool eof = false;
        std::string line = rC.NextReal(eof);
        if (eof)
            throw ReadError("Netgen: unexpected end of file in " + rSection);
        std::vector<std::string> data = split_ws(line);

        int nump = fixed_nump;
        if (dim == 2)
            nump = static_cast<int>(std::strtoll(data[4].c_str(), nullptr, 10));
        else if (dim == 3)
            nump = static_cast<int>(std::strtoll(data[1].c_str(), nullptr, 10));

        std::int64_t index = std::strtoll(data[i_index].c_str(), nullptr, 10);
        auto tit = tmap.find(nump);
        if (tit == tmap.end())
            throw ReadError("Netgen: unsupported element with " + std::to_string(nump) + " nodes");
        const std::string& type = tit->second;

        std::vector<std::int64_t> pi(nump);
        for (int j = 0; j < nump; ++j)
            pi[j] = std::strtoll(data[pi0 + j].c_str(), nullptr, 10);

        if (rBlocks.empty() || rBlocks.back().mType != type) {
            rBlocks.push_back(RawBlock{type, {}, {}});
        }
        rBlocks.back().mRows.push_back(std::move(pi));
        rBlocks.back().mIndex.push_back(index);
    }
}

}  // namespace

Mesh read_netgen(const std::string& rPath) {
    if (rPath.size() >= 7 && rPath.compare(rPath.size() - 7, 7, ".vol.gz") == 0)
        throw ReadError("Netgen: gzip container handled by Python fallback");

    std::ifstream in(rPath, std::ios::binary);
    if (!in)
        throw ReadError("Could not open file: " + rPath);
    LineCursor c(in);

    bool eof = false;
    std::string line = c.NextReal(eof);
    if (line != "mesh3d")
        throw ReadError("Not a valid Netgen mesh");

    int dimension = 3;
    std::vector<double> raw_points;  // flat, 3 per point
    std::int64_t num_points = 0;
    std::vector<RawBlock> blocks;

    while (true) {
        line = c.NextReal(eof);
        if (eof)
            break;
        if (line == "dimension") {
            dimension = static_cast<int>(std::strtoll(c.NextCount().c_str(), nullptr, 10));
        } else if (line == "geomtype") {
            c.NextCount();  // value; ignored
        } else if (line == "points") {
            num_points = std::strtoll(c.NextCount().c_str(), nullptr, 10);
            raw_points.resize(static_cast<std::size_t>(num_points) * 3, 0.0);
            for (std::int64_t i = 0; i < num_points; ++i) {
                std::string pl = c.NextReal(eof);
                if (eof)
                    throw ReadError("Netgen: unexpected EOF in points");
                std::vector<std::string> toks = split_ws(pl);
                for (int j = 0; j < 3 && j < static_cast<int>(toks.size()); ++j)
                    raw_points[i * 3 + j] = std::strtod(toks[j].c_str(), nullptr);
            }
        } else if (line == "pointelements" || line == "edgesegments" || line == "edgesegmentsgi" ||
                   line == "surfaceelements" || line == "surfaceelementsgi" ||
                   line == "surfaceelementsuv" || line == "volumeelements") {
            read_cells(c, line, blocks);
        } else if (line == "edgesegmentsgi2") {
            // Single-line variant (meshio's own output). The two-line variant
            // is signalled by a "surf1 surf2 p1 p2" header, handled below.
            read_cells(c, line, blocks);
        } else if (line == "endmesh") {
            break;
        } else {
            // identifications, materials/bcnames/cd*names, face_colours,
            // singular_*, the two-line edgesegmentsgi2 header, etc.
            throw ReadError("Netgen: token '" + line + "' handled by Python fallback");
        }
    }

    Mesh mesh;
    NDArray pts(DType::Float64,
                {static_cast<std::size_t>(num_points), static_cast<std::size_t>(dimension)});
    double* pp = pts.As<double>();
    for (std::int64_t i = 0; i < num_points; ++i)
        for (int j = 0; j < dimension; ++j)
            pp[i * dimension + j] = raw_points[i * 3 + j];
    mesh.AssignPoints(std::move(pts));

    std::vector<NDArray> index_blocks;
    for (auto& b : blocks) {
        const std::vector<int>& pmap = n2m_pmap().at(b.mType);
        std::size_t n = b.mRows.size();
        std::size_t k = pmap.size();
        NDArray data(DType::Int64, {n, k});
        std::int64_t* dp = data.As<std::int64_t>();
        for (std::size_t r = 0; r < n; ++r)
            for (std::size_t j = 0; j < k; ++j)
                dp[r * k + j] = b.mRows[r][pmap[j]] - 1;
        mesh.AddCellBlock(b.mType, std::move(data));

        NDArray idx(DType::Int64, {n});
        for (std::size_t r = 0; r < n; ++r)
            idx.As<std::int64_t>()[r] = b.mIndex[r];
        index_blocks.push_back(std::move(idx));
    }
    mesh.AddCellData("netgen:index", std::move(index_blocks));

    return mesh;
}

namespace {

void write_block(std::ostream& rOs, Mesh::CellView cb, const NDArray* pIndex) {
    if (cb.NumCells() == 0)
        return;
    int dim = topo_dim(cb.Type());
    const std::vector<int>& pmap = m2n_pmap().at(cb.Type());
    const int np = static_cast<int>(pmap.size());

    std::vector<std::int64_t> pre, post;
    int i_index = 0;
    if (dim == 0) {
        post = {1};
        i_index = 1;
    } else if (dim == 1) {
        pre = {1, 0};
        post = {-1, -1, 0, 0, 1, 0, 1, 0};
    } else if (dim == 2) {
        pre = {1, 1, 0, 0, np};
        i_index = 1;
    } else {  // dim == 3
        pre = {1, np};
    }

    const NDArray& conn = cb.Conn();
    const std::size_t n = cb.NumCells();
    for (std::size_t r = 0; r < n; ++r) {
        std::vector<std::int64_t> cols;
        cols.reserve(pre.size() + np + post.size());
        for (auto v : pre)
            cols.push_back(v);
        for (int j = 0; j < np; ++j)
            cols.push_back(detail::read_int(conn, r * np + pmap[j]) + 1);
        for (auto v : post)
            cols.push_back(v);
        if (pIndex)
            cols[i_index] = detail::read_int(*pIndex, r);

        for (std::size_t j = 0; j < cols.size(); ++j)
            rOs << cols[j] << (j + 1 == cols.size() ? '\n' : ' ');
    }
}

}  // namespace

void write_netgen(const std::string& rPath, const Mesh& rMesh, const std::string& rFloatFmt) {
    std::ofstream f(rPath, std::ios::binary);
    if (!f)
        throw WriteError("Could not open file for writing: " + rPath);

    const NDArray& points = rMesh.Points();
    const int dimension = points.Shape().size() >= 2 ? static_cast<int>(points.Shape()[1]) : 3;

    // Pick the single integer cell index, preferring "netgen:index".
    bool have_index = false;
    std::string index_key;
    if (rMesh.HasCellData("netgen:index")) {
        have_index = true;
        index_key = "netgen:index";
    } else {
        for (const auto& name : rMesh.CellDataNames()) {
            if (rMesh.CellDataNumBlocks(name) == 0)
                continue;
            DType t = rMesh.CellData(name, 0).Dtype();
            if (t != DType::Float32 && t != DType::Float64) {
                have_index = true;
                index_key = name;
                break;
            }
        }
    }
    auto index_for = [&](std::size_t ci) -> const NDArray* {
        if (!have_index || ci >= rMesh.CellDataNumBlocks(index_key))
            return nullptr;
        return &rMesh.CellData(index_key, ci);
    };

    std::int64_t per_dim[4] = {0, 0, 0, 0};
    for (const auto cb : rMesh.CellRange()) {
        int d = topo_dim(cb.Type());
        if (d >= 0 && d <= 3)
            per_dim[d] += static_cast<std::int64_t>(cb.NumCells());
    }

    f << "# Generated by meshio++ (C++ core)\n";
    f << "mesh3d\n\n";
    f << "dimension\n" << dimension << "\n\n";
    f << "geomtype\n0\n";

    f << "\n# surfnr    bcnr   domin  domout      np      p1      p2      p3\n";
    f << "surfaceelements\n" << per_dim[2] << "\n";
    for (std::size_t ci = 0; ci < rMesh.NumCellBlocks(); ++ci)
        if (topo_dim(rMesh.Cells(ci).Type()) == 2)
            write_block(f, rMesh.Cells(ci), index_for(ci));

    f << "\n#  matnr      np      p1      p2      p3      p4\n";
    f << "volumeelements\n" << per_dim[3] << "\n";
    for (std::size_t ci = 0; ci < rMesh.NumCellBlocks(); ++ci)
        if (topo_dim(rMesh.Cells(ci).Type()) == 3)
            write_block(f, rMesh.Cells(ci), index_for(ci));

    f << "\n# surfid  0   p1   p2   trignum1    trignum2   domin/surfnr1    "
         "domout/surfnr2   ednr1   dist1   ednr2   dist2\n";
    f << "edgesegmentsgi2\n" << per_dim[1] << "\n";
    for (std::size_t ci = 0; ci < rMesh.NumCellBlocks(); ++ci)
        if (topo_dim(rMesh.Cells(ci).Type()) == 1)
            write_block(f, rMesh.Cells(ci), index_for(ci));

    f << "\n#          X             Y             Z\n";
    f << "points\n" << rMesh.NumPoints() << "\n";
    std::string fmt = "%" + rFloatFmt;
    char buf[64];
    const std::size_t npts = rMesh.NumPoints();
    for (std::size_t i = 0; i < npts; ++i) {
        for (int j = 0; j < 3; ++j) {
            double v = (j < dimension) ? detail::read_double(points, i * dimension + j) : 0.0;
            std::snprintf(buf, sizeof(buf), fmt.c_str(), v);
            f << buf << (j == 2 ? '\n' : ' ');
        }
    }

    f << "\n#          pnum             index\n";
    f << "pointelements\n" << per_dim[0] << "\n";
    for (std::size_t ci = 0; ci < rMesh.NumCellBlocks(); ++ci)
        if (topo_dim(rMesh.Cells(ci).Type()) == 0)
            write_block(f, rMesh.Cells(ci), index_for(ci));

    f << "\nendmesh\n";
}

}  // namespace meshioplusplus
