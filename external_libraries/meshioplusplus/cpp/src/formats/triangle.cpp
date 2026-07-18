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
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/triangle.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"

namespace meshioplusplus {

namespace {

// Whitespace token stream over a whole file, '#' comments stripped to EOL.
struct TriangleTokens {
    std::vector<std::string> mToks;
    std::size_t mPos = 0;

    bool AtEnd() const { return mPos >= mToks.size(); }

    const std::string& Next(const char* pWhat) {
        if (AtEnd())
            throw ReadError(std::string("Triangle: unexpected end of file reading ") + pWhat);
        return mToks[mPos++];
    }

    std::int64_t NextInt(const char* pWhat) {
        const std::string& t = Next(pWhat);
        char* end = nullptr;
        const std::int64_t v = std::strtoll(t.c_str(), &end, 10);
        if (end == t.c_str())
            throw ReadError(std::string("Triangle: expected an integer for ") + pWhat);
        return v;
    }

    double NextDouble(const char* pWhat) {
        const std::string& t = Next(pWhat);
        char* end = nullptr;
        const double v = std::strtod(t.c_str(), &end);
        if (end == t.c_str())
            throw ReadError(std::string("Triangle: expected a number for ") + pWhat);
        return v;
    }
};

TriangleTokens triangle_tokenize(const std::string& rPath, bool& rOk) {
    TriangleTokens tokens;
    std::ifstream in(rPath, std::ios::binary);
    rOk = static_cast<bool>(in);
    if (!rOk)
        return tokens;
    std::string line;
    while (std::getline(in, line)) {
        const std::size_t hash = line.find('#');
        if (hash != std::string::npos)
            line.resize(hash);
        std::istringstream iss(line);
        std::string tok;
        while (iss >> tok)
            tokens.mToks.push_back(tok);
    }
    return tokens;
}

// Suffix classification: 0 = .node/.ele pair, 1 = .poly, -1 = unsupported.
int triangle_path_kind(const std::string& rPath, std::string& rStem) {
    const std::size_t dot = rPath.find_last_of('.');
    if (dot == std::string::npos)
        return -1;
    const std::string suffix = rPath.substr(dot);
    rStem = rPath.substr(0, dot);
    if (suffix == ".node" || suffix == ".ele")
        return 0;
    if (suffix == ".poly")
        return 1;
    return -1;
}

// One parsed vertex section (from a .node file or inline in a .poly).
struct TriangleNodes {
    std::int64_t mNumPoints = 0;
    std::int64_t mBase = 0;
    std::vector<double> mXY;                    // 2 per point
    std::vector<std::vector<double>> mAttrs;    // one column per attribute
    std::vector<std::vector<double>> mMarkers;  // one column per boundary marker
};

TriangleNodes triangle_read_node_section(TriangleTokens& rTokens) {
    TriangleNodes out;
    out.mNumPoints = rTokens.NextInt("vertex count");
    const std::int64_t dim = rTokens.NextInt("dimension");
    const std::int64_t nattr = rTokens.NextInt("attribute count");
    const std::int64_t nmark = rTokens.NextInt("marker count");
    if (dim != 2)
        throw ReadError("Triangle: need 2D points");
    if (out.mNumPoints < 0 || nattr < 0 || nmark < 0)
        throw ReadError("Triangle: malformed vertex header");

    out.mXY.resize(static_cast<std::size_t>(out.mNumPoints) * 2);
    out.mAttrs.assign(static_cast<std::size_t>(nattr),
                      std::vector<double>(static_cast<std::size_t>(out.mNumPoints)));
    out.mMarkers.assign(static_cast<std::size_t>(nmark),
                        std::vector<double>(static_cast<std::size_t>(out.mNumPoints)));

    for (std::int64_t i = 0; i < out.mNumPoints; ++i) {
        const std::int64_t idx = rTokens.NextInt("vertex index");
        if (i == 0)
            out.mBase = idx;
        if (idx != out.mBase + i)
            throw ReadError("Triangle: vertices not numbered consecutively");
        out.mXY[static_cast<std::size_t>(i) * 2] = rTokens.NextDouble("x coordinate");
        out.mXY[static_cast<std::size_t>(i) * 2 + 1] = rTokens.NextDouble("y coordinate");
        for (auto& col : out.mAttrs)
            col[static_cast<std::size_t>(i)] = rTokens.NextDouble("vertex attribute");
        for (auto& col : out.mMarkers)
            col[static_cast<std::size_t>(i)] = rTokens.NextDouble("boundary marker");
    }
    return out;
}

void triangle_apply_nodes(Mesh& rMesh, const TriangleNodes& rNodes) {
    const std::size_t n = static_cast<std::size_t>(rNodes.mNumPoints);
    NDArray pts = NDArray::Uninit(DType::Float64, {n, 2});
    if (n > 0)
        std::copy(rNodes.mXY.begin(), rNodes.mXY.end(), pts.As<double>());
    rMesh.AssignPoints(std::move(pts));

    for (std::size_t k = 0; k < rNodes.mAttrs.size(); ++k) {
        NDArray a = NDArray::Uninit(DType::Float64, {n});
        std::copy(rNodes.mAttrs[k].begin(), rNodes.mAttrs[k].end(), a.As<double>());
        rMesh.AddPointData("triangle:attr" + std::to_string(k + 1), std::move(a));
    }
    for (std::size_t k = 0; k < rNodes.mMarkers.size(); ++k) {
        std::string name = "triangle:ref" + (k == 0 ? std::string() : std::to_string(k + 1));
        NDArray a = NDArray::Uninit(DType::Float64, {n});
        std::copy(rNodes.mMarkers[k].begin(), rNodes.mMarkers[k].end(), a.As<double>());
        rMesh.AddPointData(std::move(name), std::move(a));
    }
}

Mesh triangle_read_node_ele(const std::string& rStem) {
    bool have_node = false;
    TriangleTokens node_tokens = triangle_tokenize(rStem + ".node", have_node);
    if (!have_node)
        throw ReadError("Triangle: could not open file: " + rStem + ".node");
    TriangleNodes nodes = triangle_read_node_section(node_tokens);

    Mesh mesh;
    triangle_apply_nodes(mesh, nodes);

    // The .ele sibling is optional: a lone .node file is a point cloud.
    bool have_ele = false;
    TriangleTokens ele_tokens = triangle_tokenize(rStem + ".ele", have_ele);
    if (!have_ele)
        return mesh;

    const std::int64_t ne = ele_tokens.NextInt("triangle count");
    const std::int64_t npc = ele_tokens.NextInt("nodes per triangle");
    const std::int64_t nattr = ele_tokens.NextInt("attribute count");
    if (npc != 3 && npc != 6)
        throw ReadError("Triangle: only 3- or 6-node triangles are supported");
    if (ne < 0 || nattr < 0)
        throw ReadError("Triangle: malformed .ele header");

    NDArray conn = NDArray::Uninit(DType::Int64,
                                   {static_cast<std::size_t>(ne), static_cast<std::size_t>(npc)});
    std::int64_t* cp = conn.As<std::int64_t>();
    std::vector<std::vector<double>> attrs(static_cast<std::size_t>(nattr),
                                           std::vector<double>(static_cast<std::size_t>(ne)));
    for (std::int64_t i = 0; i < ne; ++i) {
        ele_tokens.NextInt("triangle index");
        for (std::int64_t c = 0; c < npc; ++c) {
            const std::int64_t v = ele_tokens.NextInt("triangle connectivity") - nodes.mBase;
            if (v < 0 || v >= nodes.mNumPoints)
                throw ReadError("Triangle: connectivity index out of range");
            cp[i * npc + c] = v;
        }
        for (auto& col : attrs)
            col[static_cast<std::size_t>(i)] = ele_tokens.NextDouble("triangle attribute");
    }
    if (ne == 0)  // an empty .ele adds no block
        return mesh;
    mesh.AddCellBlock(npc == 3 ? "triangle" : "triangle6", std::move(conn));

    for (std::size_t k = 0; k < attrs.size(); ++k) {
        std::string name = "triangle:ref" + (k == 0 ? std::string() : std::to_string(k + 1));
        NDArray a = NDArray::Uninit(DType::Float64, {static_cast<std::size_t>(ne)});
        std::copy(attrs[k].begin(), attrs[k].end(), a.As<double>());
        std::vector<NDArray> blocks;
        blocks.push_back(std::move(a));
        mesh.AddCellData(std::move(name), std::move(blocks));
    }
    return mesh;
}

Mesh triangle_read_poly(const std::string& rPath, const std::string& rStem) {
    bool ok = false;
    TriangleTokens tokens = triangle_tokenize(rPath, ok);
    if (!ok)
        throw ReadError("Triangle: could not open file: " + rPath);

    // Vertex section: inline, or the sibling .node when the count is 0.
    TriangleNodes nodes;
    const std::size_t header_pos = tokens.mPos;
    const std::int64_t nv = tokens.NextInt("vertex count");
    if (nv == 0) {
        tokens.NextInt("dimension");
        tokens.NextInt("attribute count");
        tokens.NextInt("marker count");
        bool have_node = false;
        TriangleTokens node_tokens = triangle_tokenize(rStem + ".node", have_node);
        if (!have_node)
            throw ReadError("Triangle: .poly refers to a missing sibling .node file");
        nodes = triangle_read_node_section(node_tokens);
    } else {
        tokens.mPos = header_pos;
        nodes = triangle_read_node_section(tokens);
    }

    Mesh mesh;
    triangle_apply_nodes(mesh, nodes);

    // Segment section -> one "line" cell block (+ optional marker cell_data).
    const std::int64_t ns = tokens.NextInt("segment count");
    const std::int64_t nmark = tokens.NextInt("segment marker count");
    if (ns < 0 || nmark < 0 || nmark > 1)
        throw ReadError("Triangle: malformed segment header");
    NDArray conn = NDArray::Uninit(DType::Int64, {static_cast<std::size_t>(ns), 2});
    std::int64_t* cp = conn.As<std::int64_t>();
    std::vector<std::int64_t> markers(static_cast<std::size_t>(nmark == 1 ? ns : 0));
    for (std::int64_t i = 0; i < ns; ++i) {
        tokens.NextInt("segment index");
        for (int c = 0; c < 2; ++c) {
            const std::int64_t v = tokens.NextInt("segment endpoint") - nodes.mBase;
            if (v < 0 || v >= nodes.mNumPoints)
                throw ReadError("Triangle: segment endpoint out of range");
            cp[i * 2 + c] = v;
        }
        if (nmark == 1)
            markers[static_cast<std::size_t>(i)] = tokens.NextInt("segment marker");
    }
    mesh.AddCellBlock("line", std::move(conn));
    if (nmark == 1) {
        NDArray a = NDArray::Uninit(DType::Int64, {static_cast<std::size_t>(ns)});
        std::copy(markers.begin(), markers.end(), a.As<std::int64_t>());
        std::vector<NDArray> blocks;
        blocks.push_back(std::move(a));
        mesh.AddCellData("triangle:ref", std::move(blocks));
    }

    // Holes (and optional regional attributes) are not representable — skip.
    if (!tokens.AtEnd()) {
        const std::int64_t nh = tokens.NextInt("hole count");
        if (nh > 0)
            log::warn("Triangle: skipping {} hole(s) in {}", nh, rPath);
        for (std::int64_t i = 0; i < nh * 3; ++i)
            tokens.NextDouble("hole entry");
    }
    if (!tokens.AtEnd()) {
        const std::int64_t nr = tokens.NextInt("region count");
        if (nr > 0)
            log::warn("Triangle: skipping {} regional attribute(s) in {}", nr, rPath);
    }
    return mesh;
}

// Write a marker/ref value: integral values as integers, else %.16e.
void triangle_write_value(std::ostream& rOs, double v) {
    const double r = std::nearbyint(v);
    if (v == r && std::fabs(v) < 9.2e18) {
        rOs << static_cast<std::int64_t>(r);
    } else {
        char buf[40];
        std::snprintf(buf, sizeof(buf), "%.16e", v);
        rOs << buf;
    }
}

// Split point_data keys into attribute columns and at most one marker
// column (the first ":ref" key, or else the first key), mirroring tetgen.
void triangle_split_point_keys(const Mesh& rMesh, std::vector<std::string>& rAttrKeys,
                               std::vector<std::string>& rRefKeys) {
    rAttrKeys = rMesh.PointDataNames();  // sorted: deterministic column order
    if (rAttrKeys.empty())
        return;
    for (const auto& k : rAttrKeys)
        if (k.find(":ref") != std::string::npos) {
            rRefKeys.push_back(k);
            break;
        }
    if (!rRefKeys.empty()) {
        rAttrKeys.erase(std::remove(rAttrKeys.begin(), rAttrKeys.end(), rRefKeys[0]),
                        rAttrKeys.end());
    } else {
        rRefKeys.push_back(rAttrKeys.front());
        rAttrKeys.erase(rAttrKeys.begin());
    }
}

void triangle_write_node_rows(std::ostream& rOs, const Mesh& rMesh,
                              const std::vector<std::string>& rAttrKeys,
                              const std::vector<std::string>& rRefKeys) {
    const NDArray& points = rMesh.Points();
    const std::int64_t np = static_cast<std::int64_t>(rMesh.NumPoints());
    char fbuf[40];
    for (std::int64_t i = 0; i < np; ++i) {
        rOs << i;
        for (int c = 0; c < 2; ++c) {
            std::snprintf(fbuf, sizeof(fbuf), "%.16e", detail::read_double(points, i * 2 + c));
            rOs << " " << fbuf;
        }
        for (const auto& k : rAttrKeys) {
            std::snprintf(fbuf, sizeof(fbuf), "%.16e", detail::read_double(rMesh.PointData(k), i));
            rOs << " " << fbuf;
        }
        for (const auto& k : rRefKeys) {
            rOs << " ";
            triangle_write_value(rOs, detail::read_double(rMesh.PointData(k), i));
        }
        rOs << "\n";
    }
}

void triangle_write_node_ele(const std::string& rStem, const Mesh& rMesh) {
    std::vector<std::string> attr_keys, ref_keys;
    triangle_split_point_keys(rMesh, attr_keys, ref_keys);

    {
        std::ofstream fh(rStem + ".node", std::ios::binary);
        if (!fh)
            throw WriteError("Could not open file for writing: " + rStem + ".node");
        fh << "# This file was created by meshio++ (C++ core)\n";
        fh << rMesh.NumPoints() << " 2 " << attr_keys.size() << " " << ref_keys.size() << "\n";
        triangle_write_node_rows(fh, rMesh, attr_keys, ref_keys);
    }

    // Collect the triangle blocks (must all share one type).
    std::string tri_type;
    std::size_t ne = 0;
    for (const auto cb : rMesh.CellRange()) {
        if (cb.Type() != "triangle" && cb.Type() != "triangle6")
            continue;
        if (!tri_type.empty() && tri_type != cb.Type())
            throw WriteError("Triangle: cannot mix triangle and triangle6 blocks");
        tri_type = cb.Type();
        ne += cb.NumCells();
    }

    std::vector<std::string> cell_attr_keys = rMesh.CellDataNames();
    if (!cell_attr_keys.empty()) {
        std::string ref;
        for (const auto& k : cell_attr_keys)
            if (k.find(":ref") != std::string::npos) {
                ref = k;
                break;
            }
        if (!ref.empty()) {
            cell_attr_keys.erase(std::remove(cell_attr_keys.begin(), cell_attr_keys.end(), ref),
                                 cell_attr_keys.end());
            cell_attr_keys.insert(cell_attr_keys.begin(), ref);
        }
    }

    std::ofstream fh(rStem + ".ele", std::ios::binary);
    if (!fh)
        throw WriteError("Could not open file for writing: " + rStem + ".ele");
    fh << "# This file was created by meshio++ (C++ core)\n";
    const std::size_t npc = tri_type == "triangle6" ? 6 : 3;
    fh << ne << " " << npc << " " << cell_attr_keys.size() << "\n";
    std::int64_t id = 0;
    for (std::size_t ci = 0; ci < rMesh.NumCellBlocks(); ++ci) {
        const auto cb = rMesh.Cells(ci);
        if (cb.Type() != tri_type)
            continue;
        const NDArray& conn = cb.Conn();
        const std::int64_t n = static_cast<std::int64_t>(cb.NumCells());
        for (std::int64_t i = 0; i < n; ++i) {
            fh << id++;
            for (std::size_t c = 0; c < npc; ++c)
                fh << " " << detail::read_int(conn, i * static_cast<std::int64_t>(npc) + c);
            for (const auto& k : cell_attr_keys) {
                fh << " ";
                if (ci < rMesh.CellDataNumBlocks(k))
                    triangle_write_value(fh, detail::read_double(rMesh.CellData(k, ci), i));
                else
                    fh << "0";
            }
            fh << "\n";
        }
    }
}

void triangle_write_poly(const std::string& rPath, const Mesh& rMesh) {
    std::vector<std::string> attr_keys, ref_keys;
    triangle_split_point_keys(rMesh, attr_keys, ref_keys);

    // Segment markers: the first ":ref" cell_data key, when present.
    std::string seg_ref;
    for (const auto& k : rMesh.CellDataNames())
        if (k.find(":ref") != std::string::npos) {
            seg_ref = k;
            break;
        }

    std::ofstream fh(rPath, std::ios::binary);
    if (!fh)
        throw WriteError("Could not open file for writing: " + rPath);
    fh << "# This file was created by meshio++ (C++ core)\n";
    fh << rMesh.NumPoints() << " 2 " << attr_keys.size() << " " << ref_keys.size() << "\n";
    triangle_write_node_rows(fh, rMesh, attr_keys, ref_keys);

    std::size_t ns = 0;
    for (const auto cb : rMesh.CellRange())
        if (cb.Type() == "line")
            ns += cb.NumCells();
    fh << ns << " " << (seg_ref.empty() ? 0 : 1) << "\n";
    std::int64_t id = 0;
    for (std::size_t ci = 0; ci < rMesh.NumCellBlocks(); ++ci) {
        const auto cb = rMesh.Cells(ci);
        if (cb.Type() != "line")
            continue;
        const NDArray& conn = cb.Conn();
        const std::int64_t n = static_cast<std::int64_t>(cb.NumCells());
        for (std::int64_t i = 0; i < n; ++i) {
            fh << id++ << " " << detail::read_int(conn, i * 2) << " "
               << detail::read_int(conn, i * 2 + 1);
            if (!seg_ref.empty()) {
                fh << " ";
                if (ci < rMesh.CellDataNumBlocks(seg_ref))
                    triangle_write_value(fh, detail::read_double(rMesh.CellData(seg_ref, ci), i));
                else
                    fh << "0";
            }
            fh << "\n";
        }
    }
    fh << "0\n";  // holes
}

}  // namespace

Mesh read_triangle(const std::string& rPath) {
    std::string stem;
    const int kind = triangle_path_kind(rPath, stem);
    if (kind < 0)
        throw ReadError("Triangle: expected a .node, .ele, or .poly file");
    if (kind == 1)
        return triangle_read_poly(rPath, stem);
    return triangle_read_node_ele(stem);
}

void write_triangle(const std::string& rPath, const Mesh& rMesh) {
    std::string stem;
    const int kind = triangle_path_kind(rPath, stem);
    if (kind < 0)
        throw WriteError("Triangle: must specify a .node, .ele, or .poly file");
    if (rMesh.PointDim() != 2)
        throw WriteError("Triangle: can only write 2D points");
    if (kind == 1)
        triangle_write_poly(rPath, rMesh);
    else
        triangle_write_node_ele(stem, rMesh);
}

}  // namespace meshioplusplus
