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
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/avsucd.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus {

namespace {

const std::unordered_map<std::string, std::string>& meshio_to_avsucd_type() {
    static const std::unordered_map<std::string, std::string> m = {
        {"vertex", "pt"}, {"line", "line"},   {"triangle", "tri"}, {"quad", "quad"},
        {"tetra", "tet"}, {"pyramid", "pyr"}, {"wedge", "prism"},  {"hexahedron", "hex"},
    };
    return m;
}
const std::unordered_map<std::string, std::string>& avsucd_to_meshio_type() {
    static const std::unordered_map<std::string, std::string> m = [] {
        std::unordered_map<std::string, std::string> r;
        for (const auto& kv : meshio_to_avsucd_type())
            r[kv.second] = kv.first;
        return r;
    }();
    return m;
}
// meshio -> avsucd column order (empty = identity).
const std::vector<int>& meshio_to_avsucd_order(const std::string& rT) {
    static const std::unordered_map<std::string, std::vector<int>> m = {
        {"tetra", {0, 1, 3, 2}},
        {"pyramid", {4, 0, 1, 2, 3}},
        {"wedge", {3, 4, 5, 0, 1, 2}},
        {"hexahedron", {4, 5, 6, 7, 0, 1, 2, 3}},
    };
    static const std::vector<int> empty;
    auto it = m.find(rT);
    return it == m.end() ? empty : it->second;
}
const std::vector<int>& avsucd_to_meshio_order(const std::string& rT) {
    static const std::unordered_map<std::string, std::vector<int>> m = {
        {"tetra", {0, 1, 3, 2}},
        {"pyramid", {1, 2, 3, 4, 0}},
        {"wedge", {3, 4, 5, 0, 1, 2}},
        {"hexahedron", {4, 5, 6, 7, 0, 1, 2, 3}},
    };
    static const std::vector<int> empty;
    auto it = m.find(rT);
    return it == m.end() ? empty : it->second;
}

bool is_int_dtype(DType t) {
    return t == DType::Int8 || t == DType::Int16 || t == DType::Int32 || t == DType::Int64 ||
           t == DType::UInt8 || t == DType::UInt16 || t == DType::UInt32 || t == DType::UInt64;
}

std::vector<std::string> avsucd_tokens(const std::string& rS) {
    std::vector<std::string> out;
    std::istringstream iss(rS);
    std::string t;
    while (iss >> t)
        out.push_back(t);
    return out;
}

}  // namespace

Mesh read_avsucd(const std::string& rPath) {
    std::ifstream in(rPath);
    if (!in)
        throw ReadError("Could not open file: " + rPath);
    std::vector<std::string> lines;
    std::string l;
    while (std::getline(in, l)) {
        if (!l.empty() && l.back() == '\r')
            l.pop_back();
        std::string t = l;
        std::size_t b = t.find_first_not_of(" \t");
        if (b == std::string::npos)
            continue;  // blank
        if (t[b] == '#')
            continue;  // comment
        lines.push_back(l);
    }
    std::size_t li = 0;

    auto hdr = avsucd_tokens(lines.at(li++));
    long long num_nodes = std::stoll(hdr[0]);
    long long num_cells = std::stoll(hdr[1]);
    long long num_node_data = std::stoll(hdr[2]);
    long long num_cell_data = std::stoll(hdr[3]);

    Mesh mesh;
    std::unordered_map<std::int64_t, std::int64_t> point_ids;
    NDArray pts(DType::Float64, {static_cast<std::size_t>(num_nodes), 3});
    double* pp = pts.As<double>();
    for (long long i = 0; i < num_nodes; ++i) {
        auto t = avsucd_tokens(lines.at(li++));
        point_ids[std::strtoll(t[0].c_str(), nullptr, 10)] = i;
        for (int c = 0; c < 3; ++c)
            pp[i * 3 + c] = std::strtod(t[1 + c].c_str(), nullptr);
    }
    mesh.AssignPoints(std::move(pts));

    // Cells, grouped by consecutive type.
    std::unordered_map<std::int64_t, std::int64_t> cell_ids;
    struct Blk {
        std::string mType;
        int mN;
        std::vector<std::int64_t> mConn;
        std::vector<std::int64_t> mMat;
        std::size_t mCount = 0;
    };
    std::vector<Blk> blocks;
    for (long long c = 0; c < num_cells; ++c) {
        auto t = avsucd_tokens(lines.at(li++));
        std::int64_t cid = std::strtoll(t[0].c_str(), nullptr, 10);
        std::int64_t mat = std::strtoll(t[1].c_str(), nullptr, 10);
        auto it = avsucd_to_meshio_type().find(t[2]);
        if (it == avsucd_to_meshio_type().end())
            throw ReadError("AVS-UCD: unknown cell type '" + t[2] + "'");
        const std::string& mtype = it->second;
        int n = static_cast<int>(t.size()) - 3;
        if (blocks.empty() || blocks.back().mType != mtype) {
            Blk b;
            b.mType = mtype;
            b.mN = n;
            blocks.push_back(std::move(b));
        }
        Blk& blk = blocks.back();
        for (int j = 0; j < n; ++j)
            blk.mConn.push_back(point_ids.at(std::strtoll(t[3 + j].c_str(), nullptr, 10)));
        blk.mMat.push_back(mat);
        cell_ids[cid] = c;
        ++blk.mCount;
    }

    std::vector<NDArray> material_blocks;
    for (auto& blk : blocks) {
        const std::vector<int>& perm = avsucd_to_meshio_order(blk.mType);
        NDArray data(DType::Int64, {blk.mCount, static_cast<std::size_t>(blk.mN)});
        std::int64_t* dp = data.As<std::int64_t>();
        for (std::size_t r = 0; r < blk.mCount; ++r)
            for (int j = 0; j < blk.mN; ++j) {
                int src = perm.empty() ? j : perm[j];
                dp[r * blk.mN + j] = blk.mConn[r * blk.mN + src];
            }
        mesh.AddCellBlock(blk.mType, std::move(data));
        NDArray m(DType::Int64, {blk.mCount});
        for (std::size_t r = 0; r < blk.mCount; ++r)
            m.As<std::int64_t>()[r] = blk.mMat[r];
        material_blocks.push_back(std::move(m));
    }
    mesh.AddCellData("avsucd:material", std::move(material_blocks));

    // Reads a data section into name -> (num_entities, size) arrays.
    auto read_data = [&](long long num_entities,
                         const std::unordered_map<std::int64_t, std::int64_t>& ids,
                         std::vector<std::string>& names, std::vector<NDArray>& arrays) {
        auto h = avsucd_tokens(lines.at(li++));
        int narr = std::stoi(h[0]);
        std::vector<int> sizes(narr);
        for (int i = 0; i < narr; ++i)
            sizes[i] = std::stoi(h[1 + i]);
        for (int i = 0; i < narr; ++i) {
            std::string lbl = lines.at(li++);
            std::size_t comma = lbl.find(',');
            std::string name = (comma == std::string::npos) ? lbl : lbl.substr(0, comma);
            // strip + replace spaces with underscore
            std::string clean;
            for (char ch : name) {
                if (ch == ' ')
                    clean += '_';
                else if (!std::isspace(static_cast<unsigned char>(ch)))
                    clean += ch;
            }
            names.push_back(clean);
            arrays.emplace_back(DType::Float64,
                                sizes[i] == 1 ? std::vector<std::size_t>{(std::size_t)num_entities}
                                              : std::vector<std::size_t>{(std::size_t)num_entities,
                                                                         (std::size_t)sizes[i]});
        }
        for (long long e = 0; e < num_entities; ++e) {
            auto t = avsucd_tokens(lines.at(li++));
            std::int64_t eid = ids.at(std::strtoll(t[0].c_str(), nullptr, 10));
            std::size_t j = 1;
            for (int i = 0; i < narr; ++i) {
                for (int c = 0; c < sizes[i]; ++c)
                    arrays[i].As<double>()[eid * sizes[i] + c] =
                        std::strtod(t[j++].c_str(), nullptr);
            }
        }
    };

    if (num_node_data > 0) {
        std::vector<std::string> names;
        std::vector<NDArray> arrays;
        read_data(num_nodes, point_ids, names, arrays);
        for (std::size_t i = 0; i < names.size(); ++i)
            mesh.AddPointData(names[i], std::move(arrays[i]));
    }
    if (num_cell_data > 0) {
        std::vector<std::string> names;
        std::vector<NDArray> arrays;
        read_data(num_cells, cell_ids, names, arrays);
        // split each into per-block arrays
        for (std::size_t i = 0; i < names.size(); ++i) {
            const NDArray& a = arrays[i];
            std::size_t nc = a.Shape().size() >= 2 ? a.Shape()[1] : 1;
            std::size_t isz = dtype_size(a.Dtype());
            std::vector<NDArray> per_block;
            std::size_t offset = 0;
            for (auto& blk : blocks) {
                std::vector<std::size_t> shp = (nc == 1) ? std::vector<std::size_t>{blk.mCount}
                                                         : std::vector<std::size_t>{blk.mCount, nc};
                NDArray out(DType::Float64, shp);
                std::memcpy(out.Data(), a.Data() + offset * nc * isz, blk.mCount * nc * isz);
                per_block.push_back(std::move(out));
                offset += blk.mCount;
            }
            mesh.AddCellData(names[i], std::move(per_block));
        }
    }

    return mesh;
}

void write_avsucd(const std::string& rPath, const Mesh& rMesh) {
    std::ofstream os(rPath);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);

    const std::size_t num_nodes = rMesh.NumPoints();
    const std::size_t dim = rMesh.PointDim();
    std::size_t num_cells = 0;
    for (const auto cb : rMesh.CellRange())
        num_cells += cb.NumCells();

    // Material = first int cell_data array (avsucd:material if present).
    std::string mat_key;
    for (const auto& name : rMesh.CellDataNames()) {
        if (rMesh.CellDataNumBlocks(name) > 0 && is_int_dtype(rMesh.CellData(name, 0).Dtype())) {
            mat_key = name;
            break;
        }
    }

    // Node/cell data breakdowns (excluding material).
    std::vector<std::pair<std::string, const NDArray*>> ndata;
    std::vector<int> nsize;
    std::size_t nsum = 0;
    for (const auto& name : rMesh.PointDataNames()) {
        const NDArray& d = rMesh.PointData(name);
        int sz = d.Shape().size() >= 2 ? static_cast<int>(d.Shape()[1]) : 1;
        ndata.emplace_back(name, &d);
        nsize.push_back(sz);
        nsum += sz;
    }
    std::vector<std::string> cdata;
    std::vector<int> csize;
    std::size_t csum = 0;
    for (const auto& name : rMesh.CellDataNames()) {
        if (name == mat_key)
            continue;
        int sz = (rMesh.CellDataNumBlocks(name) > 0 &&
                  rMesh.CellData(name, 0).Shape().size() >= 2)
                     ? static_cast<int>(rMesh.CellData(name, 0).Shape()[1])
                     : 1;
        cdata.push_back(name);
        csize.push_back(sz);
        csum += sz;
    }

    os << "# Written by meshio++ (C++ core)\n";
    os << num_nodes << " " << num_cells << " " << nsum << " " << csum << " 0\n";

    const NDArray& points = rMesh.Points();
    char buf[48];
    for (std::size_t i = 0; i < num_nodes; ++i) {
        os << (i + 1);
        for (int c = 0; c < 3; ++c) {
            double v = (std::size_t(c) < dim) ? detail::read_double(points, i * dim + c) : 0.0;
            std::snprintf(buf, sizeof(buf), " %.17g", v);
            os << buf;
        }
        os << "\n";
    }

    // Cells
    std::size_t gi = 0;
    for (std::size_t bi = 0; bi < rMesh.NumCellBlocks(); ++bi) {
        const auto cb = rMesh.Cells(bi);
        auto it = meshio_to_avsucd_type().find(cb.Type());
        if (it == meshio_to_avsucd_type().end())
            throw WriteError("AVS-UCD writer: unsupported cell type " + cb.Type());
        const std::vector<int>& perm = meshio_to_avsucd_order(cb.Type());
        const NDArray& conn = cb.Conn();
        std::size_t n = conn.Shape().size() >= 2 ? conn.Shape()[1] : 1;
        const NDArray* mat = nullptr;
        if (!mat_key.empty())
            mat = &rMesh.CellData(mat_key, bi);
        for (std::size_t r = 0; r < cb.NumCells(); ++r) {
            std::int64_t m = mat ? detail::read_int(*mat, r) : 0;
            os << (gi + 1) << " " << m << " " << it->second;
            for (std::size_t j = 0; j < n; ++j) {
                std::size_t src = perm.empty() ? j : static_cast<std::size_t>(perm[j]);
                os << " " << (detail::read_int(conn, r * n + src) + 1);
            }
            os << "\n";
            ++gi;
        }
    }

    // Node data section.
    auto write_section = [&](std::size_t num_entities, const std::vector<int>& sizes,
                             const std::vector<std::string>& names,
                             auto value_at /* (idx, comp) -> double */) {
        os << sizes.size();
        for (int s : sizes)
            os << " " << s;
        os << "\n";
        for (const auto& nm : names)
            os << nm << ", real\n";
        for (std::size_t e = 0; e < num_entities; ++e) {
            os << (e + 1);
            for (std::size_t a = 0; a < sizes.size(); ++a)
                for (int c = 0; c < sizes[a]; ++c) {
                    std::snprintf(buf, sizeof(buf), " %.14e", value_at(a, e, c));
                    os << buf;
                }
            os << "\n";
        }
    };

    if (nsum > 0) {
        std::vector<std::string> names;
        for (auto& p : ndata)
            names.push_back(p.first);
        write_section(num_nodes, nsize, names, [&](std::size_t a, std::size_t e, int c) {
            const NDArray* arr = ndata[a].second;
            std::size_t sz = static_cast<std::size_t>(nsize[a]);
            return detail::read_double(*arr, e * sz + c);
        });
    }
    if (csum > 0) {
        // Flatten each cell-data name across blocks for global indexing.
        write_section(num_cells, csize, cdata, [&](std::size_t a, std::size_t e, int c) {
            const std::string& name = cdata[a];
            std::size_t sz = static_cast<std::size_t>(csize[a]);
            std::size_t idx = e;
            const std::size_t nblocks = rMesh.CellDataNumBlocks(name);
            for (std::size_t bi = 0; bi < nblocks; ++bi) {
                const NDArray& blk = rMesh.CellData(name, bi);
                std::size_t bcount = blk.Shape().empty() ? 0 : blk.Shape()[0];
                if (idx < bcount)
                    return detail::read_double(blk, idx * sz + c);
                idx -= bcount;
            }
            return 0.0;
        });
    }
}

}  // namespace meshioplusplus
