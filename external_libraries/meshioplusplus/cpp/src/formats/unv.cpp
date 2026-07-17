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
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/unv.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"

namespace meshioplusplus {

namespace {

// Salome/UNV parabolic node order -> meshio position (0-based).
const std::vector<int>* nd_perm(const std::string& rT) {
    static const std::unordered_map<std::string, std::vector<int>> m = {
        {"line3", {0, 2, 1}},
        {"triangle6", {0, 3, 1, 4, 2, 5}},
        {"quad8", {0, 4, 1, 5, 2, 6, 3, 7}},
        {"tetra10", {0, 4, 1, 5, 2, 6, 7, 8, 9, 3}},
        {"wedge15", {0, 6, 1, 7, 2, 8, 9, 10, 11, 3, 12, 4, 13, 5, 14}},
        {"hexahedron20", {0, 8, 1, 9, 2, 10, 3, 11, 12, 13, 14, 15, 4, 16, 5, 17, 6, 18, 7, 19}}};
    auto it = m.find(rT);
    return it == m.end() ? nullptr : &it->second;
}

std::string unv_type(int fedesc) {
    static const std::unordered_map<int, std::string> m = {
        {11, "line"},      {21, "line"},        {22, "line3"},        {24, "line3"},
        {41, "triangle"},  {81, "triangle"},    {91, "triangle"},     {42, "triangle6"},
        {82, "triangle6"}, {92, "triangle6"},   {44, "quad"},         {84, "quad"},
        {94, "quad"},      {122, "quad"},       {45, "quad8"},        {85, "quad8"},
        {95, "quad8"},     {111, "tetra"},      {118, "tetra10"},     {112, "wedge"},
        {113, "wedge15"},  {115, "hexahedron"}, {116, "hexahedron20"}};
    auto it = m.find(fedesc);
    return it == m.end() ? std::string() : it->second;
}

// meshio type -> (descriptor, is_beam)
bool meshio_descriptor(const std::string& rT, int& rDesc, bool& rBeam) {
    static const std::unordered_map<std::string, std::pair<int, bool>> m = {
        {"line", {21, true}},       {"line3", {24, true}},        {"triangle", {91, false}},
        {"triangle6", {92, false}}, {"quad", {94, false}},        {"quad8", {95, false}},
        {"tetra", {111, false}},    {"tetra10", {118, false}},    {"wedge", {112, false}},
        {"wedge15", {113, false}},  {"hexahedron", {115, false}}, {"hexahedron20", {116, false}}};
    auto it = m.find(rT);
    if (it == m.end())
        return false;
    rDesc = it->second.first;
    rBeam = it->second.second;
    return true;
}

bool is_beam(int fedesc) {
    return fedesc == 11 || fedesc == 21 || fedesc == 22 || fedesc == 24;
}

std::vector<std::string> tokens(const std::string& rS) {
    std::vector<std::string> out;
    std::istringstream iss(rS);
    std::string t;
    while (iss >> t)
        out.push_back(t);
    return out;
}

double parse_coord(std::string s) {
    for (char& c : s)
        if (c == 'D' || c == 'd')
            c = 'E';
    return std::strtod(s.c_str(), nullptr);
}

// Component count -> UNV data-characteristic code (1 scalar, 2 3-vector,
// 4 symmetric tensor, 5 general tensor); 0 when unrecognized.
int field_char(int ncomp) {
    switch (ncomp) {
        case 1:
            return 1;
        case 3:
            return 2;
        case 6:
            return 4;
        case 9:
            return 5;
        default:
            return 0;
    }
}

// A results field parsed from dataset 2414 / 55 / 56 / 57.
struct Field {
    int mLocation = 1;  // 1 = data at nodes, 2 = data on elements
    int mNcomp = 1;
    std::string mName;
    std::unordered_map<std::int64_t, std::vector<double>> mValues;  // entity label -> values
};

// Parse a field dataset body (lines[start, end)).  Returns false for
// unsupported datasets (complex data or nodes-on-elements location).
bool parse_field(int ds, const std::vector<std::string>& lines, std::size_t start, std::size_t end,
                 Field& rField) {
    auto row = [&](std::size_t i) {
        return i < end ? tokens(lines[i]) : std::vector<std::string>{};
    };
    std::size_t header;  // first index of the per-entity data
    int data_type, ndv;
    if (ds == 2414) {
        if (start + 13 > end)
            return false;
        rField.mName = lines[start + 1];
        auto loc = row(start + 2);
        rField.mLocation = loc.empty() ? 1 : std::atoi(loc[0].c_str());
        auto r9 = row(start + 8);
        if (r9.size() < 6)
            return false;
        data_type = std::atoi(r9[4].c_str());
        ndv = std::atoi(r9[5].c_str());
        header = start + 13;
    } else {
        if (start + 10 > end)
            return false;
        rField.mName = lines[start];
        auto r6 = row(start + 5);
        if (r6.size() < 6)
            return false;
        data_type = std::atoi(r6[4].c_str());
        ndv = std::atoi(r6[5].c_str());
        rField.mLocation = (ds == 57) ? 2 : (ds == 56 ? 3 : 1);
        header = start + 10;
    }
    // trim trailing whitespace from the name
    std::size_t last = rField.mName.find_last_not_of(" \t\r");
    rField.mName = last == std::string::npos ? std::string() : rField.mName.substr(0, last + 1);
    std::size_t first = rField.mName.find_first_not_of(" \t");
    if (first != std::string::npos)
        rField.mName = rField.mName.substr(first);

    if (data_type != 2 && data_type != 4)
        return false;  // complex data unsupported
    if (rField.mLocation == 3 || ndv <= 0)
        return false;  // nodes-on-elements averaging unsupported
    rField.mNcomp = ndv;

    std::size_t k = header;
    while (k < end) {
        auto rec = tokens(lines[k]);
        if (rec.empty()) {
            ++k;
            continue;
        }
        std::int64_t label = std::strtoll(rec[0].c_str(), nullptr, 10);
        ++k;
        std::vector<double> vals;
        while (static_cast<int>(vals.size()) < ndv && k < end) {
            for (const auto& v : tokens(lines[k]))
                vals.push_back(parse_coord(v));
            ++k;
        }
        vals.resize(ndv);
        rField.mValues[label] = std::move(vals);
    }
    return true;
}

}  // namespace

Mesh read_unv(const std::string& rPath, UnvInfo& rInfo) {
    std::ifstream in(rPath, std::ios::binary);
    if (!in)
        throw ReadError("Could not open file: " + rPath);
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(in, line))
        lines.push_back(line);

    std::vector<std::vector<double>> points;
    std::unordered_map<std::int64_t, std::int64_t> label_to_index;

    struct Group {
        std::string mType;
        std::vector<std::vector<std::int64_t>> mRows;
        std::vector<std::int64_t> mPid;
    };
    std::vector<Group> groups;
    std::unordered_map<std::string, std::size_t> group_index;
    // element label -> (block index, local index within block)
    std::unordered_map<std::int64_t, std::pair<std::size_t, std::size_t>> elem_label_to_ref;
    std::vector<Field> fields;
    std::unordered_set<std::string> used_keys;
    // raw permanent groups: (name, [(entity_type, tag)]) resolved after all
    // node/element datasets have been read.
    std::vector<std::pair<std::string, std::vector<std::pair<int, std::int64_t>>>> raw_groups;
    std::size_t dim = 3;

    std::size_t i = 0, n = lines.size();
    auto strip = [](const std::string& s) {
        std::size_t a = s.find_first_not_of(" \t\r");
        std::size_t b = s.find_last_not_of(" \t\r");
        return a == std::string::npos ? std::string() : s.substr(a, b - a + 1);
    };

    while (i < n) {
        if (strip(lines[i]) != "-1") {
            ++i;
            continue;
        }
        ++i;
        if (i >= n)
            break;
        int ds = std::atoi(strip(lines[i]).c_str());
        ++i;
        std::size_t start = i;
        while (i < n && strip(lines[i]) != "-1")
            ++i;
        std::size_t end = i;  // exclusive
        ++i;                  // skip closing -1

        if (ds == 2411 || ds == 781) {
            std::size_t k = start;
            while (k + 1 < end) {
                auto r1 = tokens(lines[k]);
                if (r1.empty()) {
                    ++k;
                    continue;
                }
                std::int64_t label = std::strtoll(r1[0].c_str(), nullptr, 10);
                auto co = tokens(lines[k + 1]);
                std::vector<double> p;
                for (const auto& c : co)
                    p.push_back(parse_coord(c));
                if (points.empty() && !p.empty())
                    dim = p.size();
                label_to_index[label] = static_cast<std::int64_t>(points.size());
                points.push_back(std::move(p));
                k += 2;
            }
        } else if (ds == 2412) {
            std::size_t k = start;
            while (k < end) {
                auto r1 = tokens(lines[k]);
                if (r1.size() < 6)
                    break;
                std::int64_t elabel = std::strtoll(r1[0].c_str(), nullptr, 10);
                int fedesc = std::atoi(r1[1].c_str());
                std::int64_t pid = std::strtoll(r1[2].c_str(), nullptr, 10);
                int num_nodes = std::atoi(r1[5].c_str());
                ++k;
                if (is_beam(fedesc))
                    ++k;  // skip orientation
                std::vector<std::int64_t> nl;
                while (static_cast<int>(nl.size()) < num_nodes && k < end) {
                    for (const auto& v : tokens(lines[k]))
                        nl.push_back(std::strtoll(v.c_str(), nullptr, 10));
                    ++k;
                }
                nl.resize(num_nodes);
                std::string mtype = unv_type(fedesc);
                if (mtype.empty()) {
                    log::warn("UNV: FE descriptor {} not supported; skipping element.", fedesc);
                    continue;
                }
                std::vector<std::int64_t> unv_conn(num_nodes);
                for (int j = 0; j < num_nodes; ++j)
                    unv_conn[j] = label_to_index.at(nl[j]);
                const std::vector<int>* nd = nd_perm(mtype);
                std::vector<std::int64_t> conn(num_nodes);
                if (nd)
                    for (int j = 0; j < num_nodes; ++j)
                        conn[(*nd)[j]] = unv_conn[j];
                else
                    conn = unv_conn;

                auto git = group_index.find(mtype);
                if (git == group_index.end()) {
                    group_index[mtype] = groups.size();
                    groups.push_back({mtype, {}, {}});
                    git = group_index.find(mtype);
                }
                std::size_t blk = git->second;
                std::size_t local = groups[blk].mRows.size();
                groups[blk].mRows.push_back(std::move(conn));
                groups[blk].mPid.push_back(pid);
                elem_label_to_ref[elabel] = {blk, local};
            }
        } else if (ds == 2467 || ds == 2477 || ds == 2452 || ds == 2435 || ds == 2432 ||
                   ds == 2430) {
            // permanent groups: record1 (>=8 ints; field 7 = entity count),
            // record2 (name), then 4*n ints laid out (entity_type, tag, 0, 0).
            std::size_t k = start;
            while (k < end) {
                auto r1 = tokens(lines[k]);
                if (r1.size() < 8)
                    break;
                int n_ent = std::atoi(r1[7].c_str());
                ++k;
                std::string name = k < end ? lines[k] : std::string();
                std::size_t a = name.find_first_not_of(" \t\r");
                std::size_t b = name.find_last_not_of(" \t\r");
                name = a == std::string::npos ? std::string() : name.substr(a, b - a + 1);
                ++k;
                std::vector<std::int64_t> vals;
                while (static_cast<int>(vals.size()) < 4 * n_ent && k < end) {
                    for (const auto& v : tokens(lines[k]))
                        vals.push_back(std::strtoll(v.c_str(), nullptr, 10));
                    ++k;
                }
                std::vector<std::pair<int, std::int64_t>> ents;
                for (int e = 0; e < n_ent && 4 * e + 1 < static_cast<int>(vals.size()); ++e)
                    ents.emplace_back(static_cast<int>(vals[4 * e]), vals[4 * e + 1]);
                raw_groups.emplace_back(std::move(name), std::move(ents));
            }
        } else if (ds == 2414 || ds == 55 || ds == 56 || ds == 57) {
            Field fld;
            if (parse_field(ds, lines, start, end, fld))
                fields.push_back(std::move(fld));
        }
        // other datasets ignored
    }

    Mesh mesh;
    const std::size_t np = points.size();
    NDArray pts(DType::Float64, {np, dim});
    for (std::size_t r = 0; r < np; ++r)
        for (std::size_t c = 0; c < dim && c < points[r].size(); ++c)
            pts.As<double>()[r * dim + c] = points[r][c];
    mesh.AssignPoints(std::move(pts));

    std::vector<std::size_t> block_sizes;
    std::vector<NDArray> pids;
    for (auto& g : groups) {
        std::size_t ne = g.mRows.size();
        std::size_t k = ne ? g.mRows[0].size() : 0;
        NDArray data(DType::Int64, {ne, k});
        for (std::size_t r = 0; r < ne; ++r)
            for (std::size_t j = 0; j < k; ++j)
                data.As<std::int64_t>()[r * k + j] = g.mRows[r][j];
        mesh.AddCellBlock(g.mType, std::move(data));
        NDArray pd(DType::Int64, {ne});
        for (std::size_t r = 0; r < ne; ++r)
            pd.As<std::int64_t>()[r] = g.mPid[r];
        pids.push_back(std::move(pd));
        block_sizes.push_back(ne);
    }
    if (!pids.empty())
        mesh.AddCellData("unv:pid", std::move(pids));

    // fields -> point_data (location 1) / cell_data (location 2)
    auto unique_key = [&](const std::string& base) {
        std::string name = base.empty() ? "unv:field" : base;
        std::string key = name;
        int n = 1;
        while (used_keys.count(key))
            key = name + "_" + std::to_string(++n);
        used_keys.insert(key);
        return key;
    };
    for (auto& fld : fields) {
        std::string key = unique_key(fld.mName);
        std::size_t nc = static_cast<std::size_t>(fld.mNcomp);
        if (fld.mLocation == 1) {
            NDArray arr =
                nc == 1 ? NDArray(DType::Float64, {np}) : NDArray(DType::Float64, {np, nc});
            std::fill(arr.As<double>(), arr.As<double>() + np * nc, 0.0);
            for (auto& kv : fld.mValues) {
                auto it = label_to_index.find(kv.first);
                if (it == label_to_index.end())
                    continue;
                std::size_t idx = static_cast<std::size_t>(it->second);
                for (std::size_t c = 0; c < nc && c < kv.second.size(); ++c)
                    arr.As<double>()[idx * nc + c] = kv.second[c];
            }
            mesh.AddPointData(key, std::move(arr));
        } else if (fld.mLocation == 2) {
            std::vector<NDArray> blocks;
            for (std::size_t b = 0; b < block_sizes.size(); ++b) {
                std::size_t ne = block_sizes[b];
                NDArray arr =
                    nc == 1 ? NDArray(DType::Float64, {ne}) : NDArray(DType::Float64, {ne, nc});
                std::fill(arr.As<double>(), arr.As<double>() + ne * nc, 0.0);
                blocks.push_back(std::move(arr));
            }
            for (auto& kv : fld.mValues) {
                auto it = elem_label_to_ref.find(kv.first);
                if (it == elem_label_to_ref.end())
                    continue;
                std::size_t b = it->second.first, local = it->second.second;
                for (std::size_t c = 0; c < nc && c < kv.second.size(); ++c)
                    blocks[b].As<double>()[local * nc + c] = kv.second[c];
            }
            mesh.AddCellData(key, std::move(blocks));
        }
    }

    // resolve permanent groups -> point_sets (node groups) / cell_sets
    // (element groups, split per cell block).
    for (auto& g : raw_groups) {
        bool has_node = false, has_elem = false;
        for (auto& e : g.second) {
            if (e.first == 8)
                has_node = true;
            else if (e.first == 7)
                has_elem = true;
        }
        if (has_node) {
            std::vector<std::int64_t> idx;
            for (auto& e : g.second) {
                if (e.first != 8)
                    continue;
                auto it = label_to_index.find(e.second);
                if (it != label_to_index.end())
                    idx.push_back(it->second);
            }
            rInfo.mPointSets[g.first] = std::move(idx);
        }
        if (has_elem) {
            std::vector<std::vector<std::int64_t>> blocks(block_sizes.size());
            for (auto& e : g.second) {
                if (e.first != 7)
                    continue;
                auto it = elem_label_to_ref.find(e.second);
                if (it == elem_label_to_ref.end())
                    continue;
                blocks[it->second.first].push_back(static_cast<std::int64_t>(it->second.second));
            }
            rInfo.mCellSets[g.first] = std::move(blocks);
        }
    }
    return mesh;
}

Mesh read_unv(const std::string& rPath) {
    UnvInfo info;
    return read_unv(rPath, info);
}

void write_unv(const std::string& rPath, const Mesh& rMesh, bool code_aster, int node_dataset) {
    UnvInfo info;
    write_unv(rPath, rMesh, info, code_aster, node_dataset);
}

void write_unv(const std::string& rPath, const Mesh& rMesh, const UnvInfo& rInfo, bool code_aster,
               int node_dataset) {
    std::ofstream f(rPath, std::ios::binary);
    if (!f)
        throw WriteError("Could not open file for writing: " + rPath);

    const std::size_t np = rMesh.NumPoints();
    const std::size_t pdim = rMesh.PointDim();
    const NDArray& points = rMesh.Points();

    if (node_dataset != 2411 && node_dataset != 781)
        node_dataset = 2411;

    // 2411 / 781 nodes
    char buf[128];
    std::snprintf(buf, sizeof(buf), "    -1\n%6d\n", node_dataset);
    f << buf;
    for (std::size_t k = 0; k < np; ++k) {
        std::snprintf(buf, sizeof(buf), "%10zu%10d%10d%10d\n", k + 1, 1, 1, 11);
        f << buf;
        for (int c = 0; c < 3; ++c) {
            double v = c < static_cast<int>(pdim) ? detail::read_double(points, k * pdim + c) : 0.0;
            std::snprintf(buf, sizeof(buf), "%25.16E", v);
            f << buf;
        }
        f << "\n";
    }
    f << "    -1\n";

    // 2412 elements
    f << "    -1\n  2412\n";
    const bool has_pid = rMesh.HasCellData("unv:pid");
    std::int64_t label = 0;
    const std::size_t nblocks = rMesh.NumCellBlocks();
    // per-block 1-based element labels (empty for skipped blocks) so element
    // field data can resolve (block, local) -> label on write.
    std::vector<std::vector<std::int64_t>> block_labels(nblocks);
    for (std::size_t bi = 0; bi < nblocks; ++bi) {
        const auto cb = rMesh.Cells(bi);
        int desc;
        bool beam;
        if (!meshio_descriptor(cb.Type(), desc, beam)) {
            log::warn("UNV does not support '{}' cells. Skipping.", cb.Type());
            continue;
        }
        const NDArray& conn = cb.Conn();
        const std::vector<int>* nd = nd_perm(cb.Type());
        std::size_t ncols = detail::cols(conn);
        std::size_t nrows = cb.NumCells();
        block_labels[bi].reserve(nrows);
        const NDArray* pid = (has_pid && bi < rMesh.CellDataNumBlocks("unv:pid"))
                                 ? &rMesh.CellData("unv:pid", bi)
                                 : nullptr;
        for (std::size_t r = 0; r < nrows; ++r) {
            ++label;
            block_labels[bi].push_back(label);
            std::int64_t pval = pid ? detail::read_int(*pid, r) : 1;
            std::snprintf(buf, sizeof(buf), "%10lld%10d%10lld%10lld%10d%10zu\n",
                          static_cast<long long>(label), desc, static_cast<long long>(pval),
                          static_cast<long long>(pval), 11, ncols);
            f << buf;
            if (beam)
                f << "         0         0         0\n";
            // reorder meshio -> UNV, 1-based, 8 per line
            std::vector<std::int64_t> unv(ncols);
            for (std::size_t j = 0; j < ncols; ++j) {
                std::int64_t node = detail::read_int(conn, r * ncols + j) + 1;
                if (nd)
                    unv[j] = 0;  // filled below
                else
                    unv[j] = node;
            }
            if (nd)
                for (std::size_t j = 0; j < ncols; ++j)
                    unv[j] = detail::read_int(conn, r * ncols + (*nd)[j]) + 1;
            for (std::size_t j = 0; j < ncols; ++j) {
                std::snprintf(buf, sizeof(buf), "%10lld", static_cast<long long>(unv[j]));
                f << buf;
                if ((j + 1) % 8 == 0 || j + 1 == ncols)
                    f << "\n";
            }
        }
    }
    f << "    -1\n";

    // permanent groups (dataset 2467) from rInfo's point/cell sets
    if (!rInfo.mPointSets.empty() || !rInfo.mCellSets.empty()) {
        f << "    -1\n  2467\n";
        int gid = 0;
        auto write_group = [&](const std::string& name, int entity_type,
                               const std::vector<std::int64_t>& tags) {
            ++gid;
            std::snprintf(buf, sizeof(buf), "%10d%10d%10d%10d%10d%10d%10d%10zu\n", gid, 0, 0, 0, 0,
                          0, 0, tags.size());
            f << buf << name << "\n";
            std::size_t col = 0;
            for (std::int64_t t : tags) {
                std::snprintf(buf, sizeof(buf), "%10d%10lld%10d%10d", entity_type,
                              static_cast<long long>(t), 0, 0);
                f << buf;
                if (++col == 2) {
                    f << "\n";
                    col = 0;
                }
            }
            if (col != 0)
                f << "\n";
        };
        for (const auto& kv : rInfo.mPointSets) {
            std::vector<std::int64_t> tags;
            tags.reserve(kv.second.size());
            for (std::int64_t i : kv.second)
                tags.push_back(i + 1);  // 1-based node labels
            write_group(kv.first, 8, tags);
        }
        for (const auto& kv : rInfo.mCellSets) {
            std::vector<std::int64_t> tags;
            for (std::size_t bi = 0; bi < kv.second.size() && bi < block_labels.size(); ++bi)
                for (std::int64_t local : kv.second[bi])
                    if (local >= 0 && local < static_cast<std::int64_t>(block_labels[bi].size()))
                        tags.push_back(block_labels[bi][local]);
            write_group(kv.first, 7, tags);
        }
        f << "    -1\n";
    }

    // field datasets from point_data (nodes) / cell_data (elements)
    auto write_values = [&](const std::vector<std::int64_t>& labels,
                            const std::vector<double>& flat, std::size_t nc) {
        for (std::size_t r = 0; r < labels.size(); ++r) {
            std::snprintf(buf, sizeof(buf), "%10lld\n", static_cast<long long>(labels[r]));
            f << buf;
            for (std::size_t c = 0; c < nc; ++c) {
                std::snprintf(buf, sizeof(buf), "%13.5E", flat[r * nc + c]);
                f << buf;
            }
            f << "\n";
        }
    };
    auto write_field = [&](int field_id, const std::string& name, int location, std::size_t nc,
                           const std::vector<std::int64_t>& labels,
                           const std::vector<double>& flat) {
        int ch = field_char(static_cast<int>(nc));
        if (code_aster) {
            int ds = (location == 1) ? 55 : 57;
            std::snprintf(buf, sizeof(buf), "    -1\n%6d\n", ds);
            f << buf;
            for (int i = 0; i < 5; ++i)
                f << name << "\n";
            std::snprintf(buf, sizeof(buf), "%10d%10d%10d%10d%10d%10zu\n", 1, 0, ch, 0, 4, nc);
            f << buf;
            for (int i = 0; i < 8; ++i)
                f << "         0";
            f << "\n         0         0\n";
            for (int i = 0; i < 6; ++i)
                f << "  0.00000E+00";
            f << "\n";
            for (int i = 0; i < 6; ++i)
                f << "  0.00000E+00";
            f << "\n";
        } else {
            f << "    -1\n  2414\n";
            std::snprintf(buf, sizeof(buf), "%10d\n", field_id);
            f << buf;
            f << name << "\n";
            std::snprintf(buf, sizeof(buf), "%10d\n", location);
            f << buf;
            for (int i = 0; i < 5; ++i)
                f << "meshioplusplus\n";
            std::snprintf(buf, sizeof(buf), "%10d%10d%10d%10d%10d%10zu\n", 1, 0, ch, 0, 4, nc);
            f << buf;
            for (int i = 0; i < 8; ++i)
                f << "         0";
            f << "\n         0         0\n";
            for (int i = 0; i < 6; ++i)
                f << "  0.00000E+00";
            f << "\n";
            for (int i = 0; i < 6; ++i)
                f << "  0.00000E+00";
            f << "\n";
        }
        write_values(labels, flat, nc);
        f << "    -1\n";
    };

    int field_id = 0;
    std::vector<std::int64_t> node_labels(np);
    for (std::size_t k = 0; k < np; ++k)
        node_labels[k] = static_cast<std::int64_t>(k + 1);
    for (const auto& name : rMesh.PointDataNames()) {
        const NDArray& arr = rMesh.PointData(name);
        std::size_t nc = np ? arr.Size() / np : 0;
        if (nc == 0)
            continue;
        std::vector<double> flat(np * nc);
        for (std::size_t i = 0; i < np * nc; ++i)
            flat[i] = detail::read_double(arr, i);
        write_field(++field_id, name, 1, nc, node_labels, flat);
    }
    for (const auto& name : rMesh.CellDataNames()) {
        if (name == "unv:pid")
            continue;
        // gather labels + values across blocks
        std::vector<std::int64_t> labels;
        std::vector<double> flat;
        std::size_t nc = 0;
        for (std::size_t bi = 0; bi < nblocks; ++bi) {
            if (block_labels[bi].empty())
                continue;
            const NDArray& blk = rMesh.CellData(name, bi);
            std::size_t ne = block_labels[bi].size();
            std::size_t bnc = ne ? blk.Size() / ne : 0;
            if (bnc == 0)
                continue;
            if (nc == 0)
                nc = bnc;
            for (std::size_t r = 0; r < ne; ++r) {
                labels.push_back(block_labels[bi][r]);
                for (std::size_t c = 0; c < nc; ++c)
                    flat.push_back(detail::read_double(blk, r * nc + c));
            }
        }
        if (nc == 0 || labels.empty())
            continue;
        write_field(++field_id, name, 2, nc, labels, flat);
    }
}

}  // namespace meshioplusplus
