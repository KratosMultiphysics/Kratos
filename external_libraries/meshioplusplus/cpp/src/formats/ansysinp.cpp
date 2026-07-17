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
#include <array>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/ansysinp.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {

namespace {

// ---- Ansys element type id -> family (mirrors _FAMILY in _ansysInp.py) ----
const std::unordered_map<int, std::string>& family_map() {
    static const std::unordered_map<int, std::string> m = [] {
        std::unordered_map<int, std::string> f;
        for (int n : {5, 45, 70, 87, 90, 92, 95, 162, 185, 186, 187, 226, 227, 285})
            f[n] = "solid";
        for (int n : {28, 43, 63, 93, 131, 132, 181, 281})
            f[n] = "shell";
        for (int n : {25, 42, 77, 82, 182, 183, 223})
            f[n] = "plane";
        for (int n : {1, 3, 4, 21, 180, 188, 189, 288, 289})
            f[n] = "line";
        return f;
    }();
    return m;
}

// (family, node count) -> meshio type (mirrors _TO_MESHIO).
std::string to_meshio(const std::string& rFamily, std::size_t nnodes) {
    static const std::map<std::pair<std::string, std::size_t>, std::string> m = {
        {{"solid", 4}, "tetra"},         {{"solid", 10}, "tetra10"},   {{"solid", 8}, "hexahedron"},
        {{"solid", 20}, "hexahedron20"}, {{"solid", 6}, "wedge"},      {{"solid", 15}, "wedge15"},
        {{"solid", 5}, "pyramid"},       {{"solid", 13}, "pyramid13"}, {{"shell", 3}, "triangle"},
        {{"shell", 6}, "triangle6"},     {{"shell", 4}, "quad"},       {{"shell", 8}, "quad8"},
        {{"plane", 3}, "triangle"},      {{"plane", 6}, "triangle6"},  {{"plane", 4}, "quad"},
        {{"plane", 8}, "quad8"},         {{"line", 2}, "line"},        {{"line", 3}, "line3"},
    };
    auto it = m.find({rFamily, nnodes});
    return it == m.end() ? std::string() : it->second;
}

// meshio type -> Ansys element type id on write (mirrors _FROM_MESHIO).
int from_meshio(const std::string& rT) {
    static const std::unordered_map<std::string, int> m = {
        {"tetra", 285},    {"tetra10", 187},   {"hexahedron", 185}, {"hexahedron20", 186},
        {"wedge", 185},    {"wedge15", 186},   {"pyramid", 185},    {"pyramid13", 186},
        {"triangle", 181}, {"triangle6", 281}, {"quad", 181},       {"quad8", 281},
        {"line", 188},     {"line3", 189},
    };
    auto it = m.find(rT);
    return it == m.end() ? -1 : it->second;
}

std::string upper(std::string s) {
    for (char& c : s)
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return s;
}
std::string strip(const std::string& rS) {
    std::size_t a = rS.find_first_not_of(" \t\r\n");
    if (a == std::string::npos)
        return "";
    std::size_t b = rS.find_last_not_of(" \t\r\n");
    return rS.substr(a, b - a + 1);
}

// int field width from a Fortran format spec like "(3i9,6e20.13)" -> 9.
int int_width(const std::string& rFmt) {
    // match (\d+)i(\d+)
    for (std::size_t i = 0; i + 1 < rFmt.size(); ++i) {
        if ((rFmt[i] == 'i' || rFmt[i] == 'I') && i > 0 &&
            std::isdigit(static_cast<unsigned char>(rFmt[i - 1]))) {
            std::size_t j = i + 1;
            std::string num;
            while (j < rFmt.size() && std::isdigit(static_cast<unsigned char>(rFmt[j])))
                num += rFmt[j++];
            if (!num.empty())
                return std::stoi(num);
        }
    }
    return 0;
}
// real field width from "(3i9,6e20.13)" -> 20 (the digits after e/g, before '.').
int real_width(const std::string& rFmt) {
    for (std::size_t i = 0; i + 1 < rFmt.size(); ++i) {
        char c = static_cast<char>(std::tolower(static_cast<unsigned char>(rFmt[i])));
        if ((c == 'e' || c == 'g') && i > 0 &&
            std::isdigit(static_cast<unsigned char>(rFmt[i - 1]))) {
            std::size_t j = i + 1;
            std::string num;
            while (j < rFmt.size() && std::isdigit(static_cast<unsigned char>(rFmt[j])))
                num += rFmt[j++];
            if (j < rFmt.size() && rFmt[j] == '.' && !num.empty())
                return std::stoi(num);
        }
    }
    return 0;
}

// Slice a line into fixed-width integer fields; stop at the first non-numeric
// chunk (mirrors _slice_ints).
std::vector<std::int64_t> slice_ints(const std::string& rLineIn, int width) {
    std::vector<std::int64_t> out;
    std::string line = rLineIn;
    while (!line.empty() && (line.back() == '\n' || line.back() == '\r'))
        line.pop_back();
    for (std::size_t i = 0; i < line.size(); i += static_cast<std::size_t>(width)) {
        std::string chunk = strip(line.substr(i, static_cast<std::size_t>(width)));
        if (chunk.empty())
            continue;
        try {
            std::size_t pos = 0;
            long long v = std::stoll(chunk, &pos);
            if (pos != chunk.size())
                break;  // trailing non-numeric
            out.push_back(v);
        } catch (...) {
            break;
        }
    }
    return out;
}

std::vector<double> slice_reals(const std::string& rS, int width) {
    std::vector<double> out;
    for (std::size_t i = 0; i < rS.size(); i += static_cast<std::size_t>(width)) {
        std::string chunk = strip(rS.substr(i, static_cast<std::size_t>(width)));
        if (chunk.empty())
            continue;
        try {
            out.push_back(std::stod(chunk));
        } catch (...) {
        }
    }
    return out;
}

bool is_data_line(const std::string& rLine) {
    std::string s = strip(rLine);
    if (s.empty())
        return false;
    std::string up = upper(s);
    static const char* kws[] = {"FINISH", "NBLOCK",   "EBLOCK",  "CMBLOCK", "ETBLOCK", "/PREP7",
                                "/SOLU",  "/POST1",   "/EOF",    "KEYOPT",  "MPDATA",  "MPTEMP",
                                "LOCAL",  "SECBLOCK", "RLBLOCK", "DBLOCK",  "FBLOCK",  "SFEBLOCK"};
    for (const char* kw : kws)
        if (up.rfind(kw, 0) == 0)
            return false;
    // ^[A-Z]{1,8}, -> a command line
    std::size_t comma = up.find(',');
    if (comma != std::string::npos && comma >= 1 && comma <= 8) {
        bool all_alpha = true;
        for (std::size_t i = 0; i < comma; ++i)
            if (!std::isalpha(static_cast<unsigned char>(up[i]))) {
                all_alpha = false;
                break;
            }
        if (all_alpha)
            return false;
    }
    if (s[0] == '!' || s[0] == '/')
        return false;
    return true;
}

std::vector<std::string> read_lines_file(const std::string& rPath) {
    std::ifstream f(rPath);
    if (!f)
        throw ReadError("Could not open ansysInp file: " + rPath);
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(f, line)) {
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        lines.push_back(line);
    }
    return lines;
}

}  // namespace

Mesh read_ansysinp(const std::string& rPath, AnsysInfo& rInfo) {
    std::vector<std::string> lines = read_lines_file(rPath);

    std::unordered_map<int, int> etype_lib;  // slot -> ansys element type id
    std::vector<std::int64_t> node_id;
    std::vector<std::array<double, 3>> coords;
    // (etype_local, elem_id, node ids)
    struct Elem {
        int mEtypeLocal;
        std::int64_t mElemId;
        std::vector<std::int64_t> mNodes;
    };
    std::vector<Elem> elements;
    std::vector<std::pair<std::string, std::vector<std::int64_t>>> node_comps;
    std::vector<std::pair<std::string, std::vector<std::int64_t>>> elem_comps;
    bool saw_block = false;

    std::size_t i = 0, n = lines.size();
    while (i < n) {
        std::string line = strip(lines[i]);
        std::string up = upper(line);

        if (up.rfind("ET,", 0) == 0) {
            std::stringstream ss(line);
            std::string tok;
            std::vector<std::string> p;
            while (std::getline(ss, tok, ','))
                p.push_back(tok);
            if (p.size() >= 3) {
                try {
                    etype_lib[std::stoi(strip(p[1]))] = static_cast<int>(std::stod(strip(p[2])));
                } catch (...) {
                }
            }
            ++i;
        } else if (up.rfind("ETBLOCK", 0) == 0) {
            saw_block = true;
            std::string count_field = line.substr(line.find(',') + 1);
            count_field = count_field.substr(0, count_field.find('!'));
            int ntypes = std::stoi(strip(count_field));
            int iw = (i + 1 < n) ? int_width(lines[i + 1]) : 0;
            if (iw == 0)
                iw = 9;
            i += 2;
            int got = 0;
            while (i < n && got < ntypes) {
                if (!is_data_line(lines[i]))
                    break;
                auto v = slice_ints(lines[i], iw);
                if (!v.empty() && v[0] == -1) {
                    ++i;
                    break;
                }
                if (v.size() >= 2) {
                    etype_lib[static_cast<int>(v[0])] = static_cast<int>(v[1]);
                    ++got;
                }
                ++i;
            }
        } else if (up.rfind("NBLOCK", 0) == 0) {
            saw_block = true;
            int iw = (i + 1 < n) ? int_width(lines[i + 1]) : 0;
            if (iw == 0)
                iw = 9;
            int rw = (i + 1 < n) ? real_width(lines[i + 1]) : 0;
            if (rw == 0)
                rw = 20;
            i += 2;
            while (i < n) {
                const std::string& l = lines[i];
                std::string s = upper(strip(l));
                if (s.rfind("N,", 0) == 0 || s.rfind("-1", 0) == 0 || s.empty()) {
                    ++i;
                    break;
                }
                if (!is_data_line(l))
                    break;
                std::int64_t nid;
                try {
                    nid = std::stoll(strip(l.substr(0, static_cast<std::size_t>(iw))));
                } catch (...) {
                    ++i;
                    continue;
                }
                if (nid < 0) {
                    ++i;
                    break;
                }
                std::vector<double> rs =
                    l.size() > static_cast<std::size_t>(3 * iw)
                        ? slice_reals(l.substr(static_cast<std::size_t>(3 * iw)), rw)
                        : std::vector<double>{};
                rs.resize(3, 0.0);
                node_id.push_back(nid);
                coords.push_back({rs[0], rs[1], rs[2]});
                ++i;
            }
        } else if (up.rfind("EBLOCK", 0) == 0) {
            saw_block = true;
            int iw = (i + 1 < n) ? int_width(lines[i + 1]) : 0;
            if (iw == 0)
                iw = 9;
            i += 2;
            while (i < n) {
                const std::string& l = lines[i];
                if (strip(l).rfind("-1", 0) == 0) {
                    ++i;
                    break;
                }
                if (!is_data_line(l))
                    break;
                auto fields = slice_ints(l, iw);
                if (fields.empty()) {
                    ++i;
                    continue;
                }
                int etype_local = static_cast<int>(fields[1]);
                std::size_t nnodes = static_cast<std::size_t>(fields[8]);
                std::int64_t elem_id = fields[10];
                std::vector<std::int64_t> nodes(fields.begin() + 11, fields.end());
                ++i;
                while (nodes.size() < nnodes && i < n) {
                    if (!is_data_line(lines[i]))
                        break;
                    if (strip(lines[i]).rfind("-1", 0) == 0)
                        break;
                    auto more = slice_ints(lines[i], iw);
                    nodes.insert(nodes.end(), more.begin(), more.end());
                    ++i;
                }
                nodes.resize(std::min(nodes.size(), nnodes));
                elements.push_back({etype_local, elem_id, std::move(nodes)});
            }
        } else if (up.rfind("CMBLOCK", 0) == 0) {
            saw_block = true;
            std::stringstream ss(line);
            std::string tok;
            std::vector<std::string> p;
            while (std::getline(ss, tok, ','))
                p.push_back(tok);
            std::string cname = strip(p.at(1));
            std::string entity = upper(strip(p.at(2)));
            std::string cnt = p.at(3);
            cnt = cnt.substr(0, cnt.find('!'));
            std::size_t numitems = static_cast<std::size_t>(std::stoll(strip(cnt)));
            int iw = (i + 1 < n) ? int_width(lines[i + 1]) : 0;
            if (iw == 0)
                iw = 10;
            i += 2;
            std::vector<std::int64_t> items;
            while (i < n && items.size() < numitems) {
                if (!is_data_line(lines[i]))
                    break;
                auto more = slice_ints(lines[i], iw);
                items.insert(items.end(), more.begin(), more.end());
                ++i;
            }
            if (items.size() > numitems)
                items.resize(numitems);
            std::vector<std::int64_t> expanded;
            bool have_prev = false;
            std::int64_t prev = 0;
            for (std::int64_t it : items) {
                if (it < 0) {
                    if (!have_prev)
                        throw ReadError("Invalid CMBLOCK '" + cname +
                                        "': range marker (negative value) before any "
                                        "base value.");
                    for (std::int64_t v = prev + 1; v <= -it; ++v)
                        expanded.push_back(v);
                    prev = -it;
                } else {
                    expanded.push_back(it);
                    prev = it;
                    have_prev = true;
                }
            }
            if (entity.rfind("NODE", 0) == 0)
                node_comps.emplace_back(cname, std::move(expanded));
            else
                elem_comps.emplace_back(cname, std::move(expanded));
        } else {
            ++i;
        }
    }

    if (!saw_block)
        throw ReadError("No MAPDL block (NBLOCK/EBLOCK/CMBLOCK) found.");

    // ---- build mesh ----
    Mesh mesh;
    std::size_t npts = coords.size();
    NDArray pts(DType::Float64, {npts, 3});
    for (std::size_t k = 0; k < npts; ++k)
        for (std::size_t j = 0; j < 3; ++j)
            pts.As<double>()[k * 3 + j] = coords[k][j];
    mesh.AssignPoints(std::move(pts));

    std::unordered_map<std::int64_t, std::int64_t> nid_to_index;
    for (std::size_t k = 0; k < node_id.size(); ++k)
        nid_to_index[node_id[k]] = k;

    // blocks, in first-seen order
    std::vector<std::string> order;
    std::map<std::string, std::vector<std::vector<std::int64_t>>> blocks;
    // element id -> (block index in `order`, local index)
    std::unordered_map<std::int64_t, std::pair<std::size_t, std::size_t>> eid_to_loc;
    for (const Elem& e : elements) {
        auto fam_it =
            family_map().find(etype_lib.count(e.mEtypeLocal) ? etype_lib.at(e.mEtypeLocal) : -1);
        std::string family = fam_it == family_map().end() ? "solid" : fam_it->second;
        std::string mtype = to_meshio(family, e.mNodes.size());
        if (mtype.empty())
            throw ReadError("Unsupported type: etype " + std::to_string(e.mEtypeLocal) + " with " +
                            std::to_string(e.mNodes.size()) + " nodes.");
        if (!blocks.count(mtype))
            order.push_back(mtype);
        auto& blk = blocks[mtype];
        std::size_t bidx = 0;
        for (std::size_t o = 0; o < order.size(); ++o)
            if (order[o] == mtype) {
                bidx = o;
                break;
            }
        eid_to_loc[e.mElemId] = {bidx, blk.size()};
        std::vector<std::int64_t> row;
        row.reserve(e.mNodes.size());
        for (std::int64_t x : e.mNodes)
            row.push_back(nid_to_index.at(x));
        blk.push_back(std::move(row));
    }

    for (const std::string& t : order) {
        const auto& blk = blocks[t];
        std::size_t nc = blk.size();
        std::size_t k = nc ? blk[0].size() : 0;
        NDArray data(DType::Int64, {nc, k});
        for (std::size_t r = 0; r < nc; ++r)
            for (std::size_t c = 0; c < k; ++c)
                data.As<std::int64_t>()[r * k + c] = blk[r][c];
        mesh.AddCellBlock(t, std::move(data));
    }

    // point/cell sets (side-channel)
    for (const auto& kv : node_comps) {
        std::vector<std::int64_t> idx;
        for (std::int64_t x : kv.second)
            if (nid_to_index.count(x))
                idx.push_back(nid_to_index.at(x));
        rInfo.mPointSets[kv.first] = std::move(idx);
    }
    for (const auto& kv : elem_comps) {
        std::vector<std::vector<std::int64_t>> per(order.size());
        for (std::int64_t eid : kv.second) {
            auto it = eid_to_loc.find(eid);
            if (it != eid_to_loc.end())
                per[it->second.first].push_back(static_cast<std::int64_t>(it->second.second));
        }
        rInfo.mCellSets[kv.first] = std::move(per);
    }

    return mesh;
}

void write_ansysinp(const std::string& rPath, const Mesh& rMesh, const AnsysInfo& rInfo) {
    std::ofstream f(rPath);
    if (!f)
        throw WriteError("Could not open ansysInp file for writing: " + rPath);

    const NDArray& points = rMesh.Points();
    std::size_t npts = rMesh.NumPoints();
    std::size_t dim = points.Shape().size() > 1 ? points.Shape()[1] : 3;

    // element-type slots, first-seen order
    std::vector<std::pair<std::string, int>> type_slot;
    auto slot_of = [&](const std::string& t) -> int {
        for (auto& kv : type_slot)
            if (kv.first == t)
                return kv.second;
        int s = static_cast<int>(type_slot.size()) + 1;
        type_slot.emplace_back(t, s);
        return s;
    };
    for (const auto b : rMesh.CellRange()) {
        if (from_meshio(b.Type()) < 0)
            throw WriteError("Unhandled meshio type: " + b.Type());
        slot_of(b.Type());
    }

    f << "/PREP7\n";
    for (auto& kv : type_slot)
        f << "ET," << kv.second << "," << from_meshio(kv.first) << "\n";

    char buf[64];
    std::snprintf(buf, sizeof(buf), "NBLOCK,6,SOLID,%zu,%zu\n(3i9,6e20.13)\n", npts, npts);
    f << buf;
    {
        // Format node rows in parallel (snprintf per row, bytes unchanged),
        // then stream sequentially.
        std::vector<std::string> rows(npts);
        parallel_for(npts, [&](std::size_t k) {
            char b1[32], b2[80];
            double x = dim > 0 ? detail::read_double(points, k * dim + 0) : 0.0;
            double y = dim > 1 ? detail::read_double(points, k * dim + 1) : 0.0;
            double z = dim > 2 ? detail::read_double(points, k * dim + 2) : 0.0;
            std::snprintf(b1, sizeof(b1), "%9zu%9d%9d", k + 1, 0, 0);
            std::snprintf(b2, sizeof(b2), "% .13E% .13E% .13E\n", x, y, z);
            rows[k] = std::string(b1) + b2;
        });
        for (const auto& row : rows)
            f << row;
    }
    f << "N,R5.3,LOC,      -1,\n";

    std::size_t ntot = 0;
    for (const auto b : rMesh.CellRange())
        ntot += b.NumCells();
    std::snprintf(buf, sizeof(buf), "EBLOCK,19,SOLID,%zu,%zu\n(19i9)\n", ntot, ntot);
    f << buf;

    std::int64_t eid = 0;
    // Element ids are consecutive; block_eid_base[bi] is the (exclusive) base id
    // of block bi, so element (bi, li) has id block_eid_base[bi] + 1 + li. This
    // replaces a per-cell std::map lookup with a simple prefix sum.
    std::vector<std::int64_t> block_eid_base(rMesh.NumCellBlocks());
    for (std::size_t bi = 0; bi < rMesh.NumCellBlocks(); ++bi) {
        const auto b = rMesh.Cells(bi);
        const NDArray& conn = b.Conn();
        int slot = slot_of(b.Type());
        std::size_t nc = b.NumCells();
        std::size_t k = detail::cols(conn);
        const std::int64_t eid_base = eid;  // element ids are consecutive
        block_eid_base[bi] = eid_base;
        eid += static_cast<std::int64_t>(nc);
        // Format element rows in parallel, then stream sequentially.
        std::vector<std::string> rows(nc);
        parallel_for(nc, [&](std::size_t li) {
            char fld[32];
            std::string& row = rows[li];
            std::vector<std::int64_t> nodes(k);
            for (std::size_t c = 0; c < k; ++c)
                nodes[c] = detail::read_int(conn, li * k + c) + 1;
            std::vector<std::int64_t> first = {1,
                                               slot,
                                               1,
                                               1,
                                               0,
                                               0,
                                               0,
                                               0,
                                               static_cast<std::int64_t>(k),
                                               0,
                                               eid_base + 1 + static_cast<std::int64_t>(li)};
            for (std::size_t c = 0; c < std::min<std::size_t>(8, k); ++c)
                first.push_back(nodes[c]);
            for (std::int64_t v : first) {
                std::snprintf(fld, sizeof(fld), "%9lld", static_cast<long long>(v));
                row += fld;
            }
            row += '\n';
            if (k > 8) {
                for (std::size_t c = 8; c < k; ++c) {
                    std::snprintf(fld, sizeof(fld), "%9lld", static_cast<long long>(nodes[c]));
                    row += fld;
                }
                row += '\n';
            }
        });
        for (const auto& row : rows)
            f << row;
    }
    std::snprintf(buf, sizeof(buf), "%9d\n", -1);
    f << buf;

    auto write_items = [&](const std::vector<std::int64_t>& vals) {
        for (std::size_t i = 0; i < vals.size(); i += 8) {
            for (std::size_t j = i; j < std::min(i + 8, vals.size()); ++j) {
                std::snprintf(buf, sizeof(buf), "%10lld", static_cast<long long>(vals[j]));
                f << buf;
            }
            f << "\n";
        }
    };

    for (const auto& kv : rInfo.mPointSets) {
        std::vector<std::int64_t> vals;
        for (std::int64_t x : kv.second)
            vals.push_back(x + 1);
        std::snprintf(buf, sizeof(buf), "CMBLOCK,%s,NODE,%9zu\n(8i10)\n", kv.first.c_str(),
                      vals.size());
        f << buf;
        write_items(vals);
    }
    for (const auto& kv : rInfo.mCellSets) {
        std::vector<std::int64_t> vals;
        for (std::size_t bi = 0; bi < kv.second.size(); ++bi)
            for (std::int64_t li : kv.second[bi]) {
                if (bi < block_eid_base.size() &&
                    static_cast<std::size_t>(li) < rMesh.Cells(bi).NumCells())
                    vals.push_back(block_eid_base[bi] + 1 + li);
            }
        std::sort(vals.begin(), vals.end());
        std::snprintf(buf, sizeof(buf), "CMBLOCK,%s,ELEM,%9zu\n(8i10)\n", kv.first.c_str(),
                      vals.size());
        f << buf;
        write_items(vals);
    }
    f << "FINISH\n";
}

}  // namespace meshioplusplus
