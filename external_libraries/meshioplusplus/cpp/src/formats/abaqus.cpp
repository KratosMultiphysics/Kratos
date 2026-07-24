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
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/abaqus.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/types.hpp"

namespace meshioplusplus {

namespace {

// (abaqus type, meshio type) in source order; the meshio->abaqus inverse keeps
// the last entry per meshio type (matching the Python dict comprehension).
const std::vector<std::pair<std::string, std::string>>& type_table() {
    static const std::vector<std::pair<std::string, std::string>> t = {
        {"T2D2", "line"},
        {"T2D2H", "line"},
        {"T2D3", "line3"},
        {"T2D3H", "line3"},
        {"T3D2", "line"},
        {"T3D2H", "line"},
        {"T3D3", "line3"},
        {"T3D3H", "line3"},
        {"B21", "line"},
        {"B21H", "line"},
        {"B22", "line3"},
        {"B22H", "line3"},
        {"B31", "line"},
        {"B31H", "line"},
        {"B32", "line3"},
        {"B32H", "line3"},
        {"B33", "line3"},
        {"B33H", "line3"},
        {"CPS4", "quad"},
        {"CPS4R", "quad"},
        {"S4", "quad"},
        {"S4R", "quad"},
        {"S4RS", "quad"},
        {"S4RSW", "quad"},
        {"S4R5", "quad"},
        {"S8R", "quad8"},
        {"S8R5", "quad8"},
        {"S9R5", "quad9"},
        {"CPS3", "triangle"},
        {"STRI3", "triangle"},
        {"S3", "triangle"},
        {"S3R", "triangle"},
        {"S3RS", "triangle"},
        {"R3D3", "triangle"},
        {"STRI65", "triangle6"},
        {"C3D8", "hexahedron"},
        {"C3D8H", "hexahedron"},
        {"C3D8I", "hexahedron"},
        {"C3D8IH", "hexahedron"},
        {"C3D8R", "hexahedron"},
        {"C3D8RH", "hexahedron"},
        {"C3D20", "hexahedron20"},
        {"C3D20H", "hexahedron20"},
        {"C3D20R", "hexahedron20"},
        {"C3D20RH", "hexahedron20"},
        {"C3D4", "tetra"},
        {"C3D4H", "tetra4"},
        {"C3D10", "tetra10"},
        {"C3D10H", "tetra10"},
        {"C3D10I", "tetra10"},
        {"C3D10M", "tetra10"},
        {"C3D10MH", "tetra10"},
        {"C3D6", "wedge"},
        {"C3D15", "wedge15"},
        {"CAX4P", "quad"},
        {"CPE6", "triangle6"},
    };
    return t;
}

const std::unordered_map<std::string, std::string>& abaqus_to_meshio() {
    static const std::unordered_map<std::string, std::string> m = [] {
        std::unordered_map<std::string, std::string> r;
        for (const auto& kv : type_table())
            r[kv.first] = kv.second;
        return r;
    }();
    return m;
}

const std::unordered_map<std::string, std::string>& meshio_to_abaqus() {
    static const std::unordered_map<std::string, std::string> m = [] {
        std::unordered_map<std::string, std::string> r;
        for (const auto& kv : type_table())
            r[kv.second] = kv.first;  // last wins
        return r;
    }();
    return m;
}

std::string abaqus_upper(std::string s) {
    for (auto& c : s)
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return s;
}
std::string abaqus_trim(const std::string& rS) {
    std::size_t b = 0, e = rS.size();
    while (b < e && std::isspace(static_cast<unsigned char>(rS[b])))
        ++b;
    while (e > b && std::isspace(static_cast<unsigned char>(rS[e - 1])))
        --e;
    return rS.substr(b, e - b);
}
std::vector<std::string> split(const std::string& rS, char sep) {
    std::vector<std::string> out;
    std::string cur;
    std::istringstream iss(rS);
    while (std::getline(iss, cur, sep))
        out.push_back(abaqus_trim(cur));
    return out;
}

}  // namespace

Mesh read_abaqus(const std::string& rPath) {
    std::ifstream in(rPath);
    if (!in)
        throw ReadError("Could not open file: " + rPath);
    std::vector<std::string> lines;
    std::string l;
    while (std::getline(in, l)) {
        if (!l.empty() && l.back() == '\r')
            l.pop_back();
        lines.push_back(l);
    }

    Mesh mesh;
    std::unordered_map<std::int64_t, std::int64_t> point_ids;  // file id -> index
    std::vector<std::vector<double>> pts;
    std::size_t dim = 3;
    const auto& a2m = abaqus_to_meshio();

    std::size_t i = 0;
    while (i < lines.size()) {
        const std::string& line = lines[i];
        if (line.rfind("**", 0) == 0) {  // comment
            ++i;
            continue;
        }
        std::string kw = abaqus_upper(abaqus_trim(split(line, ',')[0]));
        if (!kw.empty() && kw[0] == '*')
            kw = kw.substr(1);

        if (kw == "NODE") {
            ++i;
            while (i < lines.size() && (lines[i].empty() || lines[i][0] != '*')) {
                std::string row = abaqus_trim(lines[i]);
                ++i;
                if (row.empty())
                    continue;
                std::vector<std::string> tok = split(row, ',');
                std::int64_t id = std::strtoll(tok[0].c_str(), nullptr, 10);
                point_ids[id] = static_cast<std::int64_t>(pts.size());
                std::vector<double> c;
                for (std::size_t k = 1; k < tok.size(); ++k)
                    if (!tok[k].empty())
                        c.push_back(std::strtod(tok[k].c_str(), nullptr));
                pts.push_back(std::move(c));
            }
        } else if (kw == "ELEMENT") {
            // TYPE= parameter
            std::string etype;
            for (const auto& p : split(line, ',')) {
                std::vector<std::string> kv = split(p, '=');
                if (kv.size() == 2 && abaqus_upper(kv[0]) == "TYPE")
                    etype = kv[1];
            }
            if (etype.empty())
                throw ReadError("Abaqus ELEMENT without TYPE");
            auto it = a2m.find(abaqus_upper(etype));
            // abaqus types are case-sensitive in file; try as-is too
            if (it == a2m.end())
                it = a2m.find(etype);
            if (it == a2m.end())
                throw ReadError("Abaqus element type not supported: " + etype);
            std::string mtype = it->second;
            int n = num_nodes_per_cell().count(mtype) ? num_nodes_per_cell().at(mtype) : 0;
            if (n == 0)
                throw ReadError("Abaqus: unknown node count for " + mtype);
            ++i;
            std::vector<std::int64_t> vals;
            while (i < lines.size() && (lines[i].empty() || lines[i][0] != '*')) {
                std::string row = abaqus_trim(lines[i]);
                ++i;
                if (row.empty())
                    continue;
                for (const auto& t : split(row, ','))
                    if (!t.empty())
                        vals.push_back(std::strtoll(t.c_str(), nullptr, 10));
            }
            std::size_t stride = static_cast<std::size_t>(n) + 1;
            if (vals.size() % stride != 0)
                throw ReadError("Abaqus: bad element data");
            std::size_t ncells = vals.size() / stride;
            NDArray data(DType::Int64, {ncells, static_cast<std::size_t>(n)});
            std::int64_t* dp = data.As<std::int64_t>();
            for (std::size_t r = 0; r < ncells; ++r)
                for (int j = 0; j < n; ++j) {
                    std::int64_t node = vals[r * stride + 1 + j];
                    auto pit = point_ids.find(node);
                    if (pit == point_ids.end())
                        throw ReadError("Abaqus: unknown node id");
                    dp[r * n + j] = pit->second;
                }
            mesh.AddCellBlock(mtype, std::move(data));
        } else if (kw == "NSET" || kw == "ELSET" || kw == "INCLUDE") {
            throw ReadError("Abaqus " + kw + " not supported by the C++ reader");
        } else {
            ++i;  // skip unknown keyword line; its data lines are skipped below
            while (i < lines.size() && (lines[i].empty() || lines[i][0] != '*'))
                ++i;
        }
    }

    if (!pts.empty()) {
        dim = pts[0].size();
        if (dim == 0)
            dim = 3;
    }
    NDArray points(DType::Float64, {pts.size(), dim});
    double* pp = points.As<double>();
    for (std::size_t r = 0; r < pts.size(); ++r)
        for (std::size_t c = 0; c < dim; ++c)
            pp[r * dim + c] = (c < pts[r].size()) ? pts[r][c] : 0.0;
    mesh.AssignPoints(std::move(points));

    return mesh;
}

void write_abaqus(const std::string& rPath, const Mesh& rMesh) {
    std::ofstream os(rPath);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);

    const std::size_t n = rMesh.NumPoints();
    const NDArray& points = rMesh.Points();
    const std::size_t dim = points.Shape().size() >= 2 ? points.Shape()[1] : 0;

    os << "*HEADING\n";
    os << "Abaqus DataFile Version 6.14\n";
    os << "written by meshio++ (C++ core)\n";
    os << "*NODE\n";
    {
        // Format node rows in parallel (snprintf per row, bytes unchanged),
        // then stream sequentially.
        std::vector<std::string> rows(n);
        parallel_for(n, [&](std::size_t i) {
            char buf[48];
            std::string& row = rows[i];
            row = std::to_string(i + 1);
            for (std::size_t c = 0; c < dim; ++c) {
                std::snprintf(buf, sizeof(buf), ", %.16e",
                              detail::read_double(points, i * dim + c));
                row += buf;
            }
            row += '\n';
        });
        for (const auto& row : rows)
            os << row;
    }

    const auto& m2a = meshio_to_abaqus();
    std::size_t eid = 0;
    for (const auto cb : rMesh.CellRange()) {
        auto it = m2a.find(cb.Type());
        if (it == m2a.end())
            throw WriteError("Abaqus writer: unsupported cell type " + cb.Type());
        const NDArray& conn = cb.Conn();
        std::size_t k = conn.Shape().size() >= 2 ? conn.Shape()[1] : 1;
        os << "*ELEMENT, TYPE=" << it->second << "\n";
        for (std::size_t r = 0; r < cb.NumCells(); ++r) {
            os << (++eid);
            for (std::size_t j = 0; j < k; ++j)
                os << "," << (detail::read_int(conn, r * k + j) + 1);
            os << "\n";
        }
    }
}

}  // namespace meshioplusplus
