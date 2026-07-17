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
#include <vector>

// Project includes
#include "meshioplusplus/formats/permas.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus {

namespace {

const std::unordered_map<std::string, std::string>& permas_to_meshio() {
    static const std::unordered_map<std::string, std::string> m = {
        {"PLOT1", "vertex"},      {"PLOTL2", "line"},         {"FLA2", "line"},
        {"FLA3", "line3"},        {"PLOTL3", "line3"},        {"BECOS", "line"},
        {"BECOC", "line"},        {"BETAC", "line"},          {"BECOP", "line"},
        {"BETOP", "line"},        {"BEAM2", "line"},          {"FSCPIPE2", "line"},
        {"LOADA4", "quad"},       {"PLOTA4", "quad"},         {"QUAD4", "quad"},
        {"QUAD4S", "quad"},       {"QUAMS4", "quad"},         {"SHELL4", "quad"},
        {"PLOTA8", "quad8"},      {"LOADA8", "quad8"},        {"QUAMS8", "quad8"},
        {"PLOTA9", "quad9"},      {"LOADA9", "quad9"},        {"QUAMS9", "quad9"},
        {"PLOTA3", "triangle"},   {"SHELL3", "triangle"},     {"TRIA3", "triangle"},
        {"TRIA3K", "triangle"},   {"TRIA3S", "triangle"},     {"TRIMS3", "triangle"},
        {"LOADA6", "triangle6"},  {"TRIMS6", "triangle6"},    {"HEXE8", "hexahedron"},
        {"HEXFO8", "hexahedron"}, {"HEXE20", "hexahedron20"}, {"HEXE27", "hexahedron27"},
        {"TET4", "tetra"},        {"TET10", "tetra10"},       {"PYRA5", "pyramid"},
        {"PENTA6", "wedge"},      {"PENTA15", "wedge15"}};
    return m;
}

// meshio -> permas (last-wins over insertion order, matching the Python reverse map).
const std::unordered_map<std::string, std::string>& meshio_to_permas() {
    static const std::unordered_map<std::string, std::string> m = {
        {"vertex", "PLOT1"},        {"line", "FSCPIPE2"},       {"line3", "PLOTL3"},
        {"quad", "SHELL4"},         {"quad8", "QUAMS8"},        {"quad9", "QUAMS9"},
        {"triangle", "TRIMS3"},     {"triangle6", "TRIMS6"},    {"hexahedron", "HEXFO8"},
        {"hexahedron20", "HEXE20"}, {"hexahedron27", "HEXE27"}, {"tetra", "TET4"},
        {"tetra10", "TET10"},       {"pyramid", "PYRA5"},       {"wedge", "PENTA6"},
        {"wedge15", "PENTA15"}};
    return m;
}

// write-side meshio -> permas node reorders for second-order elements
const std::vector<int>* write_reorder(const std::string& rType) {
    static const std::vector<int> tria6 = {0, 3, 1, 4, 2, 5};
    static const std::vector<int> tet10 = {0, 4, 1, 5, 2, 6, 7, 8, 9, 3};
    static const std::vector<int> quad9 = {0, 4, 1, 7, 8, 5, 3, 6, 2};
    static const std::vector<int> wedge15 = {0, 6, 1, 7, 2, 8, 9, 10, 11, 3, 12, 4, 13, 5, 14};
    if (rType == "triangle6")
        return &tria6;
    if (rType == "tetra10")
        return &tet10;
    if (rType == "quad9")
        return &quad9;
    if (rType == "wedge15")
        return &wedge15;
    return nullptr;
}

std::vector<std::string> split_ws(const std::string& rS) {
    std::vector<std::string> out;
    std::istringstream iss(rS);
    std::string t;
    while (iss >> t)
        out.push_back(t);
    return out;
}

std::string upper(std::string s) {
    for (char& c : s)
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return s;
}

// "$COOR" -> "COOR", "$ELEMENT TYPE=QUAD4" -> "ELEMENT TYPE=QUAD4" (uppercased).
std::string keyword_of(const std::string& rLine) {
    std::size_t a = 0, b = rLine.size();
    while (a < b && (rLine[a] == '$' || std::isspace(static_cast<unsigned char>(rLine[a]))))
        ++a;
    while (b > a && (rLine[b - 1] == '$' || std::isspace(static_cast<unsigned char>(rLine[b - 1]))))
        --b;
    return upper(rLine.substr(a, b - a));
}

}  // namespace

Mesh read_permas(const std::string& rPath) {
    std::ifstream in(rPath, std::ios::binary);
    if (!in)
        throw ReadError("Could not open file: " + rPath);
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(in, line))
        lines.push_back(line);

    Mesh mesh;
    std::vector<double> points;
    std::size_t ncoord = 3;
    std::unordered_map<std::int64_t, std::int64_t> point_gids;
    std::int64_t pindex = 0;

    std::size_t pos = 0;
    const std::size_t n = lines.size();
    while (pos < n) {
        const std::string& cur = lines[pos];
        if (!cur.empty() && cur[0] == '!') {
            ++pos;
            continue;
        }
        std::string kw = keyword_of(cur);
        ++pos;
        if (kw.rfind("COOR", 0) == 0) {
            while (pos < n) {
                const std::string& l = lines[pos];
                if (!l.empty() && (l[0] == '!' || l[0] == '$'))
                    break;
                std::vector<std::string> e = split_ws(l);
                if (e.empty()) {
                    ++pos;
                    continue;
                }
                std::int64_t gid = std::strtoll(e[0].c_str(), nullptr, 10);
                point_gids[gid] = pindex++;
                if (points.empty())
                    ncoord = e.size() - 1;
                for (std::size_t j = 1; j < e.size(); ++j)
                    points.push_back(std::strtod(e[j].c_str(), nullptr));
                ++pos;
            }
        } else if (kw.rfind("ELEMENT", 0) == 0) {
            // parse TYPE=<etype>
            std::size_t eq = kw.find('=');
            if (eq == std::string::npos)
                throw ReadError("PERMAS: $ELEMENT without TYPE=");
            std::string etype =
                upper(split_ws(kw.substr(eq + 1)).empty() ? std::string()
                                                          : split_ws(kw.substr(eq + 1))[0]);
            auto tit = permas_to_meshio().find(etype);
            if (tit == permas_to_meshio().end())
                throw ReadError("PERMAS: element type not available: " + etype);
            const std::string& cell_type = tit->second;

            std::vector<std::vector<std::int64_t>> rows;
            std::vector<std::int64_t> acc;  // accumulates across "!" continuation lines
            while (pos < n) {
                const std::string& l = lines[pos];
                if (!l.empty() && l[0] == '$')
                    break;
                std::vector<std::string> e = split_ws(l);
                if (e.empty()) {
                    ++pos;
                    continue;
                }
                // A trailing "!" marks a continuation; the standalone "!"
                // separator line between blocks just yields no nodes.
                bool continued = (e.back() == "!");
                std::size_t last = continued ? e.size() - 1 : e.size();
                for (std::size_t j = 1; j < last; ++j)
                    acc.push_back(point_gids.at(std::strtoll(e[j].c_str(), nullptr, 10)));
                if (!continued) {
                    rows.push_back(std::move(acc));
                    acc.clear();
                }
                ++pos;
            }
            std::size_t k = rows.empty() ? 0 : rows.front().size();
            NDArray data(DType::Int64, {rows.size(), k});
            std::int64_t* dp = data.As<std::int64_t>();
            for (std::size_t r = 0; r < rows.size(); ++r)
                for (std::size_t j = 0; j < k; ++j)
                    dp[r * k + j] = rows[r][j];
            mesh.AddCellBlock(cell_type, std::move(data));
        }
        // all other keywords (NSET/ESET/...) are ignored
    }

    std::int64_t npoints = static_cast<std::int64_t>(point_gids.size());
    NDArray pts(DType::Float64, {static_cast<std::size_t>(npoints), ncoord});
    double* pp = pts.As<double>();
    for (std::size_t i = 0; i < points.size(); ++i)
        pp[i] = points[i];
    mesh.AssignPoints(std::move(pts));

    return mesh;
}

void write_permas(const std::string& rPath, const Mesh& rMesh) {
    std::ofstream f(rPath, std::ios::binary);
    if (!f)
        throw WriteError("Could not open file for writing: " + rPath);

    const std::size_t npts = rMesh.NumPoints();
    const std::size_t pdim = rMesh.PointDim();
    const NDArray& points = rMesh.Points();

    f << "!PERMAS DataFile Version 18.0\n";
    f << "!written by meshio++ (C++ core)\n";
    f << "$ENTER COMPONENT NAME=DFLT_COMP\n";
    f << "$STRUCTURE\n";
    f << "$COOR\n";
    char buf[32];
    for (std::size_t i = 0; i < npts; ++i) {
        f << (i + 1);
        for (int c = 0; c < 3; ++c) {
            double v =
                c < static_cast<int>(pdim) ? detail::read_double(points, i * pdim + c) : 0.0;
            std::snprintf(buf, sizeof(buf), "%.17g", v);
            f << " " << buf;
        }
        f << "\n";
    }

    std::int64_t eid = 0;
    for (const auto cb : rMesh.CellRange()) {
        auto tit = meshio_to_permas().find(cb.Type());
        if (tit == meshio_to_permas().end())
            throw WriteError("PERMAS: unsupported cell type " + cb.Type());
        f << "!\n";
        f << "$ELEMENT TYPE=" << tit->second << "\n";
        const std::vector<int>* reorder = write_reorder(cb.Type());
        const NDArray& conn = cb.Conn();
        const std::size_t ncols = detail::cols(conn);
        const std::size_t nc = cb.NumCells();
        for (std::size_t r = 0; r < nc; ++r) {
            ++eid;
            f << eid;
            if (reorder) {
                for (int local : *reorder)
                    f << " " << (detail::read_int(conn, r * ncols + local) + 1);
            } else {
                for (std::size_t j = 0; j < ncols; ++j)
                    f << " " << (detail::read_int(conn, r * ncols + j) + 1);
            }
            f << "\n";
        }
    }

    f << "$END STRUCTURE\n";
    f << "$EXIT COMPONENT\n";
    f << "$FIN\n";
}

}  // namespace meshioplusplus
