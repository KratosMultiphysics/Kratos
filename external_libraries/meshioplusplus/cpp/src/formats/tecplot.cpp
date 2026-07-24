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
#include "meshioplusplus/formats/tecplot.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus {

namespace {

std::string tecplot_upper(std::string s) {
    for (auto& c : s)
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return s;
}
std::string tecplot_strip(const std::string& rS) {
    std::size_t b = rS.find_first_not_of(" \t\r\n");
    if (b == std::string::npos)
        return "";
    std::size_t e = rS.find_last_not_of(" \t\r\n");
    return rS.substr(b, e - b + 1);
}
std::vector<std::string> tecplot_tokens(const std::string& rS) {
    std::vector<std::string> out;
    std::istringstream iss(rS);
    std::string t;
    while (iss >> t)
        out.push_back(t);
    return out;
}
bool is_float_token(const std::string& rS) {
    if (rS.empty())
        return false;
    char* endp = nullptr;
    std::strtod(rS.c_str(), &endp);
    return endp == rS.c_str() + rS.size();
}

std::string tecplot_to_meshio(const std::string& rZ) {
    std::string u = tecplot_upper(rZ);
    if (u == "LINESEG" || u == "FELINESEG")
        return "line";
    if (u == "TRIANGLE" || u == "FETRIANGLE")
        return "triangle";
    if (u == "QUADRILATERAL" || u == "FEQUADRILATERAL")
        return "quad";
    if (u == "TETRAHEDRON" || u == "FETETRAHEDRON")
        return "tetra";
    if (u == "BRICK" || u == "FEBRICK")
        return "hexahedron";
    return "";
}
std::string meshio_to_tecplot(const std::string& rM) {
    if (rM == "line")
        return "FELINESEG";
    if (rM == "triangle")
        return "FETRIANGLE";
    if (rM == "quad")
        return "FEQUADRILATERAL";
    if (rM == "tetra")
        return "FETETRAHEDRON";
    if (rM == "pyramid" || rM == "wedge" || rM == "hexahedron")
        return "FEBRICK";
    return "";
}
const std::vector<int>& tecplot_order(const std::string& rM) {
    static const std::map<std::string, std::vector<int>> o = {
        {"line", {0, 1}},
        {"triangle", {0, 1, 2}},
        {"quad", {0, 1, 2, 3}},
        {"tetra", {0, 1, 2, 3}},
        {"pyramid", {0, 1, 2, 3, 4, 4, 4, 4}},
        {"wedge", {0, 1, 4, 3, 2, 2, 5, 5}},
        {"hexahedron", {0, 1, 2, 3, 4, 5, 6, 7}},
    };
    static const std::vector<int> empty;
    auto it = o.find(rM);
    return it == o.end() ? empty : it->second;
}

}  // namespace

Mesh read_tecplot(const std::string& rPath) {
    std::ifstream in(rPath);
    if (!in)
        throw ReadError("Could not open file: " + rPath);
    std::vector<std::string> lines;
    std::string l;
    while (std::getline(in, l)) {
        std::string s = tecplot_strip(l);
        if (s.empty() || s[0] == '#')
            continue;
        lines.push_back(s);
    }

    std::vector<std::string> variables;
    std::map<std::string, std::string> zone;
    std::string varloc;
    std::size_t i = 0, data_start = lines.size();
    for (; i < lines.size(); ++i) {
        std::string u = tecplot_upper(lines[i]);
        if (u.rfind("VARIABLES", 0) == 0) {
            std::string joined = lines[i];
            while (i + 1 < lines.size() && tecplot_strip(lines[i + 1])[0] == '"')
                joined += " " + lines[++i];
            std::string rhs = joined.substr(joined.find('=') + 1);
            // collect quoted names (or bare tokens)
            std::size_t p = 0;
            while (p < rhs.size()) {
                if (rhs[p] == '"') {
                    std::size_t q = rhs.find('"', p + 1);
                    variables.push_back(rhs.substr(p + 1, q - p - 1));
                    p = q + 1;
                } else if (std::isspace((unsigned char)rhs[p]) || rhs[p] == ',') {
                    ++p;
                } else {
                    std::size_t q = p;
                    while (q < rhs.size() && !std::isspace((unsigned char)rhs[q]) && rhs[q] != ',')
                        ++q;
                    variables.push_back(rhs.substr(p, q - p));
                    p = q;
                }
            }
        } else if (u.rfind("ZONE", 0) == 0) {
            std::string joined = lines[i];
            while (i + 1 < lines.size() && !is_float_token(tecplot_tokens(lines[i + 1])[0]))
                joined += " " + lines[++i];
            data_start = i + 1;
            // Extract VARLOCATION(...)
            std::string ju = joined;
            std::size_t vp = tecplot_upper(ju).find("VARLOCATION");
            if (vp != std::string::npos) {
                std::size_t p1 = ju.find('(', vp), p2 = ju.find(')', p1);
                varloc = ju.substr(p1, p2 - p1 + 1);
                varloc.erase(std::remove(varloc.begin(), varloc.end(), ' '), varloc.end());
                ju = ju.substr(0, vp) + ju.substr(p2 + 1);
            }
            // tokenize key/values (drop ZONE, replace ,/= with space)
            std::string body = ju.substr(4);
            for (auto& c : body)
                if (c == ',' || c == '=')
                    c = ' ';
            auto tk = tecplot_tokens(body);
            for (std::size_t k = 0; k + 1 < tk.size(); ++k) {
                std::string key = tecplot_upper(tk[k]);
                if (key == "NODES" || key == "N" || key == "ELEMENTS" || key == "E" ||
                    key == "DATAPACKING" || key == "ZONETYPE" || key == "F" || key == "ET" ||
                    key == "NV")
                    zone[key] = tk[k + 1];
            }
            break;
        }
    }
    if (variables.empty())
        throw ReadError("Tecplot: no VARIABLES");

    auto getz = [&](const char* a, const char* b) -> std::string {
        if (zone.count(a))
            return zone[a];
        if (zone.count(b))
            return zone[b];
        return "";
    };
    std::size_t num_nodes = std::stoull(getz("NODES", "N"));
    std::size_t num_cells = std::stoull(getz("ELEMENTS", "E"));
    std::string fmt, ztype;
    if (zone.count("F")) {
        fmt = tecplot_upper(zone["F"]);
        ztype = zone.count("ET") ? zone["ET"] : "";
    } else {
        fmt = "FE" + tecplot_upper(getz("DATAPACKING", ""));
        ztype = getz("ZONETYPE", "");
    }
    bool feblock = (fmt == "FEBLOCK");

    std::vector<int> cell_centered(variables.size(), 0);
    if (feblock) {
        if (zone.count("NV")) {
            int nv = std::stoi(zone["NV"]);
            for (std::size_t k = nv; k < variables.size(); ++k)
                cell_centered[k] = 1;
        } else if (!varloc.empty()) {
            std::string vc = varloc.substr(1, varloc.size() - 2);  // strip ()
            for (const auto& entry : [&] {
                     std::vector<std::string> es;
                     std::string cur;
                     for (char c : vc) {
                         if (c == ',') {
                             es.push_back(cur);
                             cur.clear();
                         } else
                             cur += c;
                     }
                     if (!cur.empty())
                         es.push_back(cur);
                     return es;
                 }()) {
                std::size_t eq = entry.find('=');
                if (eq == std::string::npos)
                    continue;
                std::string rng = entry.substr(0, eq), loc = tecplot_upper(entry.substr(eq + 1));
                if (loc != "CELLCENTERED")
                    continue;
                rng = rng.substr(1, rng.size() - 2);  // strip []
                std::size_t dash = rng.find('-');
                if (dash == std::string::npos) {
                    cell_centered[std::stoi(rng) - 1] = 1;
                } else {
                    int a = std::stoi(rng.substr(0, dash)), b = std::stoi(rng.substr(dash + 1));
                    for (int k = a; k <= b; ++k)
                        cell_centered[k - 1] = 1;
                }
            }
        }
    }

    // Read data values.
    std::vector<std::size_t> ndata(variables.size());
    std::size_t total = 0;
    for (std::size_t k = 0; k < variables.size(); ++k) {
        ndata[k] = cell_centered[k] ? num_cells : num_nodes;
        total += ndata[k];
    }
    std::size_t want = feblock ? total : num_nodes * variables.size();

    std::vector<double> flat;
    flat.reserve(want);
    std::size_t li = data_start;
    while (flat.size() < want && li < lines.size()) {
        for (const auto& t : tecplot_tokens(lines[li]))
            flat.push_back(std::strtod(t.c_str(), nullptr));
        ++li;
    }

    // Per-variable columns.
    std::vector<std::vector<double>> cols(variables.size());
    if (feblock) {
        std::size_t off = 0;
        for (std::size_t k = 0; k < variables.size(); ++k) {
            cols[k].assign(flat.begin() + off, flat.begin() + off + ndata[k]);
            off += ndata[k];
        }
    } else {
        std::size_t nv = variables.size();
        for (std::size_t k = 0; k < nv; ++k)
            cols[k].resize(num_nodes);
        for (std::size_t r = 0; r < num_nodes; ++r)
            for (std::size_t k = 0; k < nv; ++k)
                cols[k][r] = flat[r * nv + k];
    }

    // Cells.
    std::string mtype = tecplot_to_meshio(ztype);
    if (mtype.empty())
        throw ReadError("Tecplot: unsupported zone type " + ztype);
    std::size_t nn;
    if (mtype == "line")
        nn = 2;
    else if (mtype == "triangle")
        nn = 3;
    else if (mtype == "quad" || mtype == "tetra")
        nn = 4;
    else
        nn = 8;
    NDArray celldata(DType::Int64, {num_cells, nn});
    std::int64_t* cp = celldata.As<std::int64_t>();
    for (std::size_t c = 0; c < num_cells; ++c) {
        auto t = tecplot_tokens(lines.at(li++));
        for (std::size_t j = 0; j < nn; ++j)
            cp[c * nn + j] = std::strtoll(t[j].c_str(), nullptr, 10) - 1;
    }

    // Assemble.
    Mesh mesh;
    int xi = -1, yi = -1, zi = -1;
    for (std::size_t k = 0; k < variables.size(); ++k) {
        std::string v = tecplot_upper(variables[k]);
        if (v == "X")
            xi = (int)k;
        else if (v == "Y")
            yi = (int)k;
        else if (v == "Z")
            zi = (int)k;
    }
    std::size_t ndim = (zi >= 0) ? 3 : 2;
    NDArray pts(DType::Float64, {num_nodes, ndim});
    double* pp = pts.As<double>();
    for (std::size_t r = 0; r < num_nodes; ++r) {
        pp[r * ndim + 0] = cols[xi][r];
        pp[r * ndim + 1] = cols[yi][r];
        if (zi >= 0)
            pp[r * ndim + 2] = cols[zi][r];
    }
    mesh.AssignPoints(std::move(pts));
    for (std::size_t k = 0; k < variables.size(); ++k) {
        if ((int)k == xi || (int)k == yi || (int)k == zi)
            continue;
        NDArray arr(DType::Float64, {cols[k].size()});
        std::memcpy(arr.Data(), cols[k].data(), cols[k].size() * sizeof(double));
        if (cell_centered[k]) {
            std::vector<NDArray> blk;
            blk.push_back(std::move(arr));
            mesh.AddCellData(variables[k], std::move(blk));
        } else {
            mesh.AddPointData(variables[k], std::move(arr));
        }
    }
    mesh.AddCellBlock(mtype, std::move(celldata));
    return mesh;
}

void write_tecplot(const std::string& rPath, const Mesh& rMesh) {
    // Gather supported cell blocks; require a single unique type.
    std::vector<std::size_t> blocks;
    std::set<std::string> types;
    for (std::size_t i = 0; i < rMesh.NumCellBlocks(); ++i) {
        const auto cb = rMesh.Cells(i);
        if (!meshio_to_tecplot(cb.Type()).empty()) {
            blocks.push_back(i);
            types.insert(cb.Type());
        }
    }
    if (types.size() != 1)
        throw WriteError("C++ Tecplot writer supports a single cell type");
    std::string mtype = *types.begin();
    std::string ztype = meshio_to_tecplot(mtype);
    const std::vector<int>& order = tecplot_order(mtype);

    std::ofstream os(rPath);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);

    const std::size_t dim = rMesh.PointDim();
    const std::size_t num_nodes = rMesh.NumPoints();
    std::size_t num_cells = 0;
    for (std::size_t b : blocks)
        num_cells += rMesh.Cells(b).NumCells();

    // Variables + data columns.
    const NDArray& points = rMesh.Points();
    std::vector<std::string> variables = {"X", "Y"};
    std::vector<std::vector<double>> data;
    auto push_point_col = [&](std::size_t comp) {
        std::vector<double> col(num_nodes);
        for (std::size_t r = 0; r < num_nodes; ++r)
            col[r] = detail::read_double(points, r * dim + comp);
        data.push_back(std::move(col));
    };
    push_point_col(0);
    push_point_col(1);
    int varrange0 = 3, varrange1 = 0;
    if (dim == 3) {
        variables.push_back("Z");
        push_point_col(2);
        varrange0 += 1;
    }

    for (const auto& k : rMesh.PointDataNames()) {
        std::string ku = tecplot_upper(k);
        if (ku == "X" || ku == "Y" || ku == "Z")
            continue;
        const NDArray& v = rMesh.PointData(k);
        std::size_t ncomp = v.Shape().size() >= 2 ? v.Shape()[1] : 1;
        for (std::size_t c = 0; c < ncomp; ++c) {
            variables.push_back(ncomp == 1 ? k : k + "_" + std::to_string(c));
            std::vector<double> col(num_nodes);
            for (std::size_t r = 0; r < num_nodes; ++r)
                col[r] = detail::read_double(v, r * ncomp + c);
            data.push_back(std::move(col));
            varrange0 += 1;
        }
    }
    bool have_cell_data = false;
    varrange1 = varrange0 - 1;
    for (const auto& k : rMesh.CellDataNames()) {
        std::string ku = tecplot_upper(k);
        if (ku == "X" || ku == "Y" || ku == "Z")
            continue;
        if (rMesh.CellDataNumBlocks(k) == 0)
            continue;
        // concatenate the (single-type) blocks
        const NDArray& first = rMesh.CellData(k, 0);
        std::size_t ncomp = first.Shape().size() >= 2 ? first.Shape()[1] : 1;
        for (std::size_t c = 0; c < ncomp; ++c) {
            variables.push_back(ncomp == 1 ? k : k + "_" + std::to_string(c));
            std::vector<double> col;
            for (std::size_t b : blocks) {
                const NDArray& vv = rMesh.CellData(k, b);
                for (std::size_t r = 0; r < vv.Shape()[0]; ++r)
                    col.push_back(detail::read_double(vv, r * ncomp + c));
            }
            data.push_back(std::move(col));
            varrange1 += 1;
            have_cell_data = true;
        }
    }

    os << "TITLE = \"Written by meshio++ (C++ core)\"\n";
    os << "VARIABLES = ";
    for (std::size_t k = 0; k < variables.size(); ++k)
        os << (k ? ", " : "") << "\"" << variables[k] << "\"";
    os << "\n";
    os << "ZONE NODES = " << num_nodes << ", ELEMENTS = " << num_cells << ",\n";
    os << "DATAPACKING = BLOCK, ZONETYPE = " << ztype;
    if (have_cell_data && varrange0 <= varrange1) {
        os << ",\n";
        std::string r = (varrange0 == varrange1)
                            ? std::to_string(varrange0)
                            : std::to_string(varrange0) + "-" + std::to_string(varrange1);
        os << "VARLOCATION = ([" << r << "] = CELLCENTERED)\n";
    } else {
        os << "\n";
    }

    char buf[40];
    for (const auto& col : data) {
        for (std::size_t i = 0; i < col.size(); ++i) {
            std::snprintf(buf, sizeof(buf), "%.17g", col[i]);
            os << buf << ((i + 1) % 20 == 0 || i + 1 == col.size() ? '\n' : ' ');
        }
        if (col.empty())
            os << "\n";
    }

    for (std::size_t b : blocks) {
        const auto cb = rMesh.Cells(b);
        const NDArray& conn = cb.Conn();
        std::size_t k = conn.Shape().size() >= 2 ? conn.Shape()[1] : 1;
        for (std::size_t r = 0; r < cb.NumCells(); ++r) {
            for (std::size_t j = 0; j < order.size(); ++j) {
                std::size_t src = static_cast<std::size_t>(order[j]);
                if (src >= k)
                    src = k - 1;
                os << (detail::read_int(conn, r * k + src) + 1)
                   << (j + 1 == order.size() ? '\n' : ' ');
            }
        }
    }
}

}  // namespace meshioplusplus
