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
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/formats/ip.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus {

namespace {

std::string ip_strip(const std::string& s) {
    std::size_t a = s.find_first_not_of(" \t\r");
    std::size_t b = s.find_last_not_of(" \t\r");
    return a == std::string::npos ? std::string() : s.substr(a, b - a + 1);
}

}  // namespace

Mesh read_ip(const std::string& rPath) {
    std::ifstream in(rPath, std::ios::binary);
    if (!in)
        throw ReadError("Could not open file: " + rPath);
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(in, line))
        lines.push_back(line);

    // header: first four non-empty lines -> version, dim, npoint, ncomp
    std::vector<int> ints;
    std::size_t idx = 0;
    while (ints.size() < 4 && idx < lines.size()) {
        std::string s = ip_strip(lines[idx++]);
        if (!s.empty()) {
            std::istringstream iss(s);
            int v;
            iss >> v;
            ints.push_back(v);
        }
    }
    if (ints.size() < 4)
        throw ReadError("IP: malformed header");
    int dim = ints[1];
    std::size_t npoint = static_cast<std::size_t>(ints[2]);
    int ncomp = ints[3];

    std::vector<std::string> names;
    while (static_cast<int>(names.size()) < ncomp && idx < lines.size()) {
        std::string s = ip_strip(lines[idx++]);
        if (!s.empty())
            names.push_back(s);
    }

    // remaining tokens (treat '(' ')' as whitespace) form (dim + ncomp)
    // column-major sections of npoint reals each.
    std::vector<double> flat;
    for (; idx < lines.size(); ++idx) {
        std::string s = lines[idx];
        for (char& c : s)
            if (c == '(' || c == ')')
                c = ' ';
            else if (c == 'D' || c == 'd')
                c = 'E';
        std::istringstream iss(s);
        std::string tok;
        while (iss >> tok)
            flat.push_back(std::strtod(tok.c_str(), nullptr));
    }

    std::size_t nsec = static_cast<std::size_t>(dim + ncomp);
    auto section = [&](std::size_t s, std::size_t i) -> double {
        std::size_t p = s * npoint + i;
        return p < flat.size() ? flat[p] : 0.0;
    };

    Mesh mesh;
    NDArray pts(DType::Float64, {npoint, static_cast<std::size_t>(dim)});
    for (std::size_t i = 0; i < npoint; ++i)
        for (int d = 0; d < dim; ++d)
            pts.As<double>()[i * dim + d] = section(static_cast<std::size_t>(d), i);
    mesh.AssignPoints(std::move(pts));

    for (int c = 0; c < ncomp; ++c) {
        NDArray vals(DType::Float64, {npoint});
        for (std::size_t i = 0; i < npoint; ++i)
            vals.As<double>()[i] = section(static_cast<std::size_t>(dim + c), i);
        mesh.AddPointData(names[c], std::move(vals));
    }
    (void)nsec;
    return mesh;
}

void write_ip(const std::string& rPath, const Mesh& rMesh) {
    std::ofstream f(rPath, std::ios::binary);
    if (!f)
        throw WriteError("Could not open file for writing: " + rPath);

    const std::size_t n = rMesh.NumPoints();
    const std::size_t dim = rMesh.PointDim();
    const NDArray& points = rMesh.Points();

    // flatten point_data into scalar component columns
    std::vector<std::string> names;
    std::vector<std::vector<double>> columns;
    for (const auto& name : rMesh.PointDataNames()) {
        const NDArray& arr = rMesh.PointData(name);
        std::size_t nc = n ? arr.Size() / n : 0;
        if (nc <= 1) {
            names.push_back(name);
            std::vector<double> col(n);
            for (std::size_t i = 0; i < n; ++i)
                col[i] = detail::read_double(arr, i);
            columns.push_back(std::move(col));
        } else {
            for (std::size_t c = 0; c < nc; ++c) {
                names.push_back(name + "_" + std::to_string(c));
                std::vector<double> col(n);
                for (std::size_t i = 0; i < n; ++i)
                    col[i] = detail::read_double(arr, i * nc + c);
                columns.push_back(std::move(col));
            }
        }
    }

    f << "3\n" << dim << "\n" << n << "\n" << columns.size() << "\n";
    for (const auto& name : names)
        f << name << "\n";
    char buf[64];
    auto write_section = [&](const std::vector<double>& col) {
        f << "(";
        for (std::size_t i = 0; i < col.size(); ++i) {
            std::snprintf(buf, sizeof(buf), "%.16g", col[i]);
            f << (i ? "\n" : "") << buf;
        }
        f << "\n)\n";
    };
    for (std::size_t d = 0; d < dim; ++d) {
        std::vector<double> col(n);
        for (std::size_t i = 0; i < n; ++i)
            col[i] = detail::read_double(points, i * dim + d);
        write_section(col);
    }
    for (const auto& col : columns)
        write_section(col);
}

}  // namespace meshioplusplus
