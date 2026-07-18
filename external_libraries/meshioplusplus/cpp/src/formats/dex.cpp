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
#include "meshioplusplus/formats/dex.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus {

namespace {

constexpr int kDim = 3;  // DEX coordinates are always x y z

// Extract `KEY = value` from a header string.
std::string header_value(const std::string& rText, const std::string& rKey) {
    std::size_t p = rText.find(rKey);
    if (p == std::string::npos)
        return {};
    p = rText.find('=', p);
    if (p == std::string::npos)
        return {};
    ++p;
    while (p < rText.size() && (rText[p] == ' ' || rText[p] == '\t'))
        ++p;
    std::size_t e = p;
    while (e < rText.size() && rText[e] != ' ' && rText[e] != '\t' && rText[e] != '#' &&
           rText[e] != '\r' && rText[e] != '\n')
        ++e;
    return rText.substr(p, e - p);
}

}  // namespace

Mesh read_dex(const std::string& rPath) {
    std::ifstream in(rPath, std::ios::binary);
    if (!in)
        throw ReadError("Could not open file: " + rPath);
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(in, line)) {
        // Files written in text mode on Windows use CRLF; the file is opened
        // in binary mode here (no newline translation) and std::getline only
        // splits on '\n', so strip a trailing '\r' explicitly.
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        lines.push_back(line);
    }

    // header = first two non-empty lines
    std::vector<std::string> header;
    std::size_t body_start = 0;
    for (std::size_t i = 0; i < lines.size(); ++i) {
        if (lines[i].find_first_not_of(" \t\r") != std::string::npos)
            header.push_back(lines[i]);
        if (header.size() == 2) {
            body_start = i + 1;
            break;
        }
    }
    std::string head = header.empty() ? std::string() : header[0];
    if (header.size() > 1)
        head += " " + header[1];

    std::string field = header_value(head, "FORMULA");
    if (field.empty())
        field = "dex:field";
    std::string ncomp_s = header_value(head, "NB_COMP");
    std::string npoint_s = header_value(head, "NB_POINT");
    int ncomp = ncomp_s.empty() ? 1 : std::atoi(ncomp_s.c_str());
    if (ncomp < 1)
        ncomp = 1;
    std::size_t npoint =
        npoint_s.empty() ? 0 : static_cast<std::size_t>(std::atoll(npoint_s.c_str()));

    std::vector<std::vector<double>> rows;
    for (std::size_t i = body_start; i < lines.size(); ++i) {
        std::istringstream iss(lines[i]);
        std::vector<double> r;
        std::string tok;
        while (iss >> tok) {
            for (char& c : tok)
                if (c == 'D' || c == 'd')
                    c = 'E';
            r.push_back(std::strtod(tok.c_str(), nullptr));
        }
        if (!r.empty())
            rows.push_back(std::move(r));
        if (npoint && rows.size() >= npoint)
            break;
    }
    std::size_t n = rows.size();

    Mesh mesh;
    NDArray pts(DType::Float64, {n, static_cast<std::size_t>(kDim)});
    for (std::size_t r = 0; r < n; ++r)
        for (int c = 0; c < kDim; ++c)
            pts.As<double>()[r * kDim + c] =
                c < static_cast<int>(rows[r].size()) ? rows[r][c] : 0.0;
    mesh.AssignPoints(std::move(pts));

    std::size_t nc = static_cast<std::size_t>(ncomp);
    NDArray vals = nc == 1 ? NDArray(DType::Float64, {n}) : NDArray(DType::Float64, {n, nc});
    for (std::size_t r = 0; r < n; ++r)
        for (std::size_t c = 0; c < nc; ++c) {
            std::size_t src = kDim + c;
            vals.As<double>()[r * nc + c] = src < rows[r].size() ? rows[r][src] : 0.0;
        }
    mesh.AddPointData(field, std::move(vals));
    return mesh;
}

void write_dex(const std::string& rPath, const Mesh& rMesh) {
    std::ofstream f(rPath, std::ios::binary);
    if (!f)
        throw WriteError("Could not open file for writing: " + rPath);

    auto names = rMesh.PointDataNames();
    if (names.empty())
        throw WriteError("DEX write needs a nodal field in point_data");
    const std::string& field = names.front();
    const NDArray& arr = rMesh.PointData(field);

    const std::size_t n = rMesh.NumPoints();
    const std::size_t pdim = rMesh.PointDim();
    const NDArray& points = rMesh.Points();
    std::size_t ncomp = n ? arr.Size() / n : 0;
    if (ncomp == 0)
        ncomp = 1;

    f << "# NAME = PIECE FORMULA = " << field << "\n";
    f << "NB_REAL = 1 NB_COMP = " << ncomp << " NB_POINT = " << n << " #\n";
    char buf[64];
    for (std::size_t r = 0; r < n; ++r) {
        for (int c = 0; c < kDim; ++c) {
            double v = c < static_cast<int>(pdim) ? detail::read_double(points, r * pdim + c) : 0.0;
            std::snprintf(buf, sizeof(buf), "%.16g", v);
            f << buf << (c + 1 < kDim ? " " : "");
        }
        for (std::size_t c = 0; c < ncomp; ++c) {
            std::snprintf(buf, sizeof(buf), " %.16g", detail::read_double(arr, r * ncomp + c));
            f << buf;
        }
        f << "\n";
    }
}

}  // namespace meshioplusplus
