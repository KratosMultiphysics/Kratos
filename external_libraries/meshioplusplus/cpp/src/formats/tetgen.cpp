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
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/tetgen.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus {

namespace {

// Split "<stem>.node" / "<stem>.ele" into the two sibling paths.
std::pair<std::string, std::string> node_ele_paths(const std::string& rPath, bool& rOk) {
    std::size_t dot = rPath.find_last_of('.');
    rOk = false;
    if (dot == std::string::npos)
        return {"", ""};
    std::string suffix = rPath.substr(dot);
    std::string stem = rPath.substr(0, dot);
    if (suffix == ".node" || suffix == ".ele") {
        rOk = true;
        return {stem + ".node", stem + ".ele"};
    }
    return {"", ""};
}

// First non-comment, non-blank line is the header; remaining non-comment
// tokens (whitespace-separated, across lines) are the data stream.
struct Parsed {
    std::vector<std::string> mHeader;
    std::vector<std::string> mData;
};

Parsed parse_file(const std::string& rPath) {
    std::ifstream in(rPath, std::ios::binary);
    if (!in)
        throw ReadError("Could not open file: " + rPath);
    Parsed p;
    bool have_header = false;
    std::string line;
    while (std::getline(in, line)) {
        // trim leading whitespace
        std::size_t s = 0;
        while (s < line.size() && std::isspace(static_cast<unsigned char>(line[s])))
            ++s;
        if (s >= line.size() || line[s] == '#')
            continue;
        std::istringstream iss(line);
        std::string tok;
        if (!have_header) {
            while (iss >> tok)
                p.mHeader.push_back(tok);
            have_header = true;
        } else {
            while (iss >> tok)
                p.mData.push_back(tok);
        }
    }
    if (!have_header)
        throw ReadError("TetGen: missing header line in " + rPath);
    return p;
}

}  // namespace

Mesh read_tetgen(const std::string& rPath) {
    bool ok = false;
    auto paths = node_ele_paths(rPath, ok);
    if (!ok)
        throw ReadError("TetGen: expected a .node or .ele file");
    const std::string& node_path = paths.first;
    const std::string& ele_path = paths.second;

    Mesh mesh;

    // ---- nodes ----
    Parsed nf = parse_file(node_path);
    if (nf.mHeader.size() < 4)
        throw ReadError("TetGen: malformed .node header");
    std::int64_t npoints = std::strtoll(nf.mHeader[0].c_str(), nullptr, 10);
    int dim = static_cast<int>(std::strtoll(nf.mHeader[1].c_str(), nullptr, 10));
    int num_attrs = static_cast<int>(std::strtoll(nf.mHeader[2].c_str(), nullptr, 10));
    int num_bmarkers = static_cast<int>(std::strtoll(nf.mHeader[3].c_str(), nullptr, 10));
    if (dim != 3)
        throw ReadError("TetGen: need 3D points");

    const int ncol = 4 + num_attrs + num_bmarkers;
    if (static_cast<std::int64_t>(nf.mData.size()) != npoints * ncol)
        throw ReadError("TetGen: .node data size mismatch");

    auto at = [&](std::int64_t r, int c) -> double {
        return std::strtod(nf.mData[r * ncol + c].c_str(), nullptr);
    };

    std::int64_t node_index_base = npoints > 0 ? static_cast<std::int64_t>(at(0, 0)) : 0;
    for (std::int64_t i = 0; i < npoints; ++i) {
        if (static_cast<std::int64_t>(at(i, 0)) != node_index_base + i)
            throw ReadError("TetGen: nodes not numbered consecutively");
    }

    NDArray pts(DType::Float64, {static_cast<std::size_t>(npoints), 3});
    double* pp = pts.As<double>();
    for (std::int64_t i = 0; i < npoints; ++i)
        for (int c = 0; c < 3; ++c)
            pp[i * 3 + c] = at(i, 1 + c);
    mesh.AssignPoints(std::move(pts));

    // point attributes
    for (int k = 0; k < num_attrs; ++k) {
        NDArray a(DType::Float64, {static_cast<std::size_t>(npoints)});
        for (std::int64_t i = 0; i < npoints; ++i)
            a.As<double>()[i] = at(i, 4 + k);
        mesh.AddPointData("tetgen:attr" + std::to_string(k + 1), std::move(a));
    }
    // boundary markers: tetgen:ref, tetgen:ref2, ...
    for (int k = 0; k < num_bmarkers; ++k) {
        std::string name = "tetgen:ref" + (k == 0 ? std::string() : std::to_string(k + 1));
        NDArray a(DType::Float64, {static_cast<std::size_t>(npoints)});
        for (std::int64_t i = 0; i < npoints; ++i)
            a.As<double>()[i] = at(i, 4 + num_attrs + k);
        mesh.AddPointData(std::move(name), std::move(a));
    }

    // ---- elements ----
    Parsed ef = parse_file(ele_path);
    if (ef.mHeader.size() < 3)
        throw ReadError("TetGen: malformed .ele header");
    std::int64_t num_tets = std::strtoll(ef.mHeader[0].c_str(), nullptr, 10);
    int npt = static_cast<int>(std::strtoll(ef.mHeader[1].c_str(), nullptr, 10));
    int ele_attrs = static_cast<int>(std::strtoll(ef.mHeader[2].c_str(), nullptr, 10));
    if (npt != 4)
        throw ReadError("TetGen: only 4-node tetrahedra supported");

    const int ecol = 5 + ele_attrs;
    if (static_cast<std::int64_t>(ef.mData.size()) != num_tets * ecol)
        throw ReadError("TetGen: .ele data size mismatch");

    auto eat = [&](std::int64_t r, int c) -> std::int64_t {
        return std::strtoll(ef.mData[r * ecol + c].c_str(), nullptr, 10);
    };

    NDArray cells(DType::Int64, {static_cast<std::size_t>(num_tets), 4});
    std::int64_t* cp = cells.As<std::int64_t>();
    for (std::int64_t i = 0; i < num_tets; ++i)
        for (int c = 0; c < 4; ++c)
            cp[i * 4 + c] = eat(i, 1 + c) - node_index_base;
    mesh.AddCellBlock("tetra", std::move(cells));

    // region attributes: tetgen:ref, tetgen:ref2, ...
    for (int k = 0; k < ele_attrs; ++k) {
        std::string name = "tetgen:ref" + (k == 0 ? std::string() : std::to_string(k + 1));
        NDArray a(DType::Int64, {static_cast<std::size_t>(num_tets)});
        for (std::int64_t i = 0; i < num_tets; ++i)
            a.As<std::int64_t>()[i] = eat(i, 5 + k);
        std::vector<NDArray> blocks;
        blocks.push_back(std::move(a));
        mesh.AddCellData(std::move(name), std::move(blocks));
    }

    return mesh;
}

namespace {

// Write a marker/ref value: integral values as integers, else %.16e.
void write_value(std::ostream& rOs, double v) {
    double r = std::nearbyint(v);
    if (v == r && std::fabs(v) < 9.2e18) {
        rOs << static_cast<std::int64_t>(r);
    } else {
        char buf[40];
        std::snprintf(buf, sizeof(buf), "%.16e", v);
        rOs << buf;
    }
}

}  // namespace

void write_tetgen(const std::string& rPath, const Mesh& rMesh) {
    bool ok = false;
    auto paths = node_ele_paths(rPath, ok);
    if (!ok)
        throw WriteError("TetGen: must specify a .node or .ele file");
    const std::string& node_path = paths.first;
    const std::string& ele_path = paths.second;

    const NDArray& points = rMesh.Points();
    const std::size_t ncols = rMesh.PointDim();
    if (ncols != 3)
        throw WriteError("TetGen: can only write 3D points");

    const std::int64_t npoints = static_cast<std::int64_t>(rMesh.NumPoints());

    // ---- node file ----
    {
        std::ofstream fh(node_path, std::ios::binary);
        if (!fh)
            throw WriteError("Could not open file for writing: " + node_path);

        // Split point_data into one ref key and the remaining attribute keys,
        // mirroring meshioplusplus.tetgen.write.
        std::vector<std::string> attr_keys =
            rMesh.PointDataNames();  // sorted: deterministic column order
        std::vector<std::string> ref_keys;
        if (!attr_keys.empty()) {
            for (const auto& k : attr_keys)
                if (k.find(":ref") != std::string::npos) {
                    ref_keys.push_back(k);
                    break;
                }
            if (!ref_keys.empty()) {
                attr_keys.erase(std::remove(attr_keys.begin(), attr_keys.end(), ref_keys[0]),
                                attr_keys.end());
            } else {
                ref_keys.push_back(attr_keys.front());
                attr_keys.erase(attr_keys.begin());
            }
        }
        const std::size_t nattr = attr_keys.size();
        const std::size_t nref = ref_keys.size();

        fh << "# This file was created by meshio++ (C++ core)\n";
        if (nattr + nref > 0) {
            fh << "# attribute and marker names: ";
            bool first = true;
            for (const auto& k : attr_keys) {
                fh << (first ? "" : ", ") << k;
                first = false;
            }
            for (const auto& k : ref_keys) {
                fh << (first ? "" : ", ") << k;
                first = false;
            }
            fh << "\n";
        }
        fh << npoints << " 3 " << nattr << " " << nref << "\n";

        char fbuf[40];
        for (std::int64_t i = 0; i < npoints; ++i) {
            fh << i;
            for (int c = 0; c < 3; ++c) {
                std::snprintf(fbuf, sizeof(fbuf), "%.16e",
                              detail::read_double(points, i * 3 + c));
                fh << " " << fbuf;
            }
            for (const auto& k : attr_keys) {
                std::snprintf(fbuf, sizeof(fbuf), "%.16e",
                              detail::read_double(rMesh.PointData(k), i));
                fh << " " << fbuf;
            }
            for (const auto& k : ref_keys) {
                fh << " ";
                write_value(fh, detail::read_double(rMesh.PointData(k), i));
            }
            fh << "\n";
        }
    }

    // ---- ele file ----
    {
        std::ofstream fh(ele_path, std::ios::binary);
        if (!fh)
            throw WriteError("Could not open file for writing: " + ele_path);

        // Cell-data attribute keys, with the first ":ref" key moved to front.
        std::vector<std::string> attr_keys =
            rMesh.CellDataNames();  // sorted: deterministic column order
        if (!attr_keys.empty()) {
            std::string ref;
            for (const auto& k : attr_keys)
                if (k.find(":ref") != std::string::npos) {
                    ref = k;
                    break;
                }
            if (!ref.empty()) {
                attr_keys.erase(std::remove(attr_keys.begin(), attr_keys.end(), ref),
                                attr_keys.end());
                attr_keys.insert(attr_keys.begin(), ref);
            }
        }
        const std::size_t nattr = attr_keys.size();

        fh << "# This file was created by meshio++ (C++ core)\n";
        if (nattr > 0) {
            fh << "# attribute names: ";
            bool first = true;
            for (const auto& k : attr_keys) {
                fh << (first ? "" : ", ") << k;
                first = false;
            }
            fh << "\n";
        }

        for (std::size_t ci = 0; ci < rMesh.NumCellBlocks(); ++ci) {
            const auto cb = rMesh.Cells(ci);
            if (cb.Type() != "tetra")
                continue;
            const NDArray& conn = cb.Conn();
            std::int64_t n = detail::rows(conn);
            fh << n << " 4 " << nattr << "\n";
            for (std::int64_t i = 0; i < n; ++i) {
                fh << i;
                for (int c = 0; c < 4; ++c)
                    fh << " " << detail::read_int(conn, i * 4 + c);
                for (const auto& k : attr_keys) {
                    if (ci < rMesh.CellDataNumBlocks(k))
                        fh << " " << detail::read_int(rMesh.CellData(k, ci), i);
                    else
                        fh << " 0";
                }
                fh << "\n";
            }
        }
    }
}

}  // namespace meshioplusplus
