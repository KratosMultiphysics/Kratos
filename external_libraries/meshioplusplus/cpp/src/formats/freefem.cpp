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
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/formats/freefem.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus {

namespace {

// Next non-blank line's whitespace tokens.
bool next_tokens(std::istream& rIn, std::vector<std::string>& rOut) {
    std::string line;
    while (std::getline(rIn, line)) {
        std::istringstream iss(line);
        std::string t;
        rOut.clear();
        while (iss >> t)
            rOut.push_back(t);
        if (!rOut.empty())
            return true;
    }
    return false;
}

}  // namespace

Mesh read_freefem(const std::string& rPath) {
    std::ifstream in(rPath, std::ios::binary);
    if (!in)
        throw ReadError("Could not open file: " + rPath);

    std::vector<std::string> tok;
    if (!next_tokens(in, tok) || tok.size() != 3)
        throw ReadError("FreeFem: expected a 3-integer header");
    const std::int64_t nver = std::strtoll(tok[0].c_str(), nullptr, 10);
    const std::int64_t n1 = std::strtoll(tok[1].c_str(), nullptr, 10);
    const std::int64_t n2 = std::strtoll(tok[2].c_str(), nullptr, 10);

    if (!next_tokens(in, tok))
        throw ReadError("FreeFem: missing vertices");
    const int dim = static_cast<int>(tok.size()) - 1;
    if (dim != 2 && dim != 3)
        throw ReadError("FreeFem: bad vertex dimension");

    Mesh mesh;
    NDArray pts(DType::Float64, {static_cast<std::size_t>(nver), static_cast<std::size_t>(dim)});
    NDArray pref(DType::Int64, {static_cast<std::size_t>(nver)});
    for (std::int64_t i = 0; i < nver; ++i) {
        if (i > 0 && !next_tokens(in, tok))
            throw ReadError("FreeFem: truncated vertices");
        for (int c = 0; c < dim; ++c)
            pts.As<double>()[i * dim + c] = std::strtod(tok[c].c_str(), nullptr);
        pref.As<std::int64_t>()[i] = std::strtoll(tok[dim].c_str(), nullptr, 10);
    }
    mesh.AssignPoints(std::move(pts));
    mesh.AddPointData("freefem:ref", std::move(pref));

    const char* t1 = dim == 2 ? "triangle" : "tetra";
    const int lnv1 = dim == 2 ? 3 : 4;
    const char* t2 = dim == 2 ? "line" : "triangle";
    const int lnv2 = dim == 2 ? 2 : 3;

    std::vector<NDArray> cell_refs;
    auto read_block = [&](std::int64_t n, const char* type, int lnv) {
        if (n <= 0)
            return;
        NDArray data(DType::Int64, {static_cast<std::size_t>(n), static_cast<std::size_t>(lnv)});
        NDArray ref(DType::Int64, {static_cast<std::size_t>(n)});
        for (std::int64_t k = 0; k < n; ++k) {
            if (!next_tokens(in, tok))
                throw ReadError("FreeFem: truncated elements");
            for (int j = 0; j < lnv; ++j)
                data.As<std::int64_t>()[k * lnv + j] =
                    std::strtoll(tok[j].c_str(), nullptr, 10) - 1;
            ref.As<std::int64_t>()[k] = std::strtoll(tok[lnv].c_str(), nullptr, 10);
        }
        mesh.AddCellBlock(type, std::move(data));
        cell_refs.push_back(std::move(ref));
    };
    read_block(n1, t1, lnv1);
    read_block(n2, t2, lnv2);
    if (!cell_refs.empty())
        mesh.AddCellData("freefem:ref", std::move(cell_refs));

    return mesh;
}

void write_freefem(const std::string& rPath, const Mesh& rMesh) {
    const int dim = static_cast<int>(rMesh.PointDim());
    if (dim != 2 && dim != 3)
        throw WriteError("FreeFem: can only write 2D/3D meshes");

    const std::string t1 = dim == 2 ? "triangle" : "tetra";
    const std::string t2 = dim == 2 ? "line" : "triangle";

    // Reject unsupported cell types so the shim falls back to Python (which
    // warns and skips). This keeps behaviour identical to the reference impl.
    for (const auto cb : rMesh.CellRange())
        if (cb.Type() != t1 && cb.Type() != t2)
            throw WriteError("FreeFem: unsupported cell type " + cb.Type());

    const bool has_ref = rMesh.HasCellData("freefem:ref");

    struct Row {
        Mesh::CellView mCb;
        const NDArray* mRef;
    };
    std::vector<Row> b1, b2;
    for (std::size_t i = 0; i < rMesh.NumCellBlocks(); ++i) {
        const NDArray* ref = (has_ref && i < rMesh.CellDataNumBlocks("freefem:ref"))
                                 ? &rMesh.CellData("freefem:ref", i)
                                 : nullptr;
        if (rMesh.Cells(i).Type() == t1)
            b1.push_back({rMesh.Cells(i), ref});
        else
            b2.push_back({rMesh.Cells(i), ref});
    }
    auto count = [](const std::vector<Row>& b) {
        std::size_t n = 0;
        for (const auto& r : b)
            n += r.mCb.NumCells();
        return n;
    };

    std::ofstream f(rPath, std::ios::binary);
    if (!f)
        throw WriteError("Could not open file for writing: " + rPath);

    const std::size_t nver = rMesh.NumPoints();
    f << nver << " " << count(b1) << " " << count(b2) << "\n";

    const NDArray* pref = rMesh.HasPointData("freefem:ref") ? &rMesh.PointData("freefem:ref")
                                                            : nullptr;

    const NDArray& points = rMesh.Points();
    char buf[32];
    for (std::size_t i = 0; i < nver; ++i) {
        for (int c = 0; c < dim; ++c) {
            std::snprintf(buf, sizeof(buf), "%.16e", detail::read_double(points, i * dim + c));
            f << buf << " ";
        }
        f << (pref ? detail::read_int(*pref, i) : 0) << "\n";
    }
    auto write_block = [&](const std::vector<Row>& b, int lnv) {
        for (const auto& r : b) {
            std::size_t n = r.mCb.NumCells();
            const NDArray& conn = r.mCb.Conn();
            std::size_t k = detail::cols(conn);
            for (std::size_t rr = 0; rr < n; ++rr) {
                for (int j = 0; j < lnv && static_cast<std::size_t>(j) < k; ++j)
                    f << (detail::read_int(conn, rr * k + j) + 1) << " ";
                f << (r.mRef ? detail::read_int(*r.mRef, rr) : 0) << "\n";
            }
        }
    };
    write_block(b1, dim == 2 ? 3 : 4);
    write_block(b2, dim == 2 ? 2 : 3);
}

}  // namespace meshioplusplus
