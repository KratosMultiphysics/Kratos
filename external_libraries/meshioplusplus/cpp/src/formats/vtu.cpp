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
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/formats/vtu.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/vtk_xml.hpp"
#include "meshioplusplus/detail/vtu_binary.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/vtk_common.hpp"

namespace meshioplusplus {

namespace {

using detail::cols;
using detail::read_double;
using detail::read_int;
using detail::vtu_ascii_double;
using detail::vtu_ascii_ndarray;
using detail::vtu_type_str;

}  // namespace

void write_vtu(const std::string& rPath, const Mesh& rMesh, bool binary, bool zlib) {
    for (const auto cb : rMesh.CellRange()) {
        if (cb.Type().rfind("polyhedron", 0) == 0)
            throw WriteError("C++ VTU writer does not support polyhedron cells");
    }

    std::ofstream os(rPath, std::ios::binary);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);

    const NDArray& points = rMesh.Points();
    const std::size_t num_points = rMesh.NumPoints();
    const std::size_t dim = rMesh.PointDim();
    const std::size_t pt_isz = dtype_size(points.Dtype());

    std::size_t total_cells = 0;
    for (const auto cb : rMesh.CellRange())
        total_cells += cb.NumCells();

    const char* fmt = binary ? "binary" : "ascii";

    auto da_header = [&](const char* type, const std::string& name, int ncomp) {
        os << "<DataArray type=\"" << type << "\" Name=\"" << name << "\"";
        if (ncomp > 0)
            os << " NumberOfComponents=\"" << ncomp << "\"";
        os << " format=\"" << fmt << "\">\n";
    };
    auto emit_bin = [&](const unsigned char* d, std::size_t n) {
        os << detail::vtu_encode_binary(d, n, zlib) << "\n";
    };

    os << "<?xml version=\"1.0\"?>\n";
    os << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
          "byte_order=\"LittleEndian\"";
    if (binary && zlib)
        os << " compressor=\"vtkZLibDataCompressor\"";
    os << ">\n";
    os << "<!--This file was created by meshio++ (C++ core)-->\n";
    os << "<UnstructuredGrid>\n";
    os << "<Piece NumberOfPoints=\"" << num_points << "\" NumberOfCells=\"" << total_cells
       << "\">\n";

    // Points (3 components; pad 2D with zero z).
    os << "<Points>\n";
    da_header(vtu_type_str(points.Dtype()), "Points", 3);
    if (binary) {
        // Pre-sized buffer (zero-filled -> the padded z stays 0), indexed
        // byte writes -> parallel over points.
        std::vector<unsigned char> buf(num_points * 3 * pt_isz, 0);
        const auto* src = reinterpret_cast<const unsigned char*>(points.Data());
        parallel_for(num_points, [&](std::size_t r) {
            for (std::size_t c = 0; c < dim && c < 3; ++c)
                std::memcpy(buf.data() + (r * 3 + c) * pt_isz, src + (r * dim + c) * pt_isz,
                            pt_isz);
        });
        emit_bin(buf.data(), buf.size());
    } else {
        for (std::size_t r = 0; r < num_points; ++r)
            for (std::size_t c = 0; c < 3; ++c)
                vtu_ascii_double(os, (c < dim) ? read_double(points, r * dim + c) : 0.0);
    }
    os << "</DataArray>\n</Points>\n";

    if (rMesh.NumCellBlocks() != 0) {
        // Build connectivity / offsets / types (Int64) into pre-sized arrays.
        // Per-block offsets are closed-form (conn_base + (r+1)*k), so rows are
        // independent and each block fills in parallel.
        const auto& tmap = meshio_to_vtk_type();
        std::size_t total_conn = 0, ncells = 0;
        for (const auto cb : rMesh.CellRange()) {
            total_conn += cb.NumCells() * cols(cb.Conn());
            ncells += cb.NumCells();
        }
        std::vector<std::int64_t> connectivity(total_conn), offsets(ncells), types(ncells);
        std::size_t conn_base = 0, cell_base = 0;
        for (const auto cb : rMesh.CellRange()) {
            const NDArray& conn = cb.Conn();
            const std::size_t nc = cb.NumCells();
            const std::size_t k = cols(conn);
            std::vector<int> order = meshio_to_vtk_order(cb.Type());
            auto it = tmap.find(cb.Type());
            if (it == tmap.end())
                throw WriteError("Unknown cell type for VTU: " + cb.Type());
            const std::int64_t vtk_type = it->second;
            parallel_for(nc, [&](std::size_t r) {
                for (std::size_t j = 0; j < k; ++j) {
                    std::size_t col = order.empty() ? j : static_cast<std::size_t>(order[j]);
                    connectivity[conn_base + r * k + j] = read_int(conn, r * k + col);
                }
                offsets[cell_base + r] = static_cast<std::int64_t>(conn_base + (r + 1) * k);
                types[cell_base + r] = vtk_type;
            });
            conn_base += nc * k;
            cell_base += nc;
        }

        auto emit_i64 = [&](const char* name, const std::vector<std::int64_t>& v) {
            da_header("Int64", name, 0);
            if (binary) {
                emit_bin(reinterpret_cast<const unsigned char*>(v.data()),
                         v.size() * sizeof(std::int64_t));
            } else {
                for (std::int64_t x : v)
                    os << x << '\n';
            }
            os << "</DataArray>\n";
        };

        os << "<Cells>\n";
        emit_i64("connectivity", connectivity);
        emit_i64("offsets", offsets);
        emit_i64("types", types);
        os << "</Cells>\n";
    }

    if (rMesh.NumPointData() != 0) {
        os << "<PointData>\n";
        for (const auto& name : rMesh.PointDataNames()) {
            const NDArray& d = rMesh.PointData(name);
            int ncomp = (d.Shape().size() == 2) ? static_cast<int>(cols(d)) : 0;
            da_header(vtu_type_str(d.Dtype()), name, ncomp);
            if (binary)
                emit_bin(reinterpret_cast<const unsigned char*>(d.Data()), d.Nbytes());
            else
                vtu_ascii_ndarray(os, d);
            os << "</DataArray>\n";
        }
        os << "</PointData>\n";
    }

    if (rMesh.NumCellData() != 0) {
        os << "<CellData>\n";
        for (const auto& name : rMesh.CellDataNames()) {
            const std::size_t nblocks = rMesh.CellDataNumBlocks(name);
            if (nblocks == 0)
                continue;
            const NDArray& first = rMesh.CellData(name, 0);
            int ncomp = (first.Shape().size() == 2) ? static_cast<int>(cols(first)) : 0;
            da_header(vtu_type_str(first.Dtype()), name, ncomp);
            if (binary) {
                std::vector<unsigned char> buf;
                for (std::size_t bi = 0; bi < nblocks; ++bi) {
                    const NDArray& blk = rMesh.CellData(name, bi);
                    const unsigned char* p = reinterpret_cast<const unsigned char*>(blk.Data());
                    buf.insert(buf.end(), p, p + blk.Nbytes());
                }
                emit_bin(buf.data(), buf.size());
            } else {
                for (std::size_t bi = 0; bi < nblocks; ++bi)
                    vtu_ascii_ndarray(os, rMesh.CellData(name, bi));
            }
            os << "</DataArray>\n";
        }
        os << "</CellData>\n";
    }

    os << "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
}

}  // namespace meshioplusplus
