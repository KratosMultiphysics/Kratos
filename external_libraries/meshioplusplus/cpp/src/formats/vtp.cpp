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
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/formats/vtp.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/vtk_xml.hpp"
#include "meshioplusplus/detail/vtu_binary.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {

namespace {

using detail::cols;
using detail::read_double;
using detail::read_int;
using detail::vtu_ascii_double;
using detail::vtu_ascii_ndarray;
using detail::vtu_type_str;

// PolyData section a cell block belongs to (VTK's canonical cell order).
enum class VtpSection { Verts = 0, Lines = 1, Polys = 2 };

VtpSection vtp_section_of(const Mesh::CellView& rCb) {
    const std::string& type = rCb.Type();
    if (rCb.IsPolyhedron())
        throw WriteError("VTP: PolyData cannot hold polyhedron cells");
    if (type == "vertex")
        return VtpSection::Verts;
    if (type == "line")
        return VtpSection::Lines;
    if (type == "triangle" || type == "quad" || type == "polygon")
        return VtpSection::Polys;
    throw WriteError("VTP: PolyData cannot hold '" + type + "' cells");
}

// Flat connectivity + VTK end-offsets of one section's blocks.
struct VtpSectionData {
    std::vector<std::int64_t> mConn;
    std::vector<std::int64_t> mOffsets;
};

void vtp_append_block(VtpSectionData& rSec, const Mesh::CellView& rCb) {
    if (rCb.IsRagged()) {
        for (std::size_t r = 0; r < rCb.NumCells(); ++r) {
            const std::size_t sz = rCb.RowSize(r);
            const std::int64_t* row = rCb.Row(r);
            rSec.mConn.insert(rSec.mConn.end(), row, row + sz);
            rSec.mOffsets.push_back(static_cast<std::int64_t>(rSec.mConn.size()));
        }
        return;
    }
    const NDArray& conn = rCb.Conn();
    const std::size_t k = cols(conn);
    for (std::size_t r = 0; r < rCb.NumCells(); ++r) {
        for (std::size_t j = 0; j < k; ++j)
            rSec.mConn.push_back(read_int(conn, r * k + j));
        rSec.mOffsets.push_back(static_cast<std::int64_t>(rSec.mConn.size()));
    }
}

}  // namespace

void write_vtp(const std::string& rPath, const Mesh& rMesh, bool binary, bool zlib) {
    // Classify blocks and build the VTK canonical order (Verts, Lines, Polys)
    // as a stable partition of the mesh's block order.
    const std::size_t nblocks = rMesh.NumCellBlocks();
    std::vector<VtpSection> sections(nblocks);
    for (std::size_t bi = 0; bi < nblocks; ++bi)
        sections[bi] = vtp_section_of(rMesh.Cells(bi));
    std::vector<std::size_t> block_order;
    block_order.reserve(nblocks);
    VtpSectionData verts, lines, polys;
    for (VtpSection want : {VtpSection::Verts, VtpSection::Lines, VtpSection::Polys})
        for (std::size_t bi = 0; bi < nblocks; ++bi) {
            if (sections[bi] != want)
                continue;
            block_order.push_back(bi);
            VtpSectionData& sec = want == VtpSection::Verts   ? verts
                                  : want == VtpSection::Lines ? lines
                                                              : polys;
            vtp_append_block(sec, rMesh.Cells(bi));
        }

    for (std::size_t i = 0; i < verts.mOffsets.size(); ++i)
        if (verts.mOffsets[i] - (i == 0 ? 0 : verts.mOffsets[i - 1]) != 1)
            throw WriteError("VTP: vertex cells must have exactly one node");

    std::ofstream os(rPath, std::ios::binary);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);

    const NDArray& points = rMesh.Points();
    const std::size_t num_points = rMesh.NumPoints();
    const std::size_t dim = rMesh.PointDim();
    const std::size_t pt_isz = dtype_size(points.Dtype());

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

    os << "<?xml version=\"1.0\"?>\n";
    os << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\"";
    if (binary && zlib)
        os << " compressor=\"vtkZLibDataCompressor\"";
    os << ">\n";
    os << "<!--This file was created by meshio++ (C++ core)-->\n";
    os << "<PolyData>\n";
    os << "<Piece NumberOfPoints=\"" << num_points << "\" NumberOfVerts=\"" << verts.mOffsets.size()
       << "\" NumberOfLines=\"" << lines.mOffsets.size()
       << "\" NumberOfStrips=\"0\" NumberOfPolys=\"" << polys.mOffsets.size() << "\">\n";

    // Points (3 components; pad 2D with zero z) — same layout as the VTU writer.
    os << "<Points>\n";
    da_header(vtu_type_str(points.Dtype()), "Points", 3);
    if (binary) {
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

    auto emit_section = [&](const char* tag, const VtpSectionData& rSec) {
        if (rSec.mOffsets.empty())
            return;
        os << "<" << tag << ">\n";
        emit_i64("connectivity", rSec.mConn);
        emit_i64("offsets", rSec.mOffsets);
        os << "</" << tag << ">\n";
    };
    emit_section("Verts", verts);
    emit_section("Lines", lines);
    emit_section("Polys", polys);

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
        // Cell data follows the reordered (Verts, Lines, Polys) block order.
        os << "<CellData>\n";
        for (const auto& name : rMesh.CellDataNames()) {
            const std::size_t ndblocks = rMesh.CellDataNumBlocks(name);
            if (ndblocks == 0)
                continue;
            const NDArray& first = rMesh.CellData(name, 0);
            int ncomp = (first.Shape().size() == 2) ? static_cast<int>(cols(first)) : 0;
            da_header(vtu_type_str(first.Dtype()), name, ncomp);
            if (binary) {
                std::vector<unsigned char> buf;
                for (std::size_t bi : block_order) {
                    if (bi >= ndblocks)
                        continue;
                    const NDArray& blk = rMesh.CellData(name, bi);
                    const unsigned char* p = reinterpret_cast<const unsigned char*>(blk.Data());
                    buf.insert(buf.end(), p, p + blk.Nbytes());
                }
                emit_bin(buf.data(), buf.size());
            } else {
                for (std::size_t bi : block_order) {
                    if (bi >= ndblocks)
                        continue;
                    vtu_ascii_ndarray(os, rMesh.CellData(name, bi));
                }
            }
            os << "</DataArray>\n";
        }
        os << "</CellData>\n";
    }

    os << "</Piece>\n</PolyData>\n</VTKFile>\n";
}

}  // namespace meshioplusplus
