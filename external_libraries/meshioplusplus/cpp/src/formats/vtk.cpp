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
#include "meshioplusplus/formats/vtk.hpp"
#include "meshioplusplus/detail/byteswap.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/vtk_common.hpp"

namespace meshioplusplus {

namespace {

using detail::cols;
using detail::dispatch_dtype;
using detail::is_float_dtype;
using detail::read_double;
using detail::read_int;

const char* vtk_dtype_str(DType dt) {
    switch (dt) {
        case DType::Float32:
            return "float";
        case DType::Float64:
            return "double";
        case DType::Int8:
            return "vtktypeint8";
        case DType::Int16:
            return "vtktypeint16";
        case DType::Int32:
            return "vtktypeint32";
        case DType::Int64:
            return "vtktypeint64";
        case DType::UInt8:
            return "vtktypeuint8";
        case DType::UInt16:
            return "vtktypeuint16";
        case DType::UInt32:
            return "vtktypeuint32";
        case DType::UInt64:
            return "vtktypeuint64";
    }
    return "double";
}

void vtk_ascii_double(std::ostream& rOs, double v) {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.17g", v);
    rOs << buf;
}

// Byte-swap a whole array into a big-endian buffer (elements independent ->
// parallel), for a single os.write instead of per-element stream calls.
std::vector<unsigned char> be_buffer(const NDArray& rA) {
    const int isz = static_cast<int>(dtype_size(rA.Dtype()));
    const std::size_t n = rA.Size();
    const auto* src = reinterpret_cast<const char*>(rA.Data());
    std::vector<unsigned char> buf(n * static_cast<std::size_t>(isz));
    auto* dst = reinterpret_cast<char*>(buf.data());
    parallel_for_bw(n,
                    [&](std::size_t i) { detail::bswap_copy(dst + i * isz, src + i * isz, isz); });
    return buf;
}

// Byte-swap a typed vector into a big-endian buffer and emit it in one write
// (replaces per-element os.put stream calls for the CELLS/OFFSETS/... sections).
template <class T>
void write_be(std::ostream& rOs, const std::vector<T>& rV) {
    constexpr int isz = static_cast<int>(sizeof(T));
    std::vector<unsigned char> buf(rV.size() * sizeof(T));
    const auto* src = reinterpret_cast<const char*>(rV.data());
    auto* dst = reinterpret_cast<char*>(buf.data());
    parallel_for_bw(rV.size(), [&](std::size_t i) {
        detail::bswap_copy(dst + i * sizeof(T), src + i * sizeof(T), isz);
    });
    rOs.write(reinterpret_cast<const char*>(buf.data()), static_cast<std::streamsize>(buf.size()));
}

// Store the low `bytes` bytes of v into dst in big-endian order.
inline void be_store(unsigned char* pDst, std::uint64_t v, std::size_t bytes) {
    if (bytes == 8) {
        std::uint64_t be = detail::bswap64(v);
        std::memcpy(pDst, &be, 8);
    } else if (bytes == 4) {
        std::uint32_t be = detail::bswap32(static_cast<std::uint32_t>(v));
        std::memcpy(pDst, &be, 4);
    } else {
        for (std::size_t b = 0; b < bytes; ++b)
            pDst[b] = static_cast<unsigned char>(v >> (8 * (bytes - 1 - b)));
    }
}

// Fused gather + big-endian store of one connectivity block: reads each
// (optionally reordered) index and writes it as `bytes` big-endian bytes into
// `dst` — one pass, no intermediate typed buffer (halves the memory traffic vs
// building an int64/int32 array and byte-swapping it separately).
inline void gather_be(unsigned char* pDst, const NDArray& rData, std::size_t nc, std::size_t k,
                      const int* pOrd, std::size_t bytes) {
    dispatch_dtype(rData.Dtype(), [&]<class T>() {
        const T* src = rData.As<T>();
        parallel_for_bw(nc, [&](std::size_t r) {
            for (std::size_t j = 0; j < k; ++j) {
                std::size_t col = pOrd ? static_cast<std::size_t>(pOrd[j]) : j;
                auto v = static_cast<std::uint64_t>(static_cast<std::int64_t>(src[r * k + col]));
                be_store(pDst + (r * k + j) * bytes, v, bytes);
            }
        });
    });
}

void write_field_block(std::ostream& rOs, const std::string& rName, DType dt,
                       std::size_t num_components, std::size_t num_tuples, bool binary,
                       const std::vector<const NDArray*>& rBlocks) {
    if (rName.find(' ') != std::string::npos)
        throw WriteError("VTK doesn't support spaces in field names ('" + rName + "').");
    rOs << rName << ' ' << num_components << ' ' << num_tuples << ' ' << vtk_dtype_str(dt) << '\n';
    const bool flt = is_float_dtype(dt);
    for (const NDArray* blk : rBlocks) {
        if (binary) {
            std::vector<unsigned char> buf = be_buffer(*blk);
            rOs.write(reinterpret_cast<const char*>(buf.data()),
                      static_cast<std::streamsize>(buf.size()));
        } else {
            const std::size_t n = blk->Size();
            for (std::size_t i = 0; i < n; ++i) {
                if (flt)
                    vtk_ascii_double(rOs, read_double(*blk, i));
                else
                    rOs << read_int(*blk, i);
                rOs << ' ';
            }
        }
    }
    rOs << '\n';
}

}  // namespace

void write_vtk(const std::string& rPath, const Mesh& rMesh, bool binary, bool v51) {
    for (const auto cb : rMesh.CellRange())
        if (cb.Type().rfind("polyhedron", 0) == 0)
            throw WriteError("C++ VTK writer does not support polyhedron cells");

    std::ofstream os(rPath, std::ios::binary);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);

    const NDArray& points = rMesh.Points();
    const std::size_t num_points = rMesh.NumPoints();
    const std::size_t dim = rMesh.PointDim();
    const std::size_t pt_isz = dtype_size(points.Dtype());

    std::size_t total_cells = 0, total_idx = 0;
    for (const auto cb : rMesh.CellRange()) {
        total_cells += cb.NumCells();
        total_idx += cb.Conn().Size();
    }

    os << (v51 ? "# vtk DataFile Version 5.1\n" : "# vtk DataFile Version 4.2\n");
    os << "written by meshio++ (C++ core)\n";
    os << (binary ? "BINARY\n" : "ASCII\n");
    os << "DATASET UNSTRUCTURED_GRID\n";

    // Points (3 components; pad 2D with zero z).
    os << "POINTS " << num_points << ' ' << vtk_dtype_str(points.Dtype()) << '\n';
    if (binary) {
        // Pre-sized padded buffer, parallel byte-swap, then one write.
        const auto* src = reinterpret_cast<const char*>(points.Data());
        std::vector<unsigned char> buf(num_points * 3 * pt_isz, 0);
        auto* dst = reinterpret_cast<char*>(buf.data());
        const int isz = static_cast<int>(pt_isz);
        parallel_for_bw(num_points, [&](std::size_t r) {
            for (std::size_t c = 0; c < dim && c < 3; ++c)
                detail::bswap_copy(dst + (r * 3 + c) * pt_isz, src + (r * dim + c) * pt_isz, isz);
        });
        os.write(reinterpret_cast<const char*>(buf.data()),
                 static_cast<std::streamsize>(buf.size()));
        os << '\n';
    } else {
        for (std::size_t r = 0; r < num_points; ++r)
            for (std::size_t c = 0; c < 3; ++c) {
                vtk_ascii_double(os, (c < dim) ? read_double(points, r * dim + c) : 0.0);
                os << ((r + 1 == num_points && c == 2) ? '\n' : ' ');
            }
        if (num_points == 0)
            os << '\n';
    }

    if (v51) {
        // Version 5.1: OFFSETS (num_cells + 1) and CONNECTIVITY (total_idx).
        os << "CELLS " << (total_cells + 1) << ' ' << total_idx << '\n';
        os << "OFFSETS vtktypeint64\n";
        // Cumulative offsets (sequential prefix sum, cheap).
        std::vector<std::int64_t> offs(total_cells + 1);
        offs[0] = 0;
        std::size_t oi = 1;
        std::int64_t running = 0;
        for (const auto cb : rMesh.CellRange()) {
            const std::int64_t k = static_cast<std::int64_t>(cols(cb.Conn()));
            for (std::size_t r = 0; r < cb.NumCells(); ++r)
                offs[oi++] = (running += k);
        }
        if (binary) {
            write_be(os, offs);
            os << '\n';
            os << "CONNECTIVITY vtktypeint64\n";
            // Fused gather + big-endian store, one pass into the byte buffer.
            std::vector<unsigned char> cbuf(total_idx * 8);
            std::size_t base = 0;
            for (const auto cb : rMesh.CellRange()) {
                const NDArray& conn = cb.Conn();
                const std::size_t nc = cb.NumCells();
                const std::size_t k = cols(conn);
                std::vector<int> order = meshio_to_vtk_order(cb.Type());
                const int* ord = order.empty() ? nullptr : order.data();
                gather_be(cbuf.data() + base * 8, conn, nc, k, ord, 8);
                base += nc * k;
            }
            os.write(reinterpret_cast<const char*>(cbuf.data()),
                     static_cast<std::streamsize>(cbuf.size()));
            os << '\n';
        } else {
            for (std::int64_t v : offs)
                os << v << '\n';
            os << "CONNECTIVITY vtktypeint64\n";
            for (const auto cb : rMesh.CellRange()) {
                const NDArray& conn = cb.Conn();
                const std::size_t nc = cb.NumCells();
                const std::size_t k = cols(conn);
                std::vector<int> order = meshio_to_vtk_order(cb.Type());
                for (std::size_t r = 0; r < nc; ++r)
                    for (std::size_t j = 0; j < k; ++j) {
                        std::size_t col = order.empty() ? j : static_cast<std::size_t>(order[j]);
                        os << read_int(conn, r * k + col) << '\n';
                    }
            }
        }
    } else {
        // Version 4.2: interleaved [count, nodes...] per cell, as int32.
        os << "CELLS " << total_cells << ' ' << (total_idx + total_cells) << '\n';
        if (binary) {
            // Fused: each cell -> [k, v0..v(k-1)] big-endian int32, one pass.
            std::vector<unsigned char> cbuf((total_idx + total_cells) * 4);
            std::size_t p = 0;  // element index into cbuf
            for (const auto cb : rMesh.CellRange()) {
                const NDArray& conn = cb.Conn();
                const std::size_t nc = cb.NumCells();
                const std::size_t k = cols(conn);
                const std::size_t stride = k + 1;
                std::vector<int> order = meshio_to_vtk_order(cb.Type());
                const int* ord = order.empty() ? nullptr : order.data();
                const std::size_t block_base = p;
                dispatch_dtype(conn.Dtype(), [&]<class T>() {
                    const T* src = conn.As<T>();
                    parallel_for_bw(nc, [&](std::size_t r) {
                        unsigned char* o = cbuf.data() + (block_base + r * stride) * 4;
                        be_store(o, static_cast<std::uint64_t>(k), 4);
                        for (std::size_t j = 0; j < k; ++j) {
                            std::size_t col = ord ? static_cast<std::size_t>(ord[j]) : j;
                            auto v = static_cast<std::uint64_t>(
                                static_cast<std::int64_t>(src[r * k + col]));
                            be_store(o + (j + 1) * 4, v, 4);
                        }
                    });
                });
                p += nc * stride;
            }
            os.write(reinterpret_cast<const char*>(cbuf.data()),
                     static_cast<std::streamsize>(cbuf.size()));
            os << '\n';
        } else {
            for (const auto cb : rMesh.CellRange()) {
                const NDArray& conn = cb.Conn();
                const std::size_t nc = cb.NumCells();
                const std::size_t k = cols(conn);
                std::vector<int> order = meshio_to_vtk_order(cb.Type());
                for (std::size_t r = 0; r < nc; ++r) {
                    os << k << '\n';
                    for (std::size_t j = 0; j < k; ++j) {
                        std::size_t col = order.empty() ? j : static_cast<std::size_t>(order[j]);
                        os << read_int(conn, r * k + col) << '\n';
                    }
                }
            }
        }
    }

    // Cell types.
    os << "CELL_TYPES " << total_cells << '\n';
    const auto& tmap = meshio_to_vtk_type();
    std::vector<std::int32_t> ctypes(total_cells);
    std::size_t ci = 0;
    for (const auto cb : rMesh.CellRange()) {
        auto it = tmap.find(cb.Type());
        if (it == tmap.end())
            throw WriteError("Unknown cell type for VTK: " + cb.Type());
        for (std::size_t r = 0; r < cb.NumCells(); ++r)
            ctypes[ci++] = it->second;
    }
    if (binary) {
        write_be(os, ctypes);
        os << '\n';
    } else {
        for (std::int32_t v : ctypes)
            os << v << '\n';
    }

    // Point data.
    if (rMesh.NumPointData() != 0) {
        os << "POINT_DATA " << num_points << '\n';
        os << "FIELD FieldData " << rMesh.NumPointData() << '\n';
        for (const auto& name : rMesh.PointDataNames()) {
            const NDArray& d = rMesh.PointData(name);
            std::size_t ncomp = cols(d);
            write_field_block(os, name, d.Dtype(), ncomp, d.Shape().empty() ? 0 : d.Shape()[0],
                              binary, {&d});
        }
    }

    // Cell data (concatenate per-block arrays for each name).
    if (rMesh.NumCellData() != 0) {
        os << "CELL_DATA " << total_cells << '\n';
        os << "FIELD FieldData " << rMesh.NumCellData() << '\n';
        for (const auto& name : rMesh.CellDataNames()) {
            const std::size_t nblocks = rMesh.CellDataNumBlocks(name);
            if (nblocks == 0)
                continue;
            std::vector<const NDArray*> ptrs;
            for (std::size_t bi = 0; bi < nblocks; ++bi)
                ptrs.push_back(&rMesh.CellData(name, bi));
            const NDArray& first = *ptrs.front();
            write_field_block(os, name, first.Dtype(), cols(first), total_cells, binary, ptrs);
        }
    }
}

}  // namespace meshioplusplus
