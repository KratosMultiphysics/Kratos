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
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/ply.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/skin.hpp"

namespace meshioplusplus {

namespace {

DType ply_to_dtype(const std::string& rS) {
    if (rS == "char" || rS == "int8")
        return DType::Int8;
    if (rS == "uchar" || rS == "uint8")
        return DType::UInt8;
    if (rS == "short" || rS == "int16")
        return DType::Int16;
    if (rS == "ushort" || rS == "uint16")
        return DType::UInt16;
    if (rS == "int" || rS == "int32")
        return DType::Int32;
    if (rS == "uint" || rS == "uint32")
        return DType::UInt32;
    if (rS == "int64")
        return DType::Int64;
    if (rS == "uint64")
        return DType::UInt64;
    if (rS == "float" || rS == "float32")
        return DType::Float32;
    if (rS == "double" || rS == "float64")
        return DType::Float64;
    throw ReadError("PLY: unknown property type '" + rS + "'");
}

const char* dtype_to_ply(DType dt) {
    switch (dt) {
        case DType::Int8:
            return "int8";
        case DType::Int16:
            return "int16";
        case DType::Int32:
            return "int32";
        case DType::Int64:
            return "int64";
        case DType::UInt8:
            return "uint8";
        case DType::UInt16:
            return "uint16";
        case DType::UInt32:
            return "uint32";
        case DType::UInt64:
            return "uint64";
        case DType::Float32:
            return "float";
        case DType::Float64:
            return "double";
    }
    return "double";
}

std::string cell_type_from_count(std::size_t n) {
    switch (n) {
        case 1:
            return "vertex";
        case 2:
            return "line";
        case 3:
            return "triangle";
        case 4:
            return "quad";
        default:
            return "polygon";
    }
}

std::string ply_trim(const std::string& rS) {
    std::size_t b = 0, e = rS.size();
    while (b < e && std::isspace(static_cast<unsigned char>(rS[b])))
        ++b;
    while (e > b && std::isspace(static_cast<unsigned char>(rS[e - 1])))
        --e;
    return rS.substr(b, e - b);
}

struct VProp {
    std::string mName;
    DType mDtype;
};

// Read one scalar of `dt` from the buffer at pos into NDArray element idx,
// byte-swapping when the file is big-endian.
void rd_into(NDArray& rA, std::size_t idx, const std::string& rBuf, std::size_t& rPos, bool big) {
    std::size_t isz = dtype_size(rA.Dtype());
    unsigned char* dst = reinterpret_cast<unsigned char*>(rA.Data()) + idx * isz;
    if (rPos + isz > rBuf.size())
        throw ReadError("PLY binary truncated");
    if (big)
        for (std::size_t b = 0; b < isz; ++b)
            dst[b] = static_cast<unsigned char>(rBuf[rPos + isz - 1 - b]);
    else
        std::memcpy(dst, rBuf.data() + rPos, isz);
    rPos += isz;
}

std::int64_t rd_int_val(const std::string& rBuf, std::size_t& rPos, DType dt, bool big) {
    NDArray t(dt, {1});
    rd_into(t, 0, rBuf, rPos, big);
    return detail::read_int(t, 0);
}

void store_scalar(NDArray& rA, std::size_t idx, double dval, std::int64_t ival, bool isflt) {
    switch (rA.Dtype()) {
        case DType::Float32:
            rA.As<float>()[idx] = static_cast<float>(dval);
            break;
        case DType::Float64:
            rA.As<double>()[idx] = dval;
            break;
        case DType::Int8:
            rA.As<std::int8_t>()[idx] = static_cast<std::int8_t>(ival);
            break;
        case DType::Int16:
            rA.As<std::int16_t>()[idx] = static_cast<std::int16_t>(ival);
            break;
        case DType::Int32:
            rA.As<std::int32_t>()[idx] = static_cast<std::int32_t>(ival);
            break;
        case DType::Int64:
            rA.As<std::int64_t>()[idx] = ival;
            break;
        case DType::UInt8:
            rA.As<std::uint8_t>()[idx] = static_cast<std::uint8_t>(ival);
            break;
        case DType::UInt16:
            rA.As<std::uint16_t>()[idx] = static_cast<std::uint16_t>(ival);
            break;
        case DType::UInt32:
            rA.As<std::uint32_t>()[idx] = static_cast<std::uint32_t>(ival);
            break;
        case DType::UInt64:
            rA.As<std::uint64_t>()[idx] = static_cast<std::uint64_t>(ival);
            break;
    }
    (void)isflt;
}

}  // namespace

Mesh read_ply(const std::string& rPath) {
    std::ifstream in(rPath, std::ios::binary);
    if (!in)
        throw ReadError("Could not open file: " + rPath);
    std::string buf((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    std::size_t pos = 0;

    auto read_line = [&]() -> std::string {
        std::size_t start = pos;
        while (pos < buf.size() && buf[pos] != '\n')
            ++pos;
        std::string line = buf.substr(start, pos - start);
        if (pos < buf.size())
            ++pos;
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        return line;
    };
    auto next_sig = [&]() -> std::string {
        while (true) {
            std::string l = ply_trim(read_line());
            if (!l.empty() && l.rfind("comment", 0) != 0)
                return l;
        }
    };

    if (ply_trim(read_line()) != "ply")
        throw ReadError("Expected 'ply'");
    std::string fmt = next_sig();
    bool is_binary, big = false;
    if (fmt == "format ascii 1.0")
        is_binary = false;
    else if (fmt == "format binary_big_endian 1.0") {
        is_binary = true;
        big = true;
    } else if (fmt == "format binary_little_endian 1.0")
        is_binary = true;
    else
        throw ReadError("PLY: unknown format line");

    std::size_t num_verts = 0, num_faces = 0;
    std::vector<VProp> vprops;
    DType face_count_dt = DType::UInt8, face_index_dt = DType::Int32;
    bool have_face = false;

    std::string line = next_sig();
    while (line != "end_header") {
        std::istringstream iss(line);
        std::string tok;
        iss >> tok;
        if (tok == "obj_info") {
            line = next_sig();
        } else if (tok == "element") {
            std::string ename;
            std::size_t count;
            iss >> ename >> count;
            if (ename == "vertex") {
                num_verts = count;
                line = next_sig();
                while (line.rfind("property", 0) == 0) {
                    std::istringstream ps(line);
                    std::string p, type, name;
                    ps >> p >> type >> name;
                    if (type == "list")
                        throw ReadError("PLY: list vertex property not supported by C++");
                    vprops.push_back({name, ply_to_dtype(type)});
                    line = next_sig();
                }
            } else if (ename == "face") {
                num_faces = count;
                have_face = true;
                line = next_sig();
                bool got_list = false;
                while (line.rfind("property", 0) == 0) {
                    std::istringstream ps(line);
                    std::string p, kind;
                    ps >> p >> kind;
                    if (kind == "list") {
                        std::string ct, it_, nm;
                        ps >> ct >> it_ >> nm;
                        face_count_dt = ply_to_dtype(ct);
                        face_index_dt = ply_to_dtype(it_);
                        got_list = true;
                    } else {
                        throw ReadError("PLY: extra face properties not supported by C++");
                    }
                    line = next_sig();
                }
                if (!got_list && num_faces > 0)
                    throw ReadError("PLY: face element without vertex index list");
            } else {
                throw ReadError("PLY: unsupported element '" + ename + "'");
            }
        } else {
            throw ReadError("PLY: unexpected header line '" + line + "'");
        }
    }

    // Vertex properties -> per-property arrays.
    std::vector<NDArray> vcols;
    for (const auto& vp : vprops)
        vcols.emplace_back(vp.mDtype, std::vector<std::size_t>{num_verts});

    if (is_binary) {
        // Fixed-width records: property c of vertex i sits at a closed-form
        // byte offset -> decode + byteswap in parallel over vertices.
        std::size_t stride = 0;
        std::vector<std::size_t> coff(vprops.size());
        for (std::size_t c = 0; c < vprops.size(); ++c) {
            coff[c] = stride;
            stride += dtype_size(vprops[c].mDtype);
        }
        if (pos + num_verts * stride > buf.size())
            throw ReadError("PLY binary truncated");
        const std::size_t start = pos;
        parallel_for(num_verts, [&](std::size_t i) {
            for (std::size_t c = 0; c < vprops.size(); ++c) {
                const std::size_t isz = dtype_size(vcols[c].Dtype());
                unsigned char* dst = reinterpret_cast<unsigned char*>(vcols[c].Data()) + i * isz;
                const std::size_t src = start + i * stride + coff[c];
                if (big)
                    for (std::size_t b = 0; b < isz; ++b)
                        dst[b] = static_cast<unsigned char>(buf[src + isz - 1 - b]);
                else
                    std::memcpy(dst, buf.data() + src, isz);
            }
        });
        pos = start + num_verts * stride;
    } else {
        for (std::size_t i = 0; i < num_verts; ++i) {
            std::string row = read_line();
            std::istringstream rs(row);
            for (std::size_t c = 0; c < vprops.size(); ++c) {
                std::string t;
                rs >> t;
                if (detail::is_float_dtype(vcols[c].Dtype()))
                    store_scalar(vcols[c], i, std::strtod(t.c_str(), nullptr), 0, true);
                else
                    store_scalar(vcols[c], i, 0.0, std::strtoll(t.c_str(), nullptr, 10), false);
            }
        }
    }

    Mesh mesh;
    // Assemble points from x/y/z; the rest become point_data.
    std::vector<std::size_t> xyz(3, SIZE_MAX);
    for (std::size_t c = 0; c < vprops.size(); ++c) {
        if (vprops[c].mName == "x")
            xyz[0] = c;
        else if (vprops[c].mName == "y")
            xyz[1] = c;
        else if (vprops[c].mName == "z")
            xyz[2] = c;
    }
    std::size_t ndim = 0;
    for (std::size_t k = 0; k < 3; ++k)
        if (xyz[k] != SIZE_MAX)
            ++ndim;
    DType pdt = (xyz[0] != SIZE_MAX) ? vcols[xyz[0]].Dtype() : DType::Float64;
    NDArray pts(pdt, {num_verts, ndim});
    for (std::size_t i = 0; i < num_verts; ++i)
        for (std::size_t k = 0; k < ndim; ++k)
            store_scalar(pts, i * ndim + k, detail::read_double(vcols[xyz[k]], i),
                         detail::read_int(vcols[xyz[k]], i), detail::is_float_dtype(pdt));
    mesh.AssignPoints(std::move(pts));
    for (std::size_t c = 0; c < vprops.size(); ++c) {
        const std::string& nm = vprops[c].mName;
        if (nm == "x" || nm == "y" || nm == "z")
            continue;
        mesh.AddPointData(nm, std::move(vcols[c]));
    }

    // Faces -> cell blocks grouped by consecutive vertex count.
    if (have_face) {
        std::size_t cur_n = SIZE_MAX;
        std::vector<std::int64_t> cur_conn;
        std::size_t cur_count = 0;
        auto flush = [&]() {
            if (cur_count == 0)
                return;
            NDArray data(DType::Int64, {cur_count, cur_n});
            std::memcpy(data.Data(), cur_conn.data(), cur_conn.size() * sizeof(std::int64_t));
            mesh.AddCellBlock(cell_type_from_count(cur_n), std::move(data));
            cur_conn.clear();
            cur_count = 0;
        };
        for (std::size_t f = 0; f < num_faces; ++f) {
            std::size_t n;
            std::vector<std::int64_t> idx;
            if (is_binary) {
                n = static_cast<std::size_t>(rd_int_val(buf, pos, face_count_dt, big));
                idx.resize(n);
                for (std::size_t j = 0; j < n; ++j)
                    idx[j] = rd_int_val(buf, pos, face_index_dt, big);
            } else {
                std::istringstream rs(read_line());
                long long cnt;
                rs >> cnt;
                n = static_cast<std::size_t>(cnt);
                idx.resize(n);
                for (std::size_t j = 0; j < n; ++j)
                    rs >> idx[j];
            }
            if (n != cur_n) {
                flush();
                cur_n = n;
            }
            cur_conn.insert(cur_conn.end(), idx.begin(), idx.end());
            ++cur_count;
        }
        flush();
    }

    return mesh;
}

void write_ply(const std::string& rPath, const Mesh& rMesh, bool binary, bool skin) {
    if (skin && has_skinnable_cells(rMesh)) {
        std::size_t dropped = 0;
        for (const auto cb : rMesh.CellRange())
            if (cell_type_dimension(cell_type_from_name(cb.Type())) != 3)
                ++dropped;
        if (dropped > 0)
            log::warn(
                "PLY: writing the extracted skin of the volume cells; {} pre-existing "
                "non-volume cell block(s) dropped (pass skin=false for the legacy behavior).",
                dropped);
        write_ply(rPath, extract_skin(rMesh, /*linearize=*/true), binary, /*skin=*/false);
        return;
    }

    std::ofstream os(rPath, std::ios::binary);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);

    const std::size_t num_points = rMesh.NumPoints();
    const NDArray& points = rMesh.Points();
    const std::size_t dim = points.Shape().size() >= 2 ? points.Shape()[1] : 0;
    const std::size_t ncoord = std::min<std::size_t>(dim, 3);

    // Scalar point data only (PLY can't store multidimensional vertex data here).
    std::vector<std::pair<std::string, const NDArray*>> pd;
    for (const auto& name : rMesh.PointDataNames()) {
        const NDArray& d = rMesh.PointData(name);
        if (d.Shape().size() <= 1)
            pd.emplace_back(name, &d);
    }

    const char* legal[] = {"vertex", "line", "triangle", "quad", "polygon"};
    auto is_legal = [&](const std::string& t) {
        for (auto* l : legal)
            if (t == l)
                return true;
        return false;
    };
    std::size_t num_cells = 0;
    for (const auto cb : rMesh.CellRange())
        if (is_legal(cb.Type()))
            num_cells += cb.NumCells();

    os << "ply\n";
    os << (binary ? "format binary_little_endian 1.0\n" : "format ascii 1.0\n");
    os << "comment Created by meshio++ (C++ core)\n";
    os << "element vertex " << num_points << "\n";
    const char* dim_names[3] = {"x", "y", "z"};
    for (std::size_t k = 0; k < ncoord; ++k)
        os << "property " << dtype_to_ply(points.Dtype()) << " " << dim_names[k] << "\n";
    for (auto& p : pd)
        os << "property " << dtype_to_ply(p.second->Dtype()) << " " << p.first << "\n";
    if (num_cells > 0) {
        os << "element face " << num_cells << "\n";
        os << "property list uint8 int32 vertex_indices\n";
    }
    os << "end_header\n";

    const std::size_t pisz = dtype_size(points.Dtype());
    if (binary) {
        // Interleaved vertex records: coords then scalar point data.
        for (std::size_t i = 0; i < num_points; ++i) {
            for (std::size_t k = 0; k < ncoord; ++k)
                os.write(reinterpret_cast<const char*>(points.Data()) + (i * dim + k) * pisz, pisz);
            for (auto& p : pd) {
                std::size_t isz = dtype_size(p.second->Dtype());
                os.write(reinterpret_cast<const char*>(p.second->Data()) + i * isz, isz);
            }
        }
        for (const auto cb : rMesh.CellRange()) {
            if (!is_legal(cb.Type()))
                continue;
            const NDArray& conn = cb.Conn();
            std::size_t n = conn.Shape().size() >= 2 ? conn.Shape()[1] : 1;
            for (std::size_t r = 0; r < cb.NumCells(); ++r) {
                std::uint8_t cnt = static_cast<std::uint8_t>(n);
                os.write(reinterpret_cast<const char*>(&cnt), 1);
                for (std::size_t j = 0; j < n; ++j) {
                    std::int32_t v = static_cast<std::int32_t>(detail::read_int(conn, r * n + j));
                    os.write(reinterpret_cast<const char*>(&v), 4);
                }
            }
        }
    } else {
        char buf[40];
        for (std::size_t i = 0; i < num_points; ++i) {
            std::string row;
            for (std::size_t k = 0; k < ncoord; ++k) {
                if (k)
                    row += " ";
                std::snprintf(buf, sizeof(buf), "%.17g", detail::read_double(points, i * dim + k));
                row += buf;
            }
            for (auto& p : pd) {
                row += " ";
                if (detail::is_float_dtype(p.second->Dtype())) {
                    std::snprintf(buf, sizeof(buf), "%.17g", detail::read_double(*p.second, i));
                    row += buf;
                } else {
                    row += std::to_string(detail::read_int(*p.second, i));
                }
            }
            os << row << "\n";
        }
        for (const auto cb : rMesh.CellRange()) {
            if (!is_legal(cb.Type()))
                continue;
            const NDArray& conn = cb.Conn();
            std::size_t n = conn.Shape().size() >= 2 ? conn.Shape()[1] : 1;
            for (std::size_t r = 0; r < cb.NumCells(); ++r) {
                os << n;
                for (std::size_t j = 0; j < n; ++j)
                    os << " " << detail::read_int(conn, r * n + j);
                os << "\n";
            }
        }
    }
}

}  // namespace meshioplusplus
