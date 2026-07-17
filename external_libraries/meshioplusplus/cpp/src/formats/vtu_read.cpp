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
#include <cstdlib>
#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>

// External includes
#include "pugixml.hpp"

// Project includes
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/vtk_cells.hpp"
#include "meshioplusplus/detail/vtu_binary.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/vtu.hpp"
#include "meshioplusplus/types.hpp"
#include "meshioplusplus/vtk_common.hpp"

namespace meshioplusplus {

namespace {

DType dtype_from_vtu(const std::string& rS) {
    if (rS == "Float32")
        return DType::Float32;
    if (rS == "Float64")
        return DType::Float64;
    if (rS == "Int8")
        return DType::Int8;
    if (rS == "Int16")
        return DType::Int16;
    if (rS == "Int32")
        return DType::Int32;
    if (rS == "Int64")
        return DType::Int64;
    if (rS == "UInt8")
        return DType::UInt8;
    if (rS == "UInt16")
        return DType::UInt16;
    if (rS == "UInt32")
        return DType::UInt32;
    if (rS == "UInt64")
        return DType::UInt64;
    throw ReadError("Illegal VTU data type '" + rS + "'");
}

void store(NDArray& rA, std::size_t i, double d, std::int64_t v, bool isflt) {
    switch (rA.Dtype()) {
        case DType::Float32:
            rA.As<float>()[i] = static_cast<float>(d);
            break;
        case DType::Float64:
            rA.As<double>()[i] = d;
            break;
        case DType::Int8:
            rA.As<std::int8_t>()[i] = static_cast<std::int8_t>(v);
            break;
        case DType::Int16:
            rA.As<std::int16_t>()[i] = static_cast<std::int16_t>(v);
            break;
        case DType::Int32:
            rA.As<std::int32_t>()[i] = static_cast<std::int32_t>(v);
            break;
        case DType::Int64:
            rA.As<std::int64_t>()[i] = v;
            break;
        case DType::UInt8:
            rA.As<std::uint8_t>()[i] = static_cast<std::uint8_t>(v);
            break;
        case DType::UInt16:
            rA.As<std::uint16_t>()[i] = static_cast<std::uint16_t>(v);
            break;
        case DType::UInt32:
            rA.As<std::uint32_t>()[i] = static_cast<std::uint32_t>(v);
            break;
        case DType::UInt64:
            rA.As<std::uint64_t>()[i] = static_cast<std::uint64_t>(v);
            break;
    }
    (void)isflt;
}

NDArray parse_ascii(const char* pText, DType dt) {
    const bool isflt = detail::is_float_dtype(dt);
    std::vector<double> dv;
    std::vector<std::int64_t> iv;
    const char* p = pText ? pText : "";
    while (*p) {
        while (*p && std::isspace(static_cast<unsigned char>(*p)))
            ++p;
        if (!*p)
            break;
        char* endp = nullptr;
        if (isflt) {
            double x = std::strtod(p, &endp);
            if (endp == p)
                break;
            dv.push_back(x);
        } else {
            long long x = std::strtoll(p, &endp, 10);
            if (endp == p)
                break;
            iv.push_back(static_cast<std::int64_t>(x));
        }
        p = endp;
    }
    std::size_t n = isflt ? dv.size() : iv.size();
    NDArray a(dt, {n});
    for (std::size_t i = 0; i < n; ++i)
        store(a, i, isflt ? dv[i] : 0.0, isflt ? 0 : iv[i], isflt);
    return a;
}

std::string strip(const char* pS) {
    std::string t = pS ? pS : "";
    std::size_t b = 0, e = t.size();
    while (b < e && std::isspace(static_cast<unsigned char>(t[b])))
        ++b;
    while (e > b && std::isspace(static_cast<unsigned char>(t[e - 1])))
        --e;
    return t.substr(b, e - b);
}

NDArray parse_binary(const std::string& rText, DType dt, int compression, std::size_t hsz) {
    std::vector<unsigned char> bytes;
    if (compression == 0)
        bytes = detail::vtu_decode_uncompressed(rText.c_str(), rText.size(), hsz);
    else
        bytes = detail::vtu_decode_zlib(rText.c_str(), rText.size(), hsz);
    std::size_t isz = dtype_size(dt);
    std::size_t n = isz ? bytes.size() / isz : 0;
    NDArray a(dt, {n});
    if (n)
        std::memcpy(a.Data(), bytes.data(), n * isz);
    return a;
}

// compression: 0 = none, 1 = zlib. lzma/appended raise (handled by caller).
NDArray read_data_array(const pugi::xml_node& rDa, int compression, std::size_t hsz,
                        int& rNumComponents) {
    std::string fmt = rDa.attribute("format").as_string("ascii");
    DType dt = dtype_from_vtu(rDa.attribute("type").as_string());
    rNumComponents = rDa.attribute("NumberOfComponents").as_int(0);

    if (fmt == "ascii")
        return parse_ascii(rDa.text().get(), dt);
    if (fmt == "binary")
        return parse_binary(strip(rDa.text().get()), dt, compression, hsz);
    throw ReadError("VTU '" + fmt + "' data is not supported by the C++ reader");
}

std::vector<std::int64_t> to_int64(const NDArray& rA) {
    std::vector<std::int64_t> v(rA.Size());
    for (std::size_t i = 0; i < rA.Size(); ++i)
        v[i] = detail::read_int(rA, i);
    return v;
}

}  // namespace

Mesh read_vtu(const std::string& rPath) {
    pugi::xml_document doc;
    pugi::xml_parse_result res = doc.load_file(rPath.c_str());
    if (!res)
        throw ReadError(std::string("VTU XML parse failed: ") + res.description());

    pugi::xml_node root = doc.child("VTKFile");
    if (!root)
        throw ReadError("Expected tag 'VTKFile'");
    if (std::string(root.attribute("type").as_string()) != "UnstructuredGrid")
        throw ReadError("Expected type UnstructuredGrid");

    int compression = 0;  // 0 none, 1 zlib
    std::string compressor = root.attribute("compressor").as_string("");
    if (compressor == "vtkZLibDataCompressor")
        compression = 1;
    else if (compressor == "vtkLZMADataCompressor")
        throw ReadError("lzma-compressed VTU not supported by the C++ reader");
    else if (!compressor.empty())
        throw ReadError("Unknown VTU compressor '" + compressor + "'");

    std::string header_type = root.attribute("header_type").as_string("UInt32");
    std::size_t hsz = (header_type == "UInt64") ? 8 : 4;

    pugi::xml_node grid = root.child("UnstructuredGrid");
    if (!grid)
        throw ReadError("No UnstructuredGrid found");

    // Appended data is not handled here -> let the Python reader take over.
    if (grid.parent().child("AppendedData") || root.child("AppendedData"))
        throw ReadError("appended VTU data not supported by the C++ reader");

    pugi::xml_node piece = grid.child("Piece");
    if (!piece)
        throw ReadError("No Piece found");
    // A single piece is supported; multiple pieces -> Python reader.
    if (piece.next_sibling("Piece"))
        throw ReadError("multi-piece VTU not supported by the C++ reader");

    std::size_t num_points =
        static_cast<std::size_t>(piece.attribute("NumberOfPoints").as_ullong());

    Mesh mesh;
    std::vector<std::int64_t> conn, offsets, types;
    std::unordered_map<std::string, NDArray> cell_data_raw;

    for (pugi::xml_node child : piece.children()) {
        std::string tag = child.name();
        if (tag == "Points") {
            pugi::xml_node da = child.child("DataArray");
            int nc = 0;
            NDArray pts = read_data_array(da, compression, hsz, nc);
            if (nc <= 0)
                nc = 3;
            pts.Reshape({num_points, static_cast<std::size_t>(nc)});
            mesh.AssignPoints(std::move(pts));
        } else if (tag == "Cells") {
            for (pugi::xml_node da : child.children("DataArray")) {
                int nc = 0;
                std::string name = da.attribute("Name").as_string();
                NDArray arr = read_data_array(da, compression, hsz, nc);
                if (name == "connectivity")
                    conn = to_int64(arr);
                else if (name == "offsets")
                    offsets = to_int64(arr);
                else if (name == "types")
                    types = to_int64(arr);
                else if (name == "faces" || name == "faceoffsets")
                    throw ReadError("polyhedron VTU not supported by the C++ reader");
            }
        } else if (tag == "PointData") {
            for (pugi::xml_node da : child.children("DataArray")) {
                int nc = 0;
                std::string name = da.attribute("Name").as_string();
                NDArray arr = read_data_array(da, compression, hsz, nc);
                if (nc > 1)
                    arr.Reshape({arr.Size() / nc, static_cast<std::size_t>(nc)});
                mesh.AddPointData(name, std::move(arr));
            }
        } else if (tag == "CellData") {
            for (pugi::xml_node da : child.children("DataArray")) {
                int nc = 0;
                std::string name = da.attribute("Name").as_string();
                NDArray arr = read_data_array(da, compression, hsz, nc);
                if (nc > 1)
                    arr.Reshape({arr.Size() / nc, static_cast<std::size_t>(nc)});
                cell_data_raw.emplace(name, std::move(arr));
            }
        }
    }

    detail::reconstruct_cells(conn.data(), offsets, types, cell_data_raw, mesh);
    return mesh;
}

}  // namespace meshioplusplus
