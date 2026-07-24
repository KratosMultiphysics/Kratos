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
#include "meshioplusplus/detail/vtk_xml.hpp"
#include "meshioplusplus/detail/vtu_binary.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/vtu.hpp"

namespace meshioplusplus {

namespace {

using detail::vtu_to_int64;

// compression: 0 = none, 1 = zlib. lzma/appended raise (handled by caller).
NDArray vtu_read_data_array(const pugi::xml_node& rDa, int compression, std::size_t hsz,
                            int& rNumComponents) {
    std::string fmt = rDa.attribute("format").as_string("ascii");
    DType dt = detail::dtype_from_vtu(rDa.attribute("type").as_string());
    rNumComponents = rDa.attribute("NumberOfComponents").as_int(0);

    if (fmt == "ascii")
        return detail::vtu_parse_ascii(rDa.text().get(), dt);
    if (fmt == "binary")
        return detail::vtu_parse_binary(detail::vtu_strip(rDa.text().get()), dt, compression, hsz);
    throw ReadError("VTU '" + fmt + "' data is not supported by the C++ reader");
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
            NDArray pts = vtu_read_data_array(da, compression, hsz, nc);
            if (nc <= 0)
                nc = 3;
            pts.Reshape({num_points, static_cast<std::size_t>(nc)});
            mesh.AssignPoints(std::move(pts));
        } else if (tag == "Cells") {
            for (pugi::xml_node da : child.children("DataArray")) {
                int nc = 0;
                std::string name = da.attribute("Name").as_string();
                NDArray arr = vtu_read_data_array(da, compression, hsz, nc);
                if (name == "connectivity")
                    conn = vtu_to_int64(arr);
                else if (name == "offsets")
                    offsets = vtu_to_int64(arr);
                else if (name == "types")
                    types = vtu_to_int64(arr);
                else if (name == "faces" || name == "faceoffsets")
                    throw ReadError("polyhedron VTU not supported by the C++ reader");
            }
        } else if (tag == "PointData") {
            for (pugi::xml_node da : child.children("DataArray")) {
                int nc = 0;
                std::string name = da.attribute("Name").as_string();
                NDArray arr = vtu_read_data_array(da, compression, hsz, nc);
                if (nc > 1)
                    arr.Reshape({arr.Size() / nc, static_cast<std::size_t>(nc)});
                mesh.AddPointData(name, std::move(arr));
            }
        } else if (tag == "CellData") {
            for (pugi::xml_node da : child.children("DataArray")) {
                int nc = 0;
                std::string name = da.attribute("Name").as_string();
                NDArray arr = vtu_read_data_array(da, compression, hsz, nc);
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
