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
#include <string>
#include <unordered_map>
#include <vector>

// External includes
#include "pugixml.hpp"

// Project includes
#include "meshioplusplus/detail/vtk_cells.hpp"
#include "meshioplusplus/detail/vtk_xml.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/vtp.hpp"

namespace meshioplusplus {

namespace {

using detail::vtu_to_int64;

// compression: 0 = none, 1 = zlib. lzma/appended raise (handled by caller).
NDArray vtp_read_data_array(const pugi::xml_node& rDa, int compression, std::size_t hsz,
                            int& rNumComponents) {
    std::string fmt = rDa.attribute("format").as_string("ascii");
    DType dt = detail::dtype_from_vtu(rDa.attribute("type").as_string());
    rNumComponents = rDa.attribute("NumberOfComponents").as_int(0);

    if (fmt == "ascii")
        return detail::vtu_parse_ascii(rDa.text().get(), dt);
    if (fmt == "binary")
        return detail::vtu_parse_binary(detail::vtu_strip(rDa.text().get()), dt, compression, hsz);
    throw ReadError("VTP '" + fmt + "' data is not supported by the C++ reader");
}

// One PolyData section's connectivity + VTK end-offsets.
struct VtpPiece {
    std::vector<std::int64_t> mConn;
    std::vector<std::int64_t> mOffsets;
    bool mPresent = false;
};

VtpPiece vtp_read_section(const pugi::xml_node& rSection, int compression, std::size_t hsz) {
    VtpPiece out;
    if (!rSection)
        return out;
    out.mPresent = true;
    for (pugi::xml_node da : rSection.children("DataArray")) {
        std::string name = da.attribute("Name").as_string();
        int nc = 0;
        NDArray arr = vtp_read_data_array(da, compression, hsz, nc);
        if (name == "connectivity")
            out.mConn = vtu_to_int64(arr);
        else if (name == "offsets")
            out.mOffsets = vtu_to_int64(arr);
    }
    return out;
}

}  // namespace

Mesh read_vtp(const std::string& rPath) {
    pugi::xml_document doc;
    pugi::xml_parse_result res = doc.load_file(rPath.c_str());
    if (!res)
        throw ReadError(std::string("VTP XML parse failed: ") + res.description());

    pugi::xml_node root = doc.child("VTKFile");
    if (!root)
        throw ReadError("Expected tag 'VTKFile'");
    if (std::string(root.attribute("type").as_string()) != "PolyData")
        throw ReadError("Expected type PolyData");

    int compression = 0;  // 0 none, 1 zlib
    std::string compressor = root.attribute("compressor").as_string("");
    if (compressor == "vtkZLibDataCompressor")
        compression = 1;
    else if (compressor == "vtkLZMADataCompressor")
        throw ReadError("lzma-compressed VTP not supported by the C++ reader");
    else if (!compressor.empty())
        throw ReadError("Unknown VTP compressor '" + compressor + "'");

    std::string header_type = root.attribute("header_type").as_string("UInt32");
    std::size_t hsz = (header_type == "UInt64") ? 8 : 4;

    pugi::xml_node grid = root.child("PolyData");
    if (!grid)
        throw ReadError("No PolyData found");

    // Appended data is not handled here -> let the Python reader take over.
    if (grid.parent().child("AppendedData") || root.child("AppendedData"))
        throw ReadError("appended VTP data not supported by the C++ reader");

    pugi::xml_node piece = grid.child("Piece");
    if (!piece)
        throw ReadError("No Piece found");
    // A single piece is supported; multiple pieces -> Python reader.
    if (piece.next_sibling("Piece"))
        throw ReadError("multi-piece VTP not supported by the C++ reader");

    std::size_t num_points =
        static_cast<std::size_t>(piece.attribute("NumberOfPoints").as_ullong());

    Mesh mesh;
    std::unordered_map<std::string, NDArray> cell_data_raw;

    for (pugi::xml_node child : piece.children()) {
        std::string tag = child.name();
        if (tag == "Points") {
            pugi::xml_node da = child.child("DataArray");
            int nc = 0;
            NDArray pts = vtp_read_data_array(da, compression, hsz, nc);
            if (nc <= 0)
                nc = 3;
            pts.Reshape({num_points, static_cast<std::size_t>(nc)});
            mesh.AssignPoints(std::move(pts));
        } else if (tag == "PointData") {
            for (pugi::xml_node da : child.children("DataArray")) {
                int nc = 0;
                std::string name = da.attribute("Name").as_string();
                NDArray arr = vtp_read_data_array(da, compression, hsz, nc);
                if (nc > 1)
                    arr.Reshape({arr.Size() / nc, static_cast<std::size_t>(nc)});
                mesh.AddPointData(name, std::move(arr));
            }
        } else if (tag == "CellData") {
            for (pugi::xml_node da : child.children("DataArray")) {
                int nc = 0;
                std::string name = da.attribute("Name").as_string();
                NDArray arr = vtp_read_data_array(da, compression, hsz, nc);
                if (nc > 1)
                    arr.Reshape({arr.Size() / nc, static_cast<std::size_t>(nc)});
                cell_data_raw.emplace(name, std::move(arr));
            }
        }
    }

    VtpPiece verts = vtp_read_section(piece.child("Verts"), compression, hsz);
    VtpPiece lines = vtp_read_section(piece.child("Lines"), compression, hsz);
    VtpPiece polys = vtp_read_section(piece.child("Polys"), compression, hsz);
    VtpPiece strips = vtp_read_section(piece.child("Strips"), compression, hsz);
    if (!strips.mOffsets.empty())
        throw ReadError("triangle-strip VTP cells not supported by the C++ reader");

    // Concatenate sections in VTK's canonical PolyData cell order (Verts,
    // Lines, Polys), synthesizing a VTK type id per row so the shared
    // reconstruction (detail/vtk_cells.hpp) can build the blocks and split
    // cell_data.
    std::vector<std::int64_t> conn, offsets, types;
    auto append_section = [&](const VtpPiece& rSec, int kind) {
        const std::int64_t conn_base = static_cast<std::int64_t>(conn.size());
        std::int64_t prev = 0;
        for (std::int64_t end : rSec.mOffsets) {
            const std::int64_t sz = end - prev;
            prev = end;
            std::int64_t vtk_type = 0;
            if (kind == 0) {
                if (sz != 1)
                    throw ReadError("poly-vertex VTP cells not supported by the C++ reader");
                vtk_type = 1;  // VTK_VERTEX
            } else if (kind == 1) {
                if (sz != 2)
                    throw ReadError("poly-line VTP cells not supported by the C++ reader");
                vtk_type = 3;  // VTK_LINE
            } else {
                vtk_type = sz == 3 ? 5 : sz == 4 ? 9 : 7;  // triangle / quad / polygon
            }
            types.push_back(vtk_type);
            offsets.push_back(conn_base + end);
        }
        conn.insert(conn.end(), rSec.mConn.begin(), rSec.mConn.end());
    };
    append_section(verts, 0);
    append_section(lines, 1);
    append_section(polys, 2);

    detail::reconstruct_cells(conn.data(), offsets, types, cell_data_raw, mesh);
    return mesh;
}

}  // namespace meshioplusplus
