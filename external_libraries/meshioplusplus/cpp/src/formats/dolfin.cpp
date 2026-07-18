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
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

// External includes
#include "pugixml.hpp"

// Project includes
#include "meshioplusplus/formats/dolfin.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace fs = std::filesystem;

namespace meshioplusplus {

namespace {

std::pair<std::string, int> dolfin_to_meshio(const std::string& rCt) {
    if (rCt == "triangle")
        return {"triangle", 3};
    if (rCt == "tetrahedron")
        return {"tetra", 4};
    throw ReadError("DOLFIN: unsupported cell type '" + rCt + "'");
}

const char* meshio_to_dolfin(const std::string& rT) {
    if (rT == "triangle")
        return "triangle";
    if (rT == "tetra")
        return "tetrahedron";
    throw WriteError("DOLFIN XML only supports triangles and tetrahedra");
}

}  // namespace

Mesh read_dolfin(const std::string& rPath) {
    pugi::xml_document doc;
    if (!doc.load_file(rPath.c_str()))
        throw ReadError("DOLFIN: could not parse " + rPath);

    pugi::xml_node dolfin = doc.child("dolfin");
    if (!dolfin)
        throw ReadError("DOLFIN: missing <dolfin> root");
    pugi::xml_node mesh_node = dolfin.child("mesh");
    if (!mesh_node)
        throw ReadError("DOLFIN: missing <mesh>");

    int dim = mesh_node.attribute("dim").as_int();
    auto [cell_type, npc] = dolfin_to_meshio(mesh_node.attribute("celltype").value());

    Mesh mesh;

    // Vertices (placed by index).
    pugi::xml_node verts = mesh_node.child("vertices");
    std::size_t nverts = verts.attribute("size").as_uint();
    NDArray pts(DType::Float64, {nverts, static_cast<std::size_t>(dim)});
    double* pp = pts.As<double>();
    const char* coord[3] = {"x", "y", "z"};
    for (pugi::xml_node v : verts.children("vertex")) {
        std::size_t k = v.attribute("index").as_uint();
        for (int c = 0; c < dim; ++c)
            pp[k * dim + c] = v.attribute(coord[c]).as_double();
    }
    mesh.AssignPoints(std::move(pts));

    // Cells (single block, placed by index).
    pugi::xml_node cells = mesh_node.child("cells");
    std::size_t ncells = cells.attribute("size").as_uint();
    NDArray data(DType::Int64, {ncells, static_cast<std::size_t>(npc)});
    std::int64_t* dp = data.As<std::int64_t>();
    for (pugi::xml_node c : cells.children()) {
        std::size_t k = c.attribute("index").as_uint();
        for (int j = 0; j < npc; ++j) {
            char tag[16];  // "v" + up to 11 digits (INT_MIN) + '\0'; GCC's static
                           // format-truncation analysis cannot prove j is small
            std::snprintf(tag, sizeof(tag), "v%d", j);
            dp[k * npc + j] = c.attribute(tag).as_llong();
        }
    }
    mesh.AddCellBlock(cell_type, std::move(data));

    // Cell data: sibling files "<stem>_<name>.xml".
    fs::path p(rPath);
    fs::path dir = p.has_parent_path() ? p.parent_path() : fs::path(".");
    std::string stem = p.stem().string();
    std::string prefix = stem + "_";
    if (fs::exists(dir)) {
        for (const auto& entry : fs::directory_iterator(dir)) {
            std::string fname = entry.path().filename().string();
            if (fname.size() <= prefix.size() + 4)
                continue;
            if (fname.compare(0, prefix.size(), prefix) != 0)
                continue;
            if (fname.compare(fname.size() - 4, 4, ".xml") != 0)
                continue;
            std::string name = fname.substr(prefix.size(), fname.size() - prefix.size() - 4);
            if (name.empty() || name.find('.') != std::string::npos)
                continue;  // [^.]+

            pugi::xml_document fdoc;
            if (!fdoc.load_file(entry.path().string().c_str()))
                continue;
            pugi::xml_node mf = fdoc.child("dolfin").child("mesh_function");
            if (!mf)
                continue;
            std::string type = mf.attribute("type").value();
            std::size_t size = mf.attribute("size").as_uint();
            DType dt = (type == "float") ? DType::Float64 : DType::Int64;
            NDArray arr(dt, {size});
            for (pugi::xml_node e : mf.children("entity")) {
                std::size_t idx = e.attribute("index").as_uint();
                if (dt == DType::Float64)
                    arr.As<double>()[idx] = e.attribute("value").as_double();
                else
                    arr.As<std::int64_t>()[idx] = e.attribute("value").as_llong();
            }
            std::vector<NDArray> blocks;
            blocks.push_back(std::move(arr));
            mesh.AddCellData(name, std::move(blocks));
        }
    }

    return mesh;
}

void write_dolfin(const std::string& rPath, const Mesh& rMesh) {
    // Pick the single supported cell type to write.
    std::string cell_type;
    for (const auto cb : rMesh.CellRange())
        if (cb.Type() == "tetra") {
            cell_type = "tetra";
            break;
        }
    if (cell_type.empty())
        for (const auto cb : rMesh.CellRange())
            if (cb.Type() == "triangle") {
                cell_type = "triangle";
                break;
            }
    if (cell_type.empty())
        throw WriteError("DOLFIN XML only supports triangles and tetrahedra");

    const std::size_t dim = rMesh.PointDim();
    if (dim != 2 && dim != 3)
        throw WriteError("DOLFIN: can only write dimension 2 or 3");

    std::ofstream f(rPath, std::ios::binary);
    if (!f)
        throw WriteError("Could not open file for writing: " + rPath);

    f << "<dolfin nsmap=\"{'dolfin': 'https://fenicsproject.org/'}\">\n";
    f << "  <mesh celltype=\"" << meshio_to_dolfin(cell_type) << "\" dim=\"" << dim << "\">\n";

    const std::size_t npts = rMesh.NumPoints();
    const NDArray& points = rMesh.Points();
    f << "    <vertices size=\"" << npts << "\">\n";
    char buf[32];
    const char* coord[3] = {"x", "y", "z"};
    for (std::size_t i = 0; i < npts; ++i) {
        f << "      <vertex index=\"" << i << "\"";
        for (std::size_t c = 0; c < dim; ++c) {
            std::snprintf(buf, sizeof(buf), "%.17g", detail::read_double(points, i * dim + c));
            f << " " << coord[c] << "=\"" << buf << "\"";
        }
        f << " />\n";
    }
    f << "    </vertices>\n";

    std::size_t num_cells = 0;
    for (const auto cb : rMesh.CellRange())
        if (cb.Type() == cell_type)
            num_cells += cb.NumCells();

    f << "    <cells size=\"" << num_cells << "\">\n";
    const char* ts = meshio_to_dolfin(cell_type);
    std::size_t idx = 0;
    for (const auto cb : rMesh.CellRange()) {
        if (cb.Type() != cell_type)
            continue;
        const NDArray& conn = cb.Conn();
        std::size_t ncols = detail::cols(conn);
        std::size_t n = cb.NumCells();
        for (std::size_t r = 0; r < n; ++r) {
            f << "      <" << ts << " index=\"" << idx << "\"";
            for (std::size_t j = 0; j < ncols; ++j)
                f << " v" << j << "=\"" << detail::read_int(conn, r * ncols + j) << "\"";
            f << " />\n";
            ++idx;
        }
    }
    f << "    </cells>\n";
    f << "  </mesh>\n";
    f << "</dolfin>";

    // Cell data -> sibling files "<stem>_<name>.xml".
    bool z_all_zero = true;
    if (dim == 3) {
        for (std::size_t i = 0; i < npts; ++i)
            if (detail::read_double(points, i * 3 + 2) != 0.0) {
                z_all_zero = false;
                break;
            }
    }
    int data_dim = (dim == 2 || z_all_zero) ? 2 : 3;

    fs::path p(rPath);
    std::string base = (p.parent_path() / p.stem()).string();
    for (const auto& name : rMesh.CellDataNames()) {
        const std::size_t nblocks = rMesh.CellDataNumBlocks(name);
        for (std::size_t bi = 0; bi < nblocks; ++bi) {
            const NDArray& arr = rMesh.CellData(name, bi);
            std::string fn = base + "_" + name + ".xml";
            std::ofstream cf(fn, std::ios::binary);
            if (!cf)
                throw WriteError("Could not open file for writing: " + fn);
            bool is_float = detail::is_float_dtype(arr.Dtype());
            const char* type = is_float ? "float" : "int";
            std::size_t sz = arr.Shape().empty() ? 0 : arr.Shape()[0];
            cf << "<dolfin><mesh_function type=\"" << type << "\" dim=\"" << data_dim
               << "\" size=\"" << sz << "\">";
            for (std::size_t k = 0; k < sz; ++k) {
                cf << "<entity index=\"" << k << "\" value=\"";
                if (is_float) {
                    std::snprintf(buf, sizeof(buf), "%.17g", detail::read_double(arr, k));
                    cf << buf;
                } else {
                    cf << detail::read_int(arr, k);
                }
                cf << "\" />";
            }
            cf << "</mesh_function></dolfin>";
        }
    }
}

}  // namespace meshioplusplus
