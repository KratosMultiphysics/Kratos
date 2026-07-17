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
#ifdef MESHIOPLUSPLUS_HAS_HDF5

// System includes
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/hmf.hpp"
#include "meshioplusplus/detail/hdf5_util.hpp"
#include "meshioplusplus/detail/xdmf_common.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus {

Mesh read_hmf(const std::string& rPath) {
    h5::SilenceErrors silence;
    h5::Hid f = h5::open_file_read(rPath);

    if (h5::read_attr_string(f, "type") != "hmf")
        throw ReadError("HMF: not an hmf file");
    if (h5::read_attr_string(f, "version") != "0.1-alpha")
        throw ReadError("HMF: unsupported version");

    h5::Hid domain = h5::open_group(f, "domain");
    h5::Hid grid = h5::open_group(domain, "grid");

    Mesh mesh;
    // Mirrors the Python reader's dict semantics: one entry per meshio type,
    // a repeated type replaces the earlier data; insertion order preserved.
    std::vector<std::pair<std::string, NDArray>> cells;
    std::vector<std::pair<std::string, NDArray>> cell_data_raw;

    for (const std::string& key : h5::group_links(grid)) {
        if (key.rfind("Topology", 0) == 0) {
            h5::Hid d(H5Dopen2(grid, key.c_str(), H5P_DEFAULT), H5Dclose);
            if (!d.Valid())
                throw ReadError("HMF: could not open " + key);
            std::string xt = h5::read_attr_string(d, "TopologyType");
            std::string mt = xdmfcommon::xdmf_to_meshio(xt);
            NDArray data = h5::read_dataset(grid, key);
            bool replaced = false;
            for (auto& kv : cells)
                if (kv.first == mt) {
                    kv.second = std::move(data);
                    replaced = true;
                    break;
                }
            if (!replaced)
                cells.emplace_back(mt, std::move(data));
        } else if (key == "Geometry") {
            h5::Hid d(H5Dopen2(grid, key.c_str(), H5P_DEFAULT), H5Dclose);
            std::string gt = h5::read_attr_string(d, "GeometryType");
            if (gt != "X" && gt != "XY" && gt != "XYZ")
                throw ReadError("HMF: unexpected GeometryType " + gt);
            mesh.AssignPoints(h5::read_dataset(grid, key));
        } else if (key == "CellAttributes") {
            h5::Hid g = h5::open_group(grid, key);
            for (const std::string& name : h5::group_links(g))
                cell_data_raw.emplace_back(name, h5::read_dataset(g, name));
        } else if (key == "NodeAttributes") {
            h5::Hid g = h5::open_group(grid, key);
            for (const std::string& name : h5::group_links(g))
                mesh.AddPointData(name, h5::read_dataset(g, name));
        } else {
            throw ReadError("HMF: unexpected entry " + key);
        }
    }

    for (auto& kv : cells)
        mesh.AddCellBlock(std::move(kv.first), std::move(kv.second));

    std::vector<std::size_t> sizes;
    for (const auto cb : mesh.CellRange())
        sizes.push_back(cb.NumCells());
    for (auto& kv : cell_data_raw)
        mesh.AddCellData(kv.first, xdmfcommon::split_raw_cell_data(kv.second, sizes));

    return mesh;
}

void write_hmf(const std::string& rPath, const Mesh& rMesh, int gzip_level) {
    h5::SilenceErrors silence;
    h5::Hid f = h5::create_file(rPath);

    h5::write_attr_string(f, "type", "hmf");
    h5::write_attr_string(f, "version", "0.1-alpha");

    h5::Hid domain = h5::create_group(f, "domain");
    h5::Hid grid = h5::create_group(domain, "grid");

    // Geometry
    {
        h5::write_dataset(grid, "Geometry", rMesh.Points(), gzip_level);
        h5::Hid d(H5Dopen2(grid, "Geometry", H5P_DEFAULT), H5Dclose);
        const std::size_t dim = rMesh.PointDim();
        h5::write_attr_string(d, "GeometryType", std::string("XYZ").substr(0, dim));
    }

    // Topology{k}
    for (std::size_t k = 0; k < rMesh.NumCellBlocks(); ++k) {
        const auto cb = rMesh.Cells(k);
        std::string name = "Topology" + std::to_string(k);
        h5::write_dataset(grid, name, cb.Conn(), gzip_level);
        h5::Hid d(H5Dopen2(grid, name.c_str(), H5P_DEFAULT), H5Dclose);
        h5::write_attr_string(d, "TopologyType", xdmfcommon::meshio_to_xdmf(cb.Type()));
    }

    // NodeAttributes / CellAttributes (sorted key order for deterministic output)
    h5::Hid na = h5::create_group(grid, "NodeAttributes");
    for (const auto& name : rMesh.PointDataNames())
        h5::write_dataset(na, name, rMesh.PointData(name), gzip_level);

    h5::Hid ca = h5::create_group(grid, "CellAttributes");
    for (const auto& name : rMesh.CellDataNames()) {
        if (rMesh.CellDataNumBlocks(name) == 0)
            continue;
        h5::write_dataset(ca, name, xdmfcommon::concat_cell_data(rMesh, name), gzip_level);
    }
}

}  // namespace meshioplusplus

#endif  // MESHIOPLUSPLUS_HAS_HDF5
