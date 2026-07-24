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
#include <cstdint>
#include <cstring>
#include <optional>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/formats/cgns.hpp"
#include "meshioplusplus/detail/hdf5_util.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus {

Mesh read_cgns(const std::string& rPath) {
    h5::SilenceErrors silence;
    h5::Hid f = h5::open_file_read(rPath);

    if (!h5::exists(f, "Base"))
        throw ReadError("Expected \"Base\" in file. Malformed CGNS?");
    h5::Hid base = h5::open_group(f, "Base");
    if (!h5::exists(base, "Zone1"))
        throw ReadError("Expected \"Zone1\" in \"Base\". Malformed CGNS?");
    h5::Hid zone = h5::open_group(base, "Zone1");

    h5::Hid coords = h5::open_group(zone, "GridCoordinates");
    h5::Hid gx = h5::open_group(coords, "CoordinateX");
    h5::Hid gy = h5::open_group(coords, "CoordinateY");
    h5::Hid gz = h5::open_group(coords, "CoordinateZ");
    NDArray x = h5::read_dataset(gx, " data");
    NDArray y = h5::read_dataset(gy, " data");
    NDArray z = h5::read_dataset(gz, " data");

    const std::size_t n = x.Shape().empty() ? 0 : x.Shape()[0];
    Mesh mesh;
    NDArray pts(DType::Float64, {n, 3});
    double* pp = pts.As<double>();
    for (std::size_t i = 0; i < n; ++i) {
        pp[i * 3 + 0] = detail::read_double(x, i);
        pp[i * 3 + 1] = detail::read_double(y, i);
        pp[i * 3 + 2] = detail::read_double(z, i);
    }
    mesh.AssignPoints(std::move(pts));

    h5::Hid elems = h5::open_group(zone, "GridElements");
    h5::Hid rng = h5::open_group(elems, "ElementRange");
    h5::Hid conn = h5::open_group(elems, "ElementConnectivity");
    NDArray range = h5::read_dataset(rng, " data");
    NDArray flat = h5::read_dataset(conn, " data");

    if (range.Size() < 2)
        throw ReadError("CGNS: malformed ElementRange");
    std::int64_t idx_max = detail::read_int(range, 1);
    if (idx_max <= 0 || flat.Size() % static_cast<std::size_t>(idx_max) != 0)
        throw ReadError("CGNS: malformed ElementConnectivity");
    std::size_t k = flat.Size() / static_cast<std::size_t>(idx_max);
    if (k != 4)
        throw ReadError("Can only read tetrahedra.");

    NDArray cells(flat.Dtype(), {static_cast<std::size_t>(idx_max), k});
    // shift 1-based -> 0-based, preserving the stored integer dtype
    for (std::size_t i = 0; i < flat.Size(); ++i) {
        std::int64_t v = detail::read_int(flat, i) - 1;
        switch (cells.Dtype()) {
            case DType::Int32:
                cells.As<std::int32_t>()[i] = static_cast<std::int32_t>(v);
                break;
            case DType::Int64:
                cells.As<std::int64_t>()[i] = v;
                break;
            case DType::UInt32:
                cells.As<std::uint32_t>()[i] = static_cast<std::uint32_t>(v);
                break;
            case DType::UInt64:
                cells.As<std::uint64_t>()[i] = static_cast<std::uint64_t>(v);
                break;
            default:
                throw ReadError("CGNS: unexpected connectivity dtype");
        }
    }
    mesh.AddCellBlock("tetra", std::move(cells));
    return mesh;
}

void write_cgns(const std::string& rPath, const Mesh& rMesh, int gzip_level) {
    h5::SilenceErrors silence;

    // Locate the tetra block (mirroring the Python writer, which only emits tetra).
    std::optional<Mesh::CellView> tet;
    for (const auto cb : rMesh.CellRange())
        if (cb.Type() == "tetra") {
            tet = cb;
            break;
        }

    h5::Hid f = h5::create_file(rPath);
    h5::Hid base = h5::create_group(f, "Base");
    h5::Hid zone = h5::create_group(base, "Zone1");
    h5::Hid coords = h5::create_group(zone, "GridCoordinates");

    const NDArray& points = rMesh.Points();
    const std::size_t n = rMesh.NumPoints();
    const std::size_t d = rMesh.PointDim();

    const char* names[3] = {"CoordinateX", "CoordinateY", "CoordinateZ"};
    for (int c = 0; c < 3; ++c) {
        h5::Hid g = h5::create_group(coords, names[c]);
        NDArray col(points.Dtype(), {n});
        for (std::size_t i = 0; i < n; ++i) {
            double v =
                (static_cast<std::size_t>(c) < d) ? detail::read_double(points, i * d + c) : 0.0;
            if (col.Dtype() == DType::Float32)
                col.As<float>()[i] = static_cast<float>(v);
            else
                col.As<double>()[i] = v;
        }
        h5::write_dataset(g, " data", col, gzip_level);
    }

    h5::Hid elems = h5::create_group(zone, "GridElements");
    h5::Hid rng = h5::create_group(elems, "ElementRange");
    h5::Hid conn = h5::create_group(elems, "ElementConnectivity");
    if (tet) {
        const NDArray& tconn = tet->Conn();
        const std::size_t nc = tet->NumCells();
        const std::size_t k = detail::cols(tconn);
        NDArray range(DType::Int64, {2});
        range.As<std::int64_t>()[0] = 1;
        range.As<std::int64_t>()[1] = static_cast<std::int64_t>(nc);
        h5::write_dataset(rng, " data", range, gzip_level);

        NDArray flat(tconn.Dtype(), {nc * k});
        for (std::size_t i = 0; i < nc * k; ++i) {
            std::int64_t v = detail::read_int(tconn, i) + 1;
            switch (flat.Dtype()) {
                case DType::Int32:
                    flat.As<std::int32_t>()[i] = static_cast<std::int32_t>(v);
                    break;
                case DType::Int64:
                    flat.As<std::int64_t>()[i] = v;
                    break;
                case DType::UInt32:
                    flat.As<std::uint32_t>()[i] = static_cast<std::uint32_t>(v);
                    break;
                case DType::UInt64:
                    flat.As<std::uint64_t>()[i] = static_cast<std::uint64_t>(v);
                    break;
                default:
                    throw WriteError("CGNS: unexpected connectivity dtype");
            }
        }
        h5::write_dataset(conn, " data", flat, gzip_level);
    }
}

}  // namespace meshioplusplus

#endif  // MESHIOPLUSPLUS_HAS_HDF5
