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
#include <ctime>
#include <cstdint>
#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>

// Project includes
#include "meshioplusplus/formats/h5m.hpp"
#include "meshioplusplus/detail/hdf5_util.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus {

namespace {

const std::unordered_map<std::string, std::string>& h5m_to_meshio() {
    static const std::unordered_map<std::string, std::string> m = {
        {"Edge2", "line"}, {"Hex8", "hexahedron"}, {"Prism6", "wedge"}, {"Pyramid5", "pyramid"},
        {"Quad4", "quad"}, {"Tri3", "triangle"},   {"Tet4", "tetra"}};
    return m;
}

// The MOAB element-type enum (h5py's special_dtype(enum=...)).
h5::Hid make_elem_enum() {
    h5::Hid t(H5Tenum_create(H5T_NATIVE_INT), H5Tclose);
    const std::pair<const char*, int> members[] = {
        {"Edge", 1},    {"Tri", 2},   {"Quad", 3},  {"Polygon", 4}, {"Tet", 5},
        {"Pyramid", 6}, {"Prism", 7}, {"Knife", 8}, {"Hex", 9},     {"Polyhedron", 10}};
    for (const auto& mv : members) {
        int v = mv.second;
        H5Tenum_insert(t, mv.first, &v);
    }
    return t;
}

// Fixed-length byte-string dataset (h5py's data=[b"...", ...]).
void write_history(hid_t loc, int gzip_level) {
    std::time_t now = std::time(nullptr);
    char stamp[64];
    std::strftime(stamp, sizeof(stamp), "%Y-%m-%d %H:%M:%S", std::localtime(&now));
    std::vector<std::string> items = {"meshioplusplus.h5m", "cpp-core", stamp};

    std::size_t maxlen = 1;
    for (const auto& s : items)
        maxlen = std::max(maxlen, s.size());
    h5::Hid st(H5Tcopy(H5T_C_S1), H5Tclose);
    H5Tset_size(st, maxlen);
    H5Tset_strpad(st, H5T_STR_NULLPAD);

    hsize_t dims[1] = {items.size()};
    h5::Hid space(H5Screate_simple(1, dims, nullptr), H5Sclose);
    h5::Hid dcpl(H5Pcreate(H5P_DATASET_CREATE), H5Pclose);
    if (gzip_level >= 0) {
        H5Pset_chunk(dcpl, 1, dims);
        H5Pset_deflate(dcpl, static_cast<unsigned>(gzip_level));
    }
    h5::Hid d(H5Dcreate2(loc, "history", st, space, H5P_DEFAULT, dcpl, H5P_DEFAULT), H5Dclose);
    std::vector<char> buf(items.size() * maxlen, '\0');
    for (std::size_t i = 0; i < items.size(); ++i)
        std::memcpy(buf.data() + i * maxlen, items[i].data(), items[i].size());
    H5Dwrite(d, st, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data());
}

// Write a point-data tag: 1-D directly, 2-D as (n,) of k-tuples (array dtype).
void write_tag_dataset(hid_t loc, const std::string& rName, const NDArray& rArr, int gzip_level) {
    if (rArr.Ndim() <= 1) {
        h5::write_dataset(loc, rName, rArr, gzip_level);
        return;
    }
    hsize_t n = rArr.Shape()[0];
    hsize_t k = rArr.Shape()[1];
    h5::Hid ft(H5Tarray_create2(h5::file_type(rArr.Dtype()), 1, &k), H5Tclose);
    h5::Hid mt(H5Tarray_create2(h5::native_type(rArr.Dtype()), 1, &k), H5Tclose);
    h5::Hid space(H5Screate_simple(1, &n, nullptr), H5Sclose);
    h5::Hid dcpl(H5Pcreate(H5P_DATASET_CREATE), H5Pclose);
    if (gzip_level >= 0 && n > 0) {
        H5Pset_chunk(dcpl, 1, &n);
        H5Pset_deflate(dcpl, static_cast<unsigned>(gzip_level));
    }
    h5::Hid d(H5Dcreate2(loc, rName.c_str(), ft, space, H5P_DEFAULT, dcpl, H5P_DEFAULT), H5Dclose);
    if (n > 0)
        H5Dwrite(d, mt, H5S_ALL, H5S_ALL, H5P_DEFAULT, rArr.Data());
}

}  // namespace

Mesh read_h5m(const std::string& rPath) {
    h5::SilenceErrors silence;
    h5::Hid f = h5::open_file_read(rPath);
    h5::Hid tstt = h5::open_group(f, "tstt");

    Mesh mesh;
    h5::Hid nodes = h5::open_group(tstt, "nodes");
    mesh.AssignPoints(h5::read_dataset(nodes, "coordinates"));

    if (h5::exists(nodes, "tags")) {
        h5::Hid tags = h5::open_group(nodes, "tags");
        for (const std::string& name : h5::group_links(tags))
            mesh.AddPointData(name, h5::read_dataset(tags, name));
    }

    if (h5::exists(tstt, "elements")) {
        h5::Hid elements = h5::open_group(tstt, "elements");
        for (const std::string& h5m_type : h5::group_links(elements)) {
            auto it = h5m_to_meshio().find(h5m_type);
            if (it == h5m_to_meshio().end())
                throw ReadError("H5M: unknown element type " + h5m_type);
            h5::Hid g = h5::open_group(elements, h5m_type);
            NDArray conn = h5::read_dataset(g, "connectivity");
            // h5m indices are 1-based.
            for (std::size_t i = 0; i < conn.Size(); ++i) {
                switch (conn.Dtype()) {
                    case DType::Int32:
                        conn.As<std::int32_t>()[i] -= 1;
                        break;
                    case DType::Int64:
                        conn.As<std::int64_t>()[i] -= 1;
                        break;
                    case DType::UInt32:
                        conn.As<std::uint32_t>()[i] -= 1;
                        break;
                    case DType::UInt64:
                        conn.As<std::uint64_t>()[i] -= 1;
                        break;
                    default:
                        throw ReadError("H5M: unexpected connectivity dtype");
                }
            }
            mesh.AddCellBlock(it->second, std::move(conn));
        }
    }
    // Element tags (cell data) and sets are not read (matching the Python reader).

    return mesh;
}

void write_h5m(const std::string& rPath, const Mesh& rMesh, bool add_global_ids, int gzip_level) {
    h5::SilenceErrors silence;
    h5::Hid f = h5::create_file(rPath);
    h5::Hid tstt = h5::create_group(f, "tstt");

    std::int64_t global_id = 1;  // h5m base index

    // nodes
    h5::Hid nodes = h5::create_group(tstt, "nodes");
    h5::write_dataset(nodes, "coordinates", rMesh.Points(), gzip_level);
    {
        h5::Hid d(H5Dopen2(nodes, "coordinates", H5P_DEFAULT), H5Dclose);
        h5::write_attr_int(d, "start_id", global_id);
    }
    global_id += static_cast<std::int64_t>(rMesh.NumPoints());

    h5::Hid tstt_tags = h5::create_group(tstt, "tags");

    // point data (+ auto GLOBAL_ID)
    std::vector<std::pair<std::string, const NDArray*>> pd;
    for (const auto& name : rMesh.PointDataNames())
        pd.emplace_back(name, &rMesh.PointData(name));
    NDArray gids;
    if (add_global_ids && !rMesh.HasPointData("GLOBAL_ID")) {
        gids = NDArray(DType::Int64, {rMesh.NumPoints()});
        for (std::size_t i = 0; i < rMesh.NumPoints(); ++i)
            gids.As<std::int64_t>()[i] = static_cast<std::int64_t>(i) + 1;
        pd.emplace_back("GLOBAL_ID", &gids);
    }

    if (!pd.empty()) {
        h5::Hid tags = h5::create_group(nodes, "tags");
        for (const auto& kv : pd) {
            write_tag_dataset(tags, kv.first, *kv.second, gzip_level);
            // Global tag entry: committed datatype + dense-class attribute.
            h5::Hid g = h5::create_group(tstt_tags, kv.first);
            h5::Hid t = [&]() -> h5::Hid {
                if (kv.second->Ndim() >= 2) {
                    hsize_t k = kv.second->Shape()[1];
                    return h5::Hid(H5Tarray_create2(h5::file_type(kv.second->Dtype()), 1, &k),
                                   H5Tclose);
                }
                return h5::Hid(H5Tcopy(h5::file_type(kv.second->Dtype())), H5Tclose);
            }();
            H5Tcommit2(g, "type", t, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            h5::write_attr_int(g, "class", 2);
        }
    }

    // elements
    h5::Hid elements = h5::create_group(tstt, "elements");
    h5::Hid elem_dt = make_elem_enum();
    H5Tcommit2(tstt, "elemtypes", elem_dt, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    write_history(tstt, gzip_level);

    struct H5mType {
        const char* mName;
        int mType;
    };
    static const std::unordered_map<std::string, H5mType> meshio_to_h5m = {
        {"line", {"Edge2", 1}}, {"triangle", {"Tri3", 2}}, {"tetra", {"Tet4", 5}}};

    for (const auto cb : rMesh.CellRange()) {
        auto it = meshio_to_h5m.find(cb.Type());
        if (it == meshio_to_h5m.end())
            continue;  // unsupported type: skipped with a warning in Python
        h5::Hid g = h5::create_group(elements, it->second.mName);
        {
            h5::Hid space(H5Screate(H5S_SCALAR), H5Sclose);
            h5::Hid a(H5Acreate2(g, "element_type", elem_dt, space, H5P_DEFAULT, H5P_DEFAULT),
                      H5Aclose);
            int v = it->second.mType;
            H5Awrite(a, elem_dt, &v);
        }
        // 1-based connectivity, preserving the integer dtype.
        const NDArray& cconn = cb.Conn();
        NDArray conn(cconn.Dtype(), cconn.Shape());
        for (std::size_t i = 0; i < cconn.Size(); ++i) {
            std::int64_t v = detail::read_int(cconn, i) + 1;
            switch (conn.Dtype()) {
                case DType::Int32:
                    conn.As<std::int32_t>()[i] = static_cast<std::int32_t>(v);
                    break;
                case DType::Int64:
                    conn.As<std::int64_t>()[i] = v;
                    break;
                case DType::UInt32:
                    conn.As<std::uint32_t>()[i] = static_cast<std::uint32_t>(v);
                    break;
                case DType::UInt64:
                    conn.As<std::uint64_t>()[i] = static_cast<std::uint64_t>(v);
                    break;
                default:
                    throw WriteError("H5M: unexpected connectivity dtype");
            }
        }
        h5::write_dataset(g, "connectivity", conn, gzip_level);
        {
            h5::Hid d(H5Dopen2(g, "connectivity", H5P_DEFAULT), H5Dclose);
            h5::write_attr_int(d, "start_id", global_id);
        }
        global_id += static_cast<std::int64_t>(cb.NumCells());
    }
    // Cell data is not written: the Python writer's cell-data path is broken
    // upstream (iterates a list as a dict) and the reader ignores element tags.

    // empty set group -- MOAB wants this
    h5::Hid sets = h5::create_group(tstt, "sets");
    h5::create_group(sets, "tags");

    h5::write_attr_int(tstt, "max_id", global_id, H5T_STD_U64LE);
}

}  // namespace meshioplusplus

#endif  // MESHIOPLUSPLUS_HAS_HDF5
