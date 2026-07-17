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
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <format>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

// External includes
#ifdef MESHIOPLUSPLUS_HAS_EIGEN
#include <Eigen/Dense>
#endif

// Project includes
#include "meshioplusplus/formats/med.hpp"
#include "meshioplusplus/detail/hdf5_util.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/types.hpp"

namespace meshioplusplus {

namespace {

const std::unordered_map<std::string, std::string>& meshio_to_med() {
    static const std::unordered_map<std::string, std::string> m = {
        {"vertex", "PO1"},       {"line", "SE2"},      {"line3", "SE3"},     {"triangle", "TR3"},
        {"triangle6", "TR6"},    {"triangle7", "TR7"}, {"quad", "QU4"},      {"quad8", "QU8"},
        {"quad9", "QU9"},        {"tetra", "TE4"},     {"tetra10", "T10"},   {"hexahedron", "HE8"},
        {"hexahedron20", "H20"}, {"pyramid", "PY5"},   {"pyramid13", "P13"}, {"wedge", "PE6"},
        {"wedge15", "P15"},      {"polygon", "POG"},   {"polygon2", "POG2"}};
    return m;
}

// Quadratic 3D types share the meshio <-> MED orientation difference, but
// their permutations are not implemented; warn (like the Python reference)
// when reading or writing them unconverted.
void warn_unconverted_3d(const std::string& rCellType) {
    if (rCellType == "tetra10" || rCellType == "hexahedron20" || rCellType == "pyramid13" ||
        rCellType == "wedge15") {
        log::warn(
            "MED: orientation conversion for quadratic 3D cells '{}' is not yet "
            "implemented. These cells may be mis-oriented for MED tools (Salome, "
            "code_saturne, code_aster, etc.).",
            rCellType);
    }
}

// self-inverse meshio <-> MED node permutations (linear 3D types).
const std::unordered_map<std::string, std::vector<int>>& med_node_perm() {
    static const std::unordered_map<std::string, std::vector<int>> m = {
        {"tetra", {0, 1, 3, 2}},
        {"pyramid", {0, 3, 2, 1, 4}},
        {"wedge", {3, 4, 5, 0, 1, 2}},
        {"hexahedron", {4, 5, 6, 7, 0, 1, 2, 3}}};
    return m;
}

// (The former reorder_med_cells pass is fused into flatten_f/unflatten_f via
// their optional `perm` argument — one pass instead of two on both read+write.)

const std::unordered_map<std::string, std::string>& med_to_meshio() {
    static const std::unordered_map<std::string, std::string> m = [] {
        std::unordered_map<std::string, std::string> out;
        for (const auto& kv : meshio_to_med())
            out.emplace(kv.second, kv.first);
        return out;
    }();
    return m;
}

// Fixed-length (h5py np.bytes_-style) string attribute.
void write_attr_bytes(hid_t loc, const std::string& rName, const std::string& rValue) {
    h5::Hid t(H5Tcopy(H5T_C_S1), H5Tclose);
    H5Tset_size(t, std::max<std::size_t>(1, rValue.size()));
    H5Tset_strpad(t, H5T_STR_NULLPAD);
    h5::Hid space(H5Screate(H5S_SCALAR), H5Sclose);
    h5::Hid a(H5Acreate2(loc, rName.c_str(), t, space, H5P_DEFAULT, H5P_DEFAULT), H5Aclose);
    if (!a.Valid())
        throw WriteError(std::format("MED: could not create attribute {}", rName));
    std::string buf = rValue.empty() ? std::string(1, '\0') : rValue;
    H5Awrite(a, t, buf.data());
}

void write_attr_double(hid_t loc, const std::string& rName, double v) {
    h5::Hid space(H5Screate(H5S_SCALAR), H5Sclose);
    h5::Hid a(H5Acreate2(loc, rName.c_str(), H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT),
              H5Aclose);
    H5Awrite(a, H5T_NATIVE_DOUBLE, &v);
}

// Fortran-order (n, k) -> flat column-major buffer, applying `shift` to
// integer dtypes and (fused, same pass) an optional column permutation `perm`
// (the meshio->MED node reorder). Pure index transpose (memory-bandwidth bound).
NDArray flatten_f(const NDArray& rA, std::int64_t shift, const std::vector<int>* pPerm = nullptr) {
    const std::size_t n = detail::rows(rA);
    const std::size_t k = detail::cols(rA);
    const int* p = (pPerm && pPerm->size() == k) ? pPerm->data() : nullptr;
    NDArray out(rA.Dtype(), {n * k});
    detail::dispatch_dtype(rA.Dtype(), [&]<class T>() {
        const T* src = rA.As<T>();
        T* dst = out.As<T>();
#ifdef MESHIOPLUSPLUS_HAS_EIGEN
        if (!p && shift == 0) {
            // (n,k) row-major -> (n,k) col-major = Eigen storage-order convert.
            using RM = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
            using CM = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
            Eigen::Map<CM>(dst, n, k) = Eigen::Map<const RM>(src, n, k);
            return;
        }
#endif
        const T s = static_cast<T>(shift);
        parallel_for_bw(n, [&](std::size_t i) {
            for (std::size_t c = 0; c < k; ++c) {
                std::size_t sc = p ? static_cast<std::size_t>(p[c]) : c;
                if constexpr (std::is_floating_point_v<T>)
                    dst[c * n + i] = src[i * k + sc];
                else
                    dst[c * n + i] = static_cast<T>(src[i * k + sc] + s);
            }
        });
    });
    return out;
}

// Flat column-major buffer -> (n, k) row-major, applying `shift` to integer
// dtypes and (fused, in the same pass) an optional column permutation `perm`
// (the MED->meshio node reorder). Inverse transpose of flatten_f.
NDArray unflatten_f(const NDArray& rFlat, std::size_t n, std::size_t k, std::int64_t shift,
                    const std::vector<int>* pPerm = nullptr) {
    const int* p = (pPerm && pPerm->size() == k) ? pPerm->data() : nullptr;
    NDArray out(rFlat.Dtype(), {n, k});
    detail::dispatch_dtype(rFlat.Dtype(), [&]<class T>() {
        const T* src = rFlat.As<T>();
        T* dst = out.As<T>();
#ifdef MESHIOPLUSPLUS_HAS_EIGEN
        if (!p && shift == 0) {
            using RM = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
            using CM = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
            Eigen::Map<RM>(dst, n, k) = Eigen::Map<const CM>(src, n, k);
            return;
        }
#endif
        const T s = static_cast<T>(shift);
        parallel_for_bw(n, [&](std::size_t i) {
            for (std::size_t c = 0; c < k; ++c) {
                std::size_t sc = p ? static_cast<std::size_t>(p[c]) : c;
                if constexpr (std::is_floating_point_v<T>)
                    dst[i * k + c] = src[sc * n + i];
                else
                    dst[i * k + c] = static_cast<T>(src[sc * n + i] + s);
            }
        });
    });
    return out;
}

constexpr const char* kProfile = "MED_NO_PROFILE_INTERNAL";

// ---- families (point/cell tags) ----

void read_families(hid_t fas_group, std::map<std::int64_t, std::vector<std::string>>& rFamilies,
                   std::map<std::int64_t, std::string>& rGroupNames) {
    for (const std::string& fam_name : h5::group_links(fas_group)) {
        h5::Hid fam = h5::open_group(fas_group, fam_name);
        std::int64_t set_id = h5::read_attr_int(fam, "NUM");
        rGroupNames[set_id] = fam_name;
        if (!h5::exists(fam, "GRO")) {
            rFamilies[set_id] = {};
            continue;
        }
        h5::Hid gro = h5::open_group(fam, "GRO");
        std::int64_t n_subsets = h5::read_attr_int(gro, "NBR");
        NDArray nom = h5::read_dataset(gro, "NOM");  // (n_subsets, 80) int8
        std::vector<std::string> names;
        for (std::int64_t i = 0; i < n_subsets; ++i) {
            std::string s;
            for (int c = 0; c < 80; ++c) {
                char ch = static_cast<char>(detail::read_int(nom, i * 80 + c));
                if (ch == '\0')
                    break;
                s += ch;
            }
            std::size_t b = s.find_first_not_of(' ');
            std::size_t e = s.find_last_not_of(' ');
            names.push_back(b == std::string::npos ? std::string() : s.substr(b, e - b + 1));
        }
        rFamilies.emplace(set_id, std::move(names));
    }
}

// Read a fixed-length string attribute (latin-1), stripped of spaces and NULs.
std::string read_attr_bytes(hid_t loc, const std::string& rName) {
    if (!h5::has_attr(loc, rName))
        return "";
    std::string s = h5::read_attr_string(loc, rName);
    // strip trailing NULs and surrounding spaces
    std::size_t z = s.find('\0');
    if (z != std::string::npos)
        s = s.substr(0, z);
    std::size_t b = s.find_first_not_of(' ');
    if (b == std::string::npos)
        return "";
    std::size_t e = s.find_last_not_of(' ');
    return s.substr(b, e - b + 1);
}

// Matches _write_families in _med.py: family link name from `group_names`
// (else "FAM_<id>_"), '/'->'_', capped at 64 bytes -> "FAM_<id>"; no GRO
// subgroup when the family has no named groups; GRO/NOM is an
// H5T_ARRAY{[80] char} dataset, one 80-char slot per name, space-padded.
void write_families(hid_t fm_group, const std::map<std::int64_t, std::vector<std::string>>& rTags,
                    const std::map<std::int64_t, std::string>& rGroupNames) {
    for (const auto& kv : rTags) {
        std::int64_t set_id = kv.first;
        const std::vector<std::string>& names = kv.second;
        auto git = rGroupNames.find(set_id);
        std::string gname =
            git != rGroupNames.end() ? git->second : ("FAM_" + std::to_string(set_id) + "_");
        for (char& c : gname)
            if (c == '/')
                c = '_';
        if (gname.size() > 64)
            gname = "FAM_" + std::to_string(set_id);

        h5::Hid family = h5::create_group(fm_group, gname);
        h5::write_attr_int(family, "NUM", set_id);
        if (names.empty())
            continue;

        h5::Hid gro = h5::create_group(family, "GRO");
        h5::write_attr_int(gro, "NBR", static_cast<std::int64_t>(names.size()));
        hsize_t n = names.size(), eighty = 80;
        h5::Hid at(H5Tarray_create2(H5T_STD_I8LE, 1, &eighty), H5Tclose);
        h5::Hid mt(H5Tarray_create2(H5T_NATIVE_INT8, 1, &eighty), H5Tclose);
        h5::Hid space(H5Screate_simple(1, &n, nullptr), H5Sclose);
        h5::Hid d(H5Dcreate2(gro, "NOM", at, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),
                  H5Dclose);
        std::vector<std::int8_t> buf(names.size() * 80, static_cast<std::int8_t>(' '));
        for (std::size_t i = 0; i < names.size(); ++i) {
            if (names[i].size() > 80)
                throw WriteError(std::format(
                    "Family name '{}' is too long for MED format (max 80 bytes).", names[i]));
            for (std::size_t c = 0; c < names[i].size(); ++c)
                buf[i * 80 + c] = static_cast<std::int8_t>(names[i][c]);
        }
        H5Dwrite(d, mt, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data());
    }
}

}  // namespace

Mesh read_med(const std::string& rPath, MedInfo& rInfo) {
    h5::SilenceErrors silence;
    h5::Hid f = h5::open_file_read(rPath);

    h5::Hid ens = h5::open_group(f, "ENS_MAA");
    std::vector<std::string> meshes = h5::group_links(ens);
    if (meshes.size() != 1)
        throw ReadError(std::format("Must only contain exactly 1 mesh, found {}.", meshes.size()));
    const std::string mesh_name = meshes[0];
    h5::Hid mesh_grp = h5::open_group(ens, mesh_name);

    std::int64_t dim = h5::read_attr_int(mesh_grp, "ESP");

    // Mesh-level metadata attributes.
    rInfo.mMeshName = mesh_name;
    rInfo.mDescription = read_attr_bytes(mesh_grp, "DES");
    rInfo.mUnitTime = read_attr_bytes(mesh_grp, "UNT");
    rInfo.mUnitCoords = read_attr_bytes(mesh_grp, "UNI");

    // Possible time-stepping indirection.
    h5::Hid data_grp;
    if (h5::exists(mesh_grp, "NOE")) {
        data_grp = std::move(mesh_grp);
    } else {
        std::vector<std::string> steps = h5::group_links(mesh_grp);
        if (steps.size() != 1)
            throw ReadError(
                std::format("Must only contain exactly 1 time-step, found {}.", steps.size()));
        data_grp = h5::open_group(mesh_grp, steps[0]);
    }

    Mesh mesh;

    // Points
    h5::Hid noe = h5::open_group(data_grp, "NOE");
    {
        h5::Hid coo_ds(H5Dopen2(noe, "COO", H5P_DEFAULT), H5Dclose);
        if (!coo_ds.Valid())
            throw ReadError("MED: missing NOE/COO");
        std::int64_t n_points = h5::read_attr_int(coo_ds, "NBR");
        NDArray coo = h5::read_dataset(noe, "COO");
        mesh.AssignPoints(
            unflatten_f(coo, static_cast<std::size_t>(n_points), static_cast<std::size_t>(dim), 0));
    }

    // Point tags
    if (h5::exists(noe, "FAM"))
        mesh.AddPointData("point_tags", h5::read_dataset(noe, "FAM"));

    // Families info
    h5::Hid fas = h5::exists(data_grp, "FAS") ? h5::open_group(data_grp, "FAS") : h5::Hid();
    if (!fas.Valid()) {
        h5::Hid fas_root = h5::open_group(f, "FAS");
        fas = h5::open_group(fas_root, mesh_name);
    }
    if (h5::exists(fas, "NOEUD")) {
        h5::Hid noeud = h5::open_group(fas, "NOEUD");
        read_families(noeud, rInfo.mPointTags, rInfo.mPointTagGroups);
    }

    // Cells
    std::vector<std::string> cell_types;  // meshio names, in read order
    h5::Hid mai = h5::open_group(data_grp, "MAI");
    std::vector<NDArray> cell_tag_blocks;
    bool any_cell_tags = false;
    // Cell-block order is significant (aligns cell_data / cell_sets); iterate in
    // HDF5 creation order to match the Python (h5py track_order) reader.
    for (const std::string& med_type : h5::group_links_crt(mai)) {
        auto it = med_to_meshio().find(med_type);
        if (it == med_to_meshio().end())
            throw ReadError(std::format("MED: unsupported cell type {}", med_type));
        h5::Hid g = h5::open_group(mai, med_type);

        if (med_type == "POG" || med_type == "POG2") {
            // Ragged polygons: flat 1-based NOD + 1-based INN offsets.
            NDArray nod = h5::read_dataset(g, "NOD");
            NDArray inn = h5::read_dataset(g, "INN");
            std::size_t npoly = inn.Size() > 0 ? inn.Size() - 1 : 0;
            std::vector<std::vector<std::int64_t>> rows;
            for (std::size_t i = 0; i < npoly; ++i) {
                std::int64_t a = detail::read_int(inn, i) - 1;
                std::int64_t b = detail::read_int(inn, i + 1) - 1;
                std::vector<std::int64_t> row;
                for (std::int64_t j = a; j < b; ++j)
                    row.push_back(detail::read_int(nod, static_cast<std::size_t>(j)) - 1);
                rows.push_back(std::move(row));
            }
            mesh.AddPolygonBlock(it->second, std::move(rows));
            cell_types.push_back(it->second);
        } else {
            h5::Hid nod_ds(H5Dopen2(g, "NOD", H5P_DEFAULT), H5Dclose);
            if (!nod_ds.Valid())
                throw ReadError(std::format("MED: missing NOD for {}", med_type));
            std::int64_t n_cells = h5::read_attr_int(nod_ds, "NBR");
            NDArray nod = h5::read_dataset(g, "NOD");
            std::size_t k = n_cells > 0 ? nod.Size() / static_cast<std::size_t>(n_cells) : 0;
            warn_unconverted_3d(it->second);
            // Fuse the Fortran->C transpose (shift -1) with the MED->meshio
            // node reorder into a single pass over the connectivity.
            auto pit = med_node_perm().find(it->second);
            const std::vector<int>* perm =
                (pit != med_node_perm().end() && pit->second.size() == k) ? &pit->second : nullptr;
            NDArray data = unflatten_f(nod, static_cast<std::size_t>(n_cells), k, -1, perm);
            mesh.AddCellBlock(it->second, std::move(data));
            cell_types.push_back(it->second);
        }

        if (h5::exists(g, "FAM")) {
            cell_tag_blocks.push_back(h5::read_dataset(g, "FAM"));
            any_cell_tags = true;
        }
    }
    if (any_cell_tags) {
        if (cell_tag_blocks.size() != mesh.NumCellBlocks())
            throw ReadError("MED: partial cell tags handled by Python fallback");
        mesh.AddCellData("cell_tags", std::move(cell_tag_blocks));
    }

    if (h5::exists(fas, "ELEME")) {
        h5::Hid eleme = h5::open_group(fas, "ELEME");
        read_families(eleme, rInfo.mCellTags, rInfo.mCellTagGroups);
    }

    // Fields (CHA): the enhanced Python reader attaches med:field_units /
    // med:step_meta and multi-timestep metadata that the C++ path does not
    // replicate byte-for-byte; defer any field-carrying file to Python.
    if (h5::exists(f, "CHA"))
        throw ReadError("MED: fields (CHA) handled by Python fallback");

    return mesh;
}

void write_med(const std::string& rPath, const Mesh& rMesh, const MedInfo& rInfo,
               const std::string& rMedVersion) {
    h5::SilenceErrors silence;

    // Fields (CHA) with the MED-4.1 bitmask / units / step metadata and the
    // gmsh:physical family bridging are produced by the enhanced Python writer
    // and inspected byte-for-byte by tests; defer any such mesh to Python.
    for (const auto& name : rMesh.PointDataNames())
        if (name != "point_tags")
            throw WriteError("MED: fields handled by Python fallback");
    for (const auto& name : rMesh.CellDataNames())
        if (name != "cell_tags")
            throw WriteError("MED: fields handled by Python fallback");
    if (rMesh.HasCellData("gmsh:physical"))
        throw WriteError("MED: gmsh physical groups handled by Python fallback");

    // MED cannot have two blocks of the same type.
    for (std::size_t i = 0; i < rMesh.NumCellBlocks(); ++i)
        for (std::size_t j = i + 1; j < rMesh.NumCellBlocks(); ++j)
            if (rMesh.Cells(i).Type() == rMesh.Cells(j).Type())
                throw WriteError("MED files cannot have two sections of the same cell type.");

    // Parse med_version -> MAJ.MIN.REL (default 4.1.0 on error).
    int maj = 4, min = 1, rel = 0;
    {
        int parts[3] = {4, 1, 0};
        std::size_t start = 0, idx = 0;
        bool ok = true;
        for (idx = 0; idx < 3; ++idx) {
            std::size_t dot = rMedVersion.find('.', start);
            std::string tok = rMedVersion.substr(
                start, dot == std::string::npos ? std::string::npos : dot - start);
            try {
                parts[idx] = std::stoi(tok);
            } catch (...) {
                ok = false;
                break;
            }
            if (dot == std::string::npos)
                break;
            start = dot + 1;
        }
        if (ok) {
            maj = parts[0];
            min = parts[1];
            rel = parts[2];
        }
    }

    h5::Hid f = h5::create_file(rPath);

    h5::Hid infos = h5::create_group(f, "INFOS_GENERALES");
    h5::write_attr_int(infos, "MAJ", maj);
    h5::write_attr_int(infos, "MIN", min);
    h5::write_attr_int(infos, "REL", rel);

    const std::string mesh_name = rInfo.mMeshName.empty() ? "mesh" : rInfo.mMeshName;
    const std::size_t dim = rMesh.PointDim();

    h5::Hid ens = h5::create_group(f, "ENS_MAA");
    h5::Hid med_mesh = h5::create_group(ens, mesh_name);
    h5::write_attr_int(med_mesh, "DIM", static_cast<std::int64_t>(dim));
    h5::write_attr_int(med_mesh, "ESP", static_cast<std::int64_t>(dim));
    h5::write_attr_int(med_mesh, "REP", 0);
    write_attr_bytes(med_mesh, "UNT", rInfo.mUnitTime);
    write_attr_bytes(med_mesh, "UNI", rInfo.mUnitCoords);
    h5::write_attr_int(med_mesh, "SRT", 1);
    {
        const char* names[3] = {"X", "Y", "Z"};
        std::string nom;
        for (std::size_t c = 0; c < dim && c < 3; ++c) {
            char buf[20];
            std::snprintf(buf, sizeof(buf), "%-16s", names[c]);
            nom += buf;
        }
        write_attr_bytes(med_mesh, "NOM", nom);
    }
    write_attr_bytes(
        med_mesh, "DES",
        rInfo.mDescription.empty() ? "Mesh created with meshio++" : rInfo.mDescription);
    h5::write_attr_int(med_mesh, "TYP", 0);

    h5::Hid time_step = h5::create_group(med_mesh, "-0000000000000000001-0000000000000000001");
    h5::write_attr_int(time_step, "CGT", 1);
    h5::write_attr_int(time_step, "NDT", -1);
    h5::write_attr_int(time_step, "NOR", -1);
    write_attr_double(time_step, "PDT", -1.0);

    // Points
    h5::Hid noe = h5::create_group(time_step, "NOE");
    h5::write_attr_int(noe, "CGT", 1);
    h5::write_attr_int(noe, "CGS", 1);
    write_attr_bytes(noe, "PFL", kProfile);
    {
        NDArray coo = flatten_f(rMesh.Points(), 0);
        h5::write_dataset(noe, "COO", coo);
        h5::Hid d(H5Dopen2(noe, "COO", H5P_DEFAULT), H5Dclose);
        h5::write_attr_int(d, "CGT", 1);
        h5::write_attr_int(d, "NBR", static_cast<std::int64_t>(rMesh.NumPoints()));
    }
    if (rMesh.HasPointData("point_tags")) {
        h5::write_dataset(noe, "FAM", rMesh.PointData("point_tags"));
        h5::Hid d(H5Dopen2(noe, "FAM", H5P_DEFAULT), H5Dclose);
        h5::write_attr_int(d, "CGT", 1);
        h5::write_attr_int(d, "NBR", static_cast<std::int64_t>(rMesh.NumPoints()));
    }

    // Cells
    h5::Hid mai = h5::create_group(time_step, "MAI");
    h5::write_attr_int(mai, "CGT", 1);
    const bool has_cell_tags = rMesh.HasCellData("cell_tags");
    for (std::size_t k = 0; k < rMesh.NumCellBlocks(); ++k) {
        const auto cb = rMesh.Cells(k);
        auto it = meshio_to_med().find(cb.Type());
        if (it == meshio_to_med().end())
            throw WriteError(std::format("MED: unsupported cell type {}", cb.Type()));
        h5::Hid g = h5::create_group(mai, it->second);
        h5::write_attr_int(g, "CGT", 1);
        h5::write_attr_int(g, "CGS", 1);
        write_attr_bytes(g, "PFL", kProfile);

        if (cb.Type() == "polygon" || cb.Type() == "polygon2") {
            // Ragged: flat 1-based NOD + 1-based INN offsets.
            std::vector<std::int64_t> nod;
            std::vector<std::int64_t> inn = {1};
            for (std::size_t i = 0; i < cb.NumCells(); ++i) {
                const std::int64_t* row = cb.Row(i);
                const std::size_t row_size = cb.RowSize(i);
                for (std::size_t j = 0; j < row_size; ++j)
                    nod.push_back(row[j] + 1);
                inn.push_back(inn.back() + static_cast<std::int64_t>(row_size));
            }
            NDArray nod_a(DType::Int64, {nod.size()});
            for (std::size_t i = 0; i < nod.size(); ++i)
                nod_a.As<std::int64_t>()[i] = nod[i];
            NDArray inn_a(DType::Int64, {inn.size()});
            for (std::size_t i = 0; i < inn.size(); ++i)
                inn_a.As<std::int64_t>()[i] = inn[i];
            h5::write_dataset(g, "NOD", nod_a);
            h5::write_dataset(g, "INN", inn_a);
            h5::Hid d(H5Dopen2(g, "NOD", H5P_DEFAULT), H5Dclose);
            h5::write_attr_int(d, "CGT", 1);
            h5::write_attr_int(d, "NBR", static_cast<std::int64_t>(cb.NumCells()));
        } else {
            warn_unconverted_3d(cb.Type());
            // Fuse the meshio->MED node reorder with the Fortran transpose
            // (shift +1) into a single pass (mirrors the read side).
            auto pit = med_node_perm().find(cb.Type());
            const std::vector<int>* perm = (pit != med_node_perm().end()) ? &pit->second : nullptr;
            NDArray nod = flatten_f(cb.Conn(), +1, perm);
            h5::write_dataset(g, "NOD", nod);
            h5::Hid d(H5Dopen2(g, "NOD", H5P_DEFAULT), H5Dclose);
            h5::write_attr_int(d, "CGT", 1);
            h5::write_attr_int(d, "NBR", static_cast<std::int64_t>(cb.NumCells()));
        }
        if (has_cell_tags && k < rMesh.CellDataNumBlocks("cell_tags")) {
            h5::write_dataset(g, "FAM", rMesh.CellData("cell_tags", k));
            h5::Hid d(H5Dopen2(g, "FAM", H5P_DEFAULT), H5Dclose);
            h5::write_attr_int(d, "CGT", 1);
            h5::write_attr_int(d, "NBR", static_cast<std::int64_t>(cb.NumCells()));
        }
    }

    // Families
    h5::Hid fas = h5::create_group(f, "FAS");
    h5::Hid families = h5::create_group(fas, mesh_name);
    h5::Hid family_zero = h5::create_group(families, "FAMILLE_ZERO");
    h5::write_attr_int(family_zero, "NUM", 0);
    if (!rInfo.mPointTags.empty()) {
        h5::Hid node = h5::create_group(families, "NOEUD");
        write_families(node, rInfo.mPointTags, rInfo.mPointTagGroups);
    }
    if (!rInfo.mCellTags.empty()) {
        h5::Hid element = h5::create_group(families, "ELEME");
        write_families(element, rInfo.mCellTags, rInfo.mCellTagGroups);
    }
}

}  // namespace meshioplusplus

#endif  // MESHIOPLUSPLUS_HAS_HDF5
