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
#pragma once

/**
 * @file hdf5_util.hpp
 * @brief Shared low-level HDF5 helpers used by every C++ format that stores
 * data in an HDF5 container (MED, XDMF's `Format="HDF"` DataItems).
 *
 * Provides: an RAII handle wrapper (`Hid`) so every `H5*` resource is closed
 * exactly once even under exceptions; dataset read/write helpers that
 * translate between `meshioplusplus::DType` and HDF5's native/file types
 * (matching what h5py writes on x86: little-endian file types via
 * `file_type()`); scalar/string attribute helpers matching h5py's own
 * variable-length UTF-8 convention; group link listing in both name order
 * and (where the file tracks it) creation order, the latter needed where
 * block order carries meaning (e.g. MED's `MAI` cell blocks, whose order
 * must align with `cell_data`/`cell_sets`); and `SilenceErrors`, which
 * suppresses HDF5's default stderr error-stack printing so failures surface
 * only as the C++ exceptions this codebase converts them to (`ReadError`/
 * `WriteError`). This entire header compiles to nothing when
 * `MESHIOPLUSPLUS_HAS_HDF5` is not defined, i.e. when the build has no HDF5
 * library — the HDF-dependent C++ code paths are then simply absent and
 * callers fall back to the pure-Python (h5py-based) implementation.
 */

#ifdef MESHIOPLUSPLUS_HAS_HDF5

// External includes
#include <hdf5.h>

// System includes
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/ndarray.hpp"

namespace meshioplusplus {
namespace h5 {

/**
 * @brief RAII wrapper for an HDF5 `hid_t` handle, paired with the `H5*Close`
 * function that must release it.
 *
 * Move-only (copying an `hid_t` would double-close it): moving transfers
 * ownership and leaves the source handle invalid (`mId = -1`). Implicitly
 * convertible to `hid_t` so it can be passed straight into `H5*` C API
 * calls. `Valid()` reports whether the handle is currently open (id `>= 0`).
 */
class Hid {
public:
    using Closer = herr_t (*)(hid_t);
    Hid() = default;
    Hid(hid_t id, Closer closer) : mId(id), mCloser(closer) {}
    Hid(Hid&& o) noexcept : mId(o.mId), mCloser(o.mCloser) { o.mId = -1; }
    Hid& operator=(Hid&& o) noexcept {
        Reset();
        mId = o.mId;
        mCloser = o.mCloser;
        o.mId = -1;
        return *this;
    }
    Hid(const Hid&) = delete;
    Hid& operator=(const Hid&) = delete;
    ~Hid() { Reset(); }

    void Reset() {
        if (mId >= 0 && mCloser)
            mCloser(mId);
        mId = -1;
    }
    bool Valid() const { return mId >= 0; }
    hid_t Get() const { return mId; }
    operator hid_t() const { return mId; }

private:
    hid_t mId = -1;
    Closer mCloser = nullptr;
};

/**
 * @brief Opens an existing HDF5 file read-only.
 * @param rPath Filesystem path of the file to open.
 * @return Owning `Hid` for the open file.
 * @throws ReadError if the file cannot be opened.
 */
inline Hid open_file_read(const std::string& rPath) {
    Hid f(H5Fopen(rPath.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose);
    if (!f.Valid())
        throw ReadError("HDF5: could not open file " + rPath);
    return f;
}

/**
 * @brief Creates a new HDF5 file, truncating any existing file at `path`.
 * @param rPath Filesystem path of the file to create.
 * @return Owning `Hid` for the new file.
 * @throws WriteError if the file cannot be created.
 */
inline Hid create_file(const std::string& rPath) {
    Hid f(H5Fcreate(rPath.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT), H5Fclose);
    if (!f.Valid())
        throw WriteError("HDF5: could not create file " + rPath);
    return f;
}

/**
 * @brief Whether a link named `name` exists directly under group/file `loc`.
 * @param loc Group or file handle to look under.
 * @param rName Link name to test.
 * @return `true` if the link exists.
 */
inline bool exists(hid_t loc, const std::string& rName) {
    return H5Lexists(loc, rName.c_str(), H5P_DEFAULT) > 0;
}

/**
 * @brief Opens an existing HDF5 group.
 * @param loc Parent group or file handle.
 * @param rName Name of the group to open.
 * @return Owning `Hid` for the opened group.
 * @throws ReadError if the group does not exist.
 */
inline Hid open_group(hid_t loc, const std::string& rName) {
    Hid g(H5Gopen2(loc, rName.c_str(), H5P_DEFAULT), H5Gclose);
    if (!g.Valid())
        throw ReadError("HDF5: missing group '" + rName + "'");
    return g;
}

/**
 * @brief Creates a new HDF5 group.
 * @param loc Parent group or file handle.
 * @param rName Name of the group to create.
 * @return Owning `Hid` for the new group.
 * @throws WriteError if the group cannot be created.
 */
inline Hid create_group(hid_t loc, const std::string& rName) {
    Hid g(H5Gcreate2(loc, rName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), H5Gclose);
    if (!g.Valid())
        throw WriteError("HDF5: could not create group '" + rName + "'");
    return g;
}

/**
 * @brief Maps a `meshioplusplus::DType` to the native in-memory HDF5 type used
 * for `H5Dread`/`H5Dwrite` (host byte order/representation, not the
 * on-disk file type — see `file_type()` for that).
 * @param dt The dtype to convert.
 * @return The matching `H5T_NATIVE_*` constant (defaults to
 *         `H5T_NATIVE_DOUBLE` for an unrecognized/invalid `dt`).
 */
inline hid_t native_type(DType dt) {
    switch (dt) {
        case DType::Float32:
            return H5T_NATIVE_FLOAT;
        case DType::Float64:
            return H5T_NATIVE_DOUBLE;
        case DType::Int8:
            return H5T_NATIVE_INT8;
        case DType::Int16:
            return H5T_NATIVE_INT16;
        case DType::Int32:
            return H5T_NATIVE_INT32;
        case DType::Int64:
            return H5T_NATIVE_INT64;
        case DType::UInt8:
            return H5T_NATIVE_UINT8;
        case DType::UInt16:
            return H5T_NATIVE_UINT16;
        case DType::UInt32:
            return H5T_NATIVE_UINT32;
        case DType::UInt64:
            return H5T_NATIVE_UINT64;
    }
    return H5T_NATIVE_DOUBLE;
}

/**
 * @brief Maps a `meshioplusplus::DType` to the on-disk (file) HDF5 type to use
 * when creating a dataset/attribute.
 *
 * Always little-endian (`H5T_*LE`), matching what h5py writes on x86, so
 * files produced by the C++ writer are byte-for-byte compatible with the
 * pure-Python/h5py writer's output.
 * @param dt The dtype to convert.
 * @return The matching `H5T_*LE` constant (defaults to `H5T_IEEE_F64LE`).
 */
inline hid_t file_type(DType dt) {
    switch (dt) {
        case DType::Float32:
            return H5T_IEEE_F32LE;
        case DType::Float64:
            return H5T_IEEE_F64LE;
        case DType::Int8:
            return H5T_STD_I8LE;
        case DType::Int16:
            return H5T_STD_I16LE;
        case DType::Int32:
            return H5T_STD_I32LE;
        case DType::Int64:
            return H5T_STD_I64LE;
        case DType::UInt8:
            return H5T_STD_U8LE;
        case DType::UInt16:
            return H5T_STD_U16LE;
        case DType::UInt32:
            return H5T_STD_U32LE;
        case DType::UInt64:
            return H5T_STD_U64LE;
    }
    return H5T_IEEE_F64LE;
}

/**
 * @brief Converts a stored HDF5 datatype (of a dataset or attribute) to the
 * corresponding `meshioplusplus::DType`.
 * @param type_id HDF5 type id, as returned e.g. by `H5Dget_type`.
 * @return The matching `DType`.
 * @throws ReadError if `type_id`'s class is neither float nor integer.
 */
inline DType dtype_from_h5(hid_t type_id) {
    H5T_class_t cls = H5Tget_class(type_id);
    std::size_t sz = H5Tget_size(type_id);
    if (cls == H5T_FLOAT)
        return sz == 4 ? DType::Float32 : DType::Float64;
    if (cls == H5T_INTEGER) {
        bool is_signed = H5Tget_sign(type_id) != H5T_SGN_NONE;
        switch (sz) {
            case 1:
                return is_signed ? DType::Int8 : DType::UInt8;
            case 2:
                return is_signed ? DType::Int16 : DType::UInt16;
            case 4:
                return is_signed ? DType::Int32 : DType::UInt32;
            default:
                return is_signed ? DType::Int64 : DType::UInt64;
        }
    }
    throw ReadError("HDF5: unsupported datatype class");
}

/**
 * @brief Reads a full HDF5 dataset into a freshly-allocated, owning `NDArray`.
 *
 * The output shape and dtype are taken from the file. A dataset whose
 * datatype is an `ARRAY` of a scalar type — h5py's "(n,) of k-tuples" trick,
 * used e.g. by MED's `H5M` node coordinates — is unpacked into a plain
 * `(n, k)` `NDArray` by appending the array dimensions to the dataset's
 * shape, rather than exposed as a compound/array-typed element.
 * A scalar (0-dimensional) dataset comes back with shape `{1}`.
 *
 * @param loc Group or file handle the dataset lives under.
 * @param rName Name of the dataset to read.
 * @return A new owning `NDArray` holding the dataset's contents.
 * @throws ReadError if the dataset is missing or the read fails.
 */
inline NDArray read_dataset(hid_t loc, const std::string& rName) {
    Hid d(H5Dopen2(loc, rName.c_str(), H5P_DEFAULT), H5Dclose);
    if (!d.Valid())
        throw ReadError("HDF5: missing dataset '" + rName + "'");
    Hid space(H5Dget_space(d), H5Sclose);
    int ndim = H5Sget_simple_extent_ndims(space);
    std::vector<hsize_t> hdims(ndim > 0 ? ndim : 0);
    if (ndim > 0)
        H5Sget_simple_extent_dims(space, hdims.data(), nullptr);
    Hid dt(H5Dget_type(d), H5Tclose);

    std::vector<std::size_t> shape(hdims.begin(), hdims.end());
    if (shape.empty())
        shape.push_back(1);  // scalar -> length-1

    DType mdt;
    if (H5Tget_class(dt) == H5T_ARRAY) {
        Hid base(H5Tget_super(dt), H5Tclose);
        mdt = dtype_from_h5(base);
        int arank = H5Tget_array_ndims(dt);
        std::vector<hsize_t> adims(arank > 0 ? arank : 0);
        if (arank > 0)
            H5Tget_array_dims2(dt, adims.data());
        for (hsize_t ad : adims)
            shape.push_back(static_cast<std::size_t>(ad));
    } else {
        mdt = dtype_from_h5(dt);
    }

    NDArray out(mdt, shape);
    if (out.Size() > 0) {
        // For ARRAY-typed datasets the memory type must be the matching array
        // type; for scalar types the plain native type suffices.
        if (H5Tget_class(dt) == H5T_ARRAY) {
            int arank = H5Tget_array_ndims(dt);
            std::vector<hsize_t> adims(arank > 0 ? arank : 0);
            if (arank > 0)
                H5Tget_array_dims2(dt, adims.data());
            Hid mem(H5Tarray_create2(native_type(mdt), arank, adims.data()), H5Tclose);
            if (H5Dread(d, mem, H5S_ALL, H5S_ALL, H5P_DEFAULT, out.Data()) < 0)
                throw ReadError("HDF5: failed reading dataset '" + rName + "'");
        } else if (H5Dread(d, native_type(mdt), H5S_ALL, H5S_ALL, H5P_DEFAULT, out.Data()) < 0) {
            throw ReadError("HDF5: failed reading dataset '" + rName + "'");
        }
    }
    return out;
}

/**
 * @brief Writes a full dataset in one call, optionally gzip-compressed.
 *
 * When `gzip_level >= 0` and `arr` is non-empty, the dataset is created
 * chunked with a single chunk spanning the whole shape and gzip deflate
 * filtering enabled at that level; otherwise it is a plain contiguous
 * dataset. Uses `file_type(arr.Dtype())` for the on-disk type and
 * `native_type(arr.Dtype())` for the in-memory transfer type.
 *
 * @param loc Group or file handle to create the dataset under.
 * @param rName Name for the new dataset.
 * @param rArr Data to write; its shape and dtype determine the dataset's.
 * @param gzip_level gzip compression level (0-9), or negative to disable
 *                   compression (the default).
 * @throws WriteError if the dataset cannot be created or the write fails.
 */
inline void write_dataset(hid_t loc, const std::string& rName, const NDArray& rArr,
                          int gzip_level = -1) {
    std::vector<hsize_t> hdims(rArr.Shape().begin(), rArr.Shape().end());
    if (hdims.empty())
        hdims.push_back(0);
    Hid space(H5Screate_simple(static_cast<int>(hdims.size()), hdims.data(), nullptr), H5Sclose);

    Hid dcpl(H5Pcreate(H5P_DATASET_CREATE), H5Pclose);
    if (gzip_level >= 0 && rArr.Size() > 0) {
        H5Pset_chunk(dcpl, static_cast<int>(hdims.size()), hdims.data());
        H5Pset_deflate(dcpl, static_cast<unsigned>(gzip_level));
    }

    Hid d(H5Dcreate2(loc, rName.c_str(), file_type(rArr.Dtype()), space, H5P_DEFAULT, dcpl,
                     H5P_DEFAULT),
          H5Dclose);
    if (!d.Valid())
        throw WriteError("HDF5: could not create dataset '" + rName + "'");
    if (rArr.Size() > 0) {
        if (H5Dwrite(d, native_type(rArr.Dtype()), H5S_ALL, H5S_ALL, H5P_DEFAULT, rArr.Data()) < 0)
            throw WriteError("HDF5: failed writing dataset '" + rName + "'");
    }
}

// ---- attribute helpers ----

/**
 * @brief Whether an attribute named `name` exists on `loc`.
 * @param loc Object (group/dataset/file) to check.
 * @param rName Attribute name to test.
 * @return `true` if the attribute exists.
 */
inline bool has_attr(hid_t loc, const std::string& rName) {
    return H5Aexists(loc, rName.c_str()) > 0;
}

/**
 * @brief Reads a scalar integer attribute.
 * @param loc Object the attribute is attached to.
 * @param rName Attribute name.
 * @return The attribute's value as `int64_t`.
 * @throws ReadError if the attribute is missing or unreadable.
 */
inline std::int64_t read_attr_int(hid_t loc, const std::string& rName) {
    Hid a(H5Aopen(loc, rName.c_str(), H5P_DEFAULT), H5Aclose);
    if (!a.Valid())
        throw ReadError("HDF5: missing attribute '" + rName + "'");
    std::int64_t v = 0;
    if (H5Aread(a, H5T_NATIVE_INT64, &v) < 0)
        throw ReadError("HDF5: failed reading attribute '" + rName + "'");
    return v;
}

/**
 * @brief Writes a scalar integer attribute.
 * @param loc Object to attach the attribute to.
 * @param rName Attribute name.
 * @param v Value to write.
 * @param ftype On-disk integer type to store as (default `H5T_STD_I64LE`).
 * @throws WriteError if the attribute cannot be created.
 */
inline void write_attr_int(hid_t loc, const std::string& rName, std::int64_t v,
                           hid_t ftype = H5T_STD_I64LE) {
    Hid space(H5Screate(H5S_SCALAR), H5Sclose);
    Hid a(H5Acreate2(loc, rName.c_str(), ftype, space, H5P_DEFAULT, H5P_DEFAULT), H5Aclose);
    if (!a.Valid())
        throw WriteError("HDF5: could not create attribute '" + rName + "'");
    H5Awrite(a, H5T_NATIVE_INT64, &v);
}

/**
 * @brief Reads a string attribute, handling both variable- and fixed-length
 * HDF5 string encodings.
 *
 * For a fixed-length (`NULLPAD`) string, reads into a same-sized buffer
 * (converting to `NULLTERM` would truncate the last character to make room
 * for a terminator) and then trims trailing NUL bytes and spaces.
 * @param loc Object the attribute is attached to.
 * @param rName Attribute name.
 * @return The attribute's value as a `std::string`.
 * @throws ReadError if the attribute is missing or unreadable.
 */
inline std::string read_attr_string(hid_t loc, const std::string& rName) {
    Hid a(H5Aopen(loc, rName.c_str(), H5P_DEFAULT), H5Aclose);
    if (!a.Valid())
        throw ReadError("HDF5: missing attribute '" + rName + "'");
    Hid t(H5Aget_type(a), H5Tclose);
    if (H5Tis_variable_str(t) > 0) {
        char* p = nullptr;
        Hid mt(H5Tcopy(H5T_C_S1), H5Tclose);
        H5Tset_size(mt, H5T_VARIABLE);
        H5Tset_cset(mt, H5Tget_cset(t));
        if (H5Aread(a, mt, &p) < 0 || p == nullptr)
            throw ReadError("HDF5: failed reading attribute '" + rName + "'");
        std::string out(p);
        H5free_memory(p);
        return out;
    }
    std::size_t sz = H5Tget_size(t);
    std::vector<char> buf(sz + 1, '\0');
    Hid mt(H5Tcopy(H5T_C_S1), H5Tclose);
    H5Tset_size(mt, sz);
    H5Tset_cset(mt, H5Tget_cset(t));
    // NULLPAD memory type: converting a NULLPAD file string into a NULLTERM
    // memory string of the same size would truncate the last character to
    // make room for the terminator.
    H5Tset_strpad(mt, H5T_STR_NULLPAD);
    if (H5Aread(a, mt, buf.data()) < 0)
        throw ReadError("HDF5: failed reading attribute '" + rName + "'");
    // trim trailing NULs/spaces
    std::string out(buf.data(), strnlen(buf.data(), sz));
    while (!out.empty() && out.back() == ' ')
        out.pop_back();
    return out;
}

/**
 * @brief Writes a string attribute the way h5py does by default:
 * variable-length, UTF-8-tagged.
 *
 * Matching h5py's convention keeps files produced by the C++ writer
 * byte-for-byte compatible with the Python/h5py writer's output.
 * @param loc Object to attach the attribute to.
 * @param rName Attribute name.
 * @param rValue String value to write.
 * @throws WriteError if the attribute cannot be created.
 */
inline void write_attr_string(hid_t loc, const std::string& rName, const std::string& rValue) {
    Hid space(H5Screate(H5S_SCALAR), H5Sclose);
    Hid t(H5Tcopy(H5T_C_S1), H5Tclose);
    H5Tset_size(t, H5T_VARIABLE);
    H5Tset_cset(t, H5T_CSET_UTF8);
    Hid a(H5Acreate2(loc, rName.c_str(), t, space, H5P_DEFAULT, H5P_DEFAULT), H5Aclose);
    if (!a.Valid())
        throw WriteError("HDF5: could not create attribute '" + rName + "'");
    const char* p = rValue.c_str();
    H5Awrite(a, t, &p);
}

/**
 * @brief Lists the link (child) names directly under a group, in HDF5's
 * default name-index iteration order.
 * @param loc Group handle to list.
 * @return Child link names, in name order.
 */
inline std::vector<std::string> group_links(hid_t loc) {
    H5G_info_t info;
    H5Gget_info(loc, &info);
    std::vector<std::string> names;
    names.reserve(info.nlinks);
    for (hsize_t i = 0; i < info.nlinks; ++i) {
        ssize_t len =
            H5Lget_name_by_idx(loc, ".", H5_INDEX_NAME, H5_ITER_INC, i, nullptr, 0, H5P_DEFAULT);
        std::string name(static_cast<std::size_t>(len), '\0');
        H5Lget_name_by_idx(loc, ".", H5_INDEX_NAME, H5_ITER_INC, i, name.data(),
                           static_cast<std::size_t>(len) + 1, H5P_DEFAULT);
        names.push_back(std::move(name));
    }
    return names;
}

/**
 * @brief Like `group_links`, but iterates in HDF5 link *creation* order when
 * the group tracks it (matching h5py's iteration order on
 * `track_order=True` files); silently falls back to name order otherwise.
 *
 * Needed wherever the order children were created in is semantically
 * significant rather than incidental — e.g. MED's `MAI` cell-block groups,
 * whose order must line up with the corresponding entries in `cell_data`/
 * `cell_sets`, which are positional (not keyed by group name).
 * @param loc Group handle to list.
 * @return Child link names, in creation order if indexed, else name order.
 */
inline std::vector<std::string> group_links_crt(hid_t loc) {
    H5G_info_t info;
    H5Gget_info(loc, &info);
    std::vector<std::string> names;
    names.reserve(info.nlinks);
    for (hsize_t i = 0; i < info.nlinks; ++i) {
        ssize_t len = H5Lget_name_by_idx(loc, ".", H5_INDEX_CRT_ORDER, H5_ITER_INC, i, nullptr, 0,
                                         H5P_DEFAULT);
        if (len < 0)
            return group_links(loc);  // creation order not indexed
        std::string name(static_cast<std::size_t>(len), '\0');
        H5Lget_name_by_idx(loc, ".", H5_INDEX_CRT_ORDER, H5_ITER_INC, i, name.data(),
                           static_cast<std::size_t>(len) + 1, H5P_DEFAULT);
        names.push_back(std::move(name));
    }
    return names;
}

/**
 * @brief RAII guard that silences HDF5's default stderr error-stack printing
 * for its lifetime, restoring the previous handler on destruction.
 *
 * The library's own error reporting is redundant here since every failure
 * this codebase cares about is converted to a `ReadError`/`WriteError`
 * exception; without this guard, HDF5 would additionally dump a raw error
 * stack to stderr on every recoverable failure (e.g. a probing "does this
 * attribute exist" call that's expected to fail sometimes).
 */
struct SilenceErrors {
    SilenceErrors() {
        H5Eget_auto2(H5E_DEFAULT, &mOldFunc, &mOldData);
        H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
    }
    ~SilenceErrors() { H5Eset_auto2(H5E_DEFAULT, mOldFunc, mOldData); }
    H5E_auto2_t mOldFunc = nullptr;
    void* mOldData = nullptr;
};

}  // namespace h5
}  // namespace meshioplusplus

#endif  // MESHIOPLUSPLUS_HAS_HDF5
