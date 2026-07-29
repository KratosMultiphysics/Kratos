// KRATOS  __  __          _    _                _ _           _   _
//        |  \/  | ___  __| |  / \   _ __  _ __ | (_) ___ __ _| |_(_) ___  _ ___
//        | |\/| |/ _ \/ _` | / _ \ | '_ \| '_ \| | |/ __/ _` | __| |/ _ \| '_  |
//        | |  | |  __/ (_| |/ ___ \| |_) | |_) | | | (_| (_| | |_| | (_) | | | |
//        |_|  |_|\___|\__,_/_/   \_\ .__/| .__/|_|_|\___\__,_|\__|_|\___/|_| |_|
//                                  |_|   |_|
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

// System includes
#include <algorithm>
#include <sstream>

// External includes

// Project includes
#include "med_hdf5_utilities.h"
#include "includes/define.h"
#include "includes/exception.h"

namespace Kratos {
namespace MedHdf5Utilities {

namespace {

std::string TrimNullTerminated(
    const char* pData,
    const std::size_t Size)
{
    std::size_t len = 0;
    while (len < Size && pData[len] != '\0') ++len;
    return std::string(pData, len);
}

} // anonymous namespace

bool LinkExists(
    const hid_t FileId,
    const std::string& rPath)
{
    KRATOS_TRY

    std::stringstream path_stream(rPath);
    std::string segment;
    std::string current_path;

    while (std::getline(path_stream, segment, '/')) {
        if (segment.empty()) continue;
        current_path += "/" + segment;

        const htri_t exists = H5Lexists(FileId, current_path.c_str(), H5P_DEFAULT);
        KRATOS_ERROR_IF(exists < 0) << "H5Lexists failed for path \"" << current_path << "\"!" << std::endl;
        if (!exists) return false;
    }

    return true;

    KRATOS_CATCH("")
}

std::vector<std::string> ListGroupChildren(
    const hid_t FileId,
    const std::string& rGroupPath)
{
    KRATOS_TRY

    const hid_t group_id = H5Gopen2(FileId, rGroupPath.c_str(), H5P_DEFAULT);
    KRATOS_ERROR_IF(group_id < 0) << "Could not open group \"" << rGroupPath << "\"!" << std::endl;

    H5G_info_t group_info;
    KRATOS_ERROR_IF(H5Gget_info(group_id, &group_info) < 0) << "H5Gget_info failed for group \"" << rGroupPath << "\"!" << std::endl;

    std::vector<std::string> names(group_info.nlinks);

    constexpr std::size_t buffer_size = 256;
    std::vector<char> buffer(buffer_size);

    for (hsize_t i = 0; i < group_info.nlinks; ++i) {
        const ssize_t name_size = H5Lget_name_by_idx(
            FileId, rGroupPath.c_str(), H5_INDEX_NAME, H5_ITER_INC, i,
            buffer.data(), buffer_size, H5P_DEFAULT);
        KRATOS_ERROR_IF(name_size < 0) << "H5Lget_name_by_idx failed for group \"" << rGroupPath << "\"!" << std::endl;
        KRATOS_ERROR_IF(static_cast<std::size_t>(name_size) >= buffer_size)
            << "Link name in group \"" << rGroupPath << "\" exceeds the buffer size!" << std::endl;

        names[i] = std::string(buffer.data(), name_size);
    }

    KRATOS_ERROR_IF(H5Gclose(group_id) < 0) << "H5Gclose failed for group \"" << rGroupPath << "\"!" << std::endl;

    return names;

    KRATOS_CATCH("")
}

hid_t CreateGroup(
    const hid_t FileId,
    const std::string& rPath)
{
    KRATOS_TRY

    const hid_t lcpl_id = H5Pcreate(H5P_LINK_CREATE);
    KRATOS_ERROR_IF(H5Pset_create_intermediate_group(lcpl_id, 1) < 0) << "H5Pset_create_intermediate_group failed!" << std::endl;

    const hid_t group_id = H5Gcreate2(FileId, rPath.c_str(), lcpl_id, H5P_DEFAULT, H5P_DEFAULT);
    KRATOS_ERROR_IF(group_id < 0) << "Could not create group \"" << rPath << "\"!" << std::endl;

    KRATOS_ERROR_IF(H5Pclose(lcpl_id) < 0) << "H5Pclose failed!" << std::endl;

    return group_id;

    KRATOS_CATCH("")
}

void CreateAndCloseGroup(
    const hid_t FileId,
    const std::string& rPath)
{
    KRATOS_TRY

    const hid_t group_id = CreateGroup(FileId, rPath);
    KRATOS_ERROR_IF(H5Gclose(group_id) < 0) << "H5Gclose failed for group \"" << rPath << "\"!" << std::endl;

    KRATOS_CATCH("")
}

int64_t ReadInt64Attribute(
    const hid_t FileId,
    const std::string& rObjectPath,
    const std::string& rName)
{
    KRATOS_TRY

    const hid_t attr_id = H5Aopen_by_name(FileId, rObjectPath.c_str(), rName.c_str(), H5P_DEFAULT, H5P_DEFAULT);
    KRATOS_ERROR_IF(attr_id < 0) << "Could not open attribute \"" << rName << "\" on \"" << rObjectPath << "\"!" << std::endl;

    int64_t value{};
    const herr_t err = H5Aread(attr_id, H5T_NATIVE_INT64, &value);
    KRATOS_ERROR_IF(H5Aclose(attr_id) < 0) << "H5Aclose failed!" << std::endl;
    KRATOS_ERROR_IF(err < 0) << "Could not read attribute \"" << rName << "\" on \"" << rObjectPath << "\"!" << std::endl;

    return value;

    KRATOS_CATCH("")
}

void WriteInt64Attribute(
    const hid_t FileId,
    const std::string& rObjectPath,
    const std::string& rName,
    const int64_t Value)
{
    KRATOS_TRY

    const hid_t space_id = H5Screate(H5S_SCALAR);
    const hid_t attr_id = H5Acreate_by_name(
        FileId, rObjectPath.c_str(), rName.c_str(), H5T_STD_I64LE, space_id,
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    KRATOS_ERROR_IF(attr_id < 0) << "Could not create attribute \"" << rName << "\" on \"" << rObjectPath << "\"!" << std::endl;

    KRATOS_ERROR_IF(H5Awrite(attr_id, H5T_NATIVE_INT64, &Value) < 0)
        << "Could not write attribute \"" << rName << "\" on \"" << rObjectPath << "\"!" << std::endl;

    KRATOS_ERROR_IF(H5Aclose(attr_id) < 0) << "H5Aclose failed!" << std::endl;
    KRATOS_ERROR_IF(H5Sclose(space_id) < 0) << "H5Sclose failed!" << std::endl;

    KRATOS_CATCH("")
}

void WriteFloat64Attribute(
    const hid_t FileId,
    const std::string& rObjectPath,
    const std::string& rName,
    const double Value)
{
    KRATOS_TRY

    const hid_t space_id = H5Screate(H5S_SCALAR);
    const hid_t attr_id = H5Acreate_by_name(
        FileId, rObjectPath.c_str(), rName.c_str(), H5T_IEEE_F64LE, space_id,
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    KRATOS_ERROR_IF(attr_id < 0) << "Could not create attribute \"" << rName << "\" on \"" << rObjectPath << "\"!" << std::endl;

    KRATOS_ERROR_IF(H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &Value) < 0)
        << "Could not write attribute \"" << rName << "\" on \"" << rObjectPath << "\"!" << std::endl;

    KRATOS_ERROR_IF(H5Aclose(attr_id) < 0) << "H5Aclose failed!" << std::endl;
    KRATOS_ERROR_IF(H5Sclose(space_id) < 0) << "H5Sclose failed!" << std::endl;

    KRATOS_CATCH("")
}

void WriteStringAttribute(
    const hid_t FileId,
    const std::string& rObjectPath,
    const std::string& rName,
    const std::string& rValue)
{
    KRATOS_TRY

    const hid_t type_id = H5Tcopy(H5T_C_S1);
    KRATOS_ERROR_IF(H5Tset_size(type_id, rValue.size() + 1) < 0) << "H5Tset_size failed!" << std::endl;
    KRATOS_ERROR_IF(H5Tset_strpad(type_id, H5T_STR_NULLTERM) < 0) << "H5Tset_strpad failed!" << std::endl;

    const hid_t space_id = H5Screate(H5S_SCALAR);
    const hid_t attr_id = H5Acreate_by_name(
        FileId, rObjectPath.c_str(), rName.c_str(), type_id, space_id,
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    KRATOS_ERROR_IF(attr_id < 0) << "Could not create attribute \"" << rName << "\" on \"" << rObjectPath << "\"!" << std::endl;

    KRATOS_ERROR_IF(H5Awrite(attr_id, type_id, rValue.c_str()) < 0)
        << "Could not write attribute \"" << rName << "\" on \"" << rObjectPath << "\"!" << std::endl;

    KRATOS_ERROR_IF(H5Aclose(attr_id) < 0) << "H5Aclose failed!" << std::endl;
    KRATOS_ERROR_IF(H5Sclose(space_id) < 0) << "H5Sclose failed!" << std::endl;
    KRATOS_ERROR_IF(H5Tclose(type_id) < 0) << "H5Tclose failed!" << std::endl;

    KRATOS_CATCH("")
}

std::vector<int64_t> ReadInt64Dataset(
    const hid_t FileId,
    const std::string& rPath)
{
    KRATOS_TRY

    const hid_t dset_id = H5Dopen2(FileId, rPath.c_str(), H5P_DEFAULT);
    KRATOS_ERROR_IF(dset_id < 0) << "Could not open dataset \"" << rPath << "\"!" << std::endl;

    const hid_t space_id = H5Dget_space(dset_id);
    const hssize_t num_points = H5Sget_simple_extent_npoints(space_id);
    KRATOS_ERROR_IF(num_points < 0) << "H5Sget_simple_extent_npoints failed for dataset \"" << rPath << "\"!" << std::endl;

    std::vector<int64_t> data(static_cast<std::size_t>(num_points));
    if (num_points > 0) {
        KRATOS_ERROR_IF(H5Dread(dset_id, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data()) < 0)
            << "Could not read dataset \"" << rPath << "\"!" << std::endl;
    }

    KRATOS_ERROR_IF(H5Sclose(space_id) < 0) << "H5Sclose failed!" << std::endl;
    KRATOS_ERROR_IF(H5Dclose(dset_id) < 0) << "H5Dclose failed!" << std::endl;

    return data;

    KRATOS_CATCH("")
}

std::vector<double> ReadFloat64Dataset(
    const hid_t FileId,
    const std::string& rPath)
{
    KRATOS_TRY

    const hid_t dset_id = H5Dopen2(FileId, rPath.c_str(), H5P_DEFAULT);
    KRATOS_ERROR_IF(dset_id < 0) << "Could not open dataset \"" << rPath << "\"!" << std::endl;

    const hid_t space_id = H5Dget_space(dset_id);
    const hssize_t num_points = H5Sget_simple_extent_npoints(space_id);
    KRATOS_ERROR_IF(num_points < 0) << "H5Sget_simple_extent_npoints failed for dataset \"" << rPath << "\"!" << std::endl;

    std::vector<double> data(static_cast<std::size_t>(num_points));
    if (num_points > 0) {
        KRATOS_ERROR_IF(H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data()) < 0)
            << "Could not read dataset \"" << rPath << "\"!" << std::endl;
    }

    KRATOS_ERROR_IF(H5Sclose(space_id) < 0) << "H5Sclose failed!" << std::endl;
    KRATOS_ERROR_IF(H5Dclose(dset_id) < 0) << "H5Dclose failed!" << std::endl;

    return data;

    KRATOS_CATCH("")
}

void WriteInt64Dataset(
    const hid_t FileId,
    const std::string& rPath,
    const std::vector<int64_t>& rData)
{
    KRATOS_TRY

    // Mirrors the observed behavior of libmed's "MEDmeshEntityFamilyNumberWr" /
    // "MEDmeshGlobalNumberWr": writing zero entities is a no-op (no dataset created).
    if (rData.empty()) return;

    const hsize_t dims[1] = { rData.size() };
    const hid_t space_id = H5Screate_simple(1, dims, nullptr);

    const hid_t lcpl_id = H5Pcreate(H5P_LINK_CREATE);
    KRATOS_ERROR_IF(H5Pset_create_intermediate_group(lcpl_id, 1) < 0) << "H5Pset_create_intermediate_group failed!" << std::endl;

    const hid_t dset_id = H5Dcreate2(FileId, rPath.c_str(), H5T_STD_I64LE, space_id, lcpl_id, H5P_DEFAULT, H5P_DEFAULT);
    KRATOS_ERROR_IF(dset_id < 0) << "Could not create dataset \"" << rPath << "\"!" << std::endl;

    KRATOS_ERROR_IF(H5Dwrite(dset_id, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, rData.data()) < 0)
        << "Could not write dataset \"" << rPath << "\"!" << std::endl;

    KRATOS_ERROR_IF(H5Dclose(dset_id) < 0) << "H5Dclose failed!" << std::endl;
    KRATOS_ERROR_IF(H5Sclose(space_id) < 0) << "H5Sclose failed!" << std::endl;
    KRATOS_ERROR_IF(H5Pclose(lcpl_id) < 0) << "H5Pclose failed!" << std::endl;

    KRATOS_CATCH("")
}

void WriteFloat64Dataset(
    const hid_t FileId,
    const std::string& rPath,
    const std::vector<double>& rData)
{
    KRATOS_TRY

    const hsize_t dims[1] = { rData.size() };
    const hid_t space_id = H5Screate_simple(1, dims, nullptr);

    const hid_t lcpl_id = H5Pcreate(H5P_LINK_CREATE);
    KRATOS_ERROR_IF(H5Pset_create_intermediate_group(lcpl_id, 1) < 0) << "H5Pset_create_intermediate_group failed!" << std::endl;

    const hid_t dset_id = H5Dcreate2(FileId, rPath.c_str(), H5T_IEEE_F64LE, space_id, lcpl_id, H5P_DEFAULT, H5P_DEFAULT);
    KRATOS_ERROR_IF(dset_id < 0) << "Could not create dataset \"" << rPath << "\"!" << std::endl;

    if (!rData.empty()) {
        KRATOS_ERROR_IF(H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rData.data()) < 0)
            << "Could not write dataset \"" << rPath << "\"!" << std::endl;
    }

    KRATOS_ERROR_IF(H5Dclose(dset_id) < 0) << "H5Dclose failed!" << std::endl;
    KRATOS_ERROR_IF(H5Sclose(space_id) < 0) << "H5Sclose failed!" << std::endl;
    KRATOS_ERROR_IF(H5Pclose(lcpl_id) < 0) << "H5Pclose failed!" << std::endl;

    KRATOS_CATCH("")
}

std::vector<std::string> ReadFixedWidthStringDataset(
    const hid_t FileId,
    const std::string& rPath,
    const std::size_t Width)
{
    KRATOS_TRY

    const hid_t dset_id = H5Dopen2(FileId, rPath.c_str(), H5P_DEFAULT);
    KRATOS_ERROR_IF(dset_id < 0) << "Could not open dataset \"" << rPath << "\"!" << std::endl;

    const hid_t space_id = H5Dget_space(dset_id);
    const hssize_t num_points = H5Sget_simple_extent_npoints(space_id);
    KRATOS_ERROR_IF(num_points < 0) << "H5Sget_simple_extent_npoints failed for dataset \"" << rPath << "\"!" << std::endl;

    const hsize_t array_dims[1] = { Width };
    const hid_t array_type_id = H5Tarray_create2(H5T_NATIVE_CHAR, 1, array_dims);

    std::vector<char> buffer(static_cast<std::size_t>(num_points) * Width);
    if (num_points > 0) {
        KRATOS_ERROR_IF(H5Dread(dset_id, array_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data()) < 0)
            << "Could not read dataset \"" << rPath << "\"!" << std::endl;
    }

    KRATOS_ERROR_IF(H5Tclose(array_type_id) < 0) << "H5Tclose failed!" << std::endl;
    KRATOS_ERROR_IF(H5Sclose(space_id) < 0) << "H5Sclose failed!" << std::endl;
    KRATOS_ERROR_IF(H5Dclose(dset_id) < 0) << "H5Dclose failed!" << std::endl;

    std::vector<std::string> names(static_cast<std::size_t>(num_points));
    for (std::size_t i = 0; i < names.size(); ++i) {
        names[i] = TrimNullTerminated(&buffer[i * Width], Width);
    }

    return names;

    KRATOS_CATCH("")
}

void WriteFixedWidthStringDataset(
    const hid_t FileId,
    const std::string& rPath,
    const std::vector<std::string>& rNames,
    const std::size_t Width)
{
    KRATOS_TRY

    const hsize_t array_dims[1] = { Width };
    const hid_t array_type_id = H5Tarray_create2(H5T_NATIVE_CHAR, 1, array_dims);

    std::vector<char> buffer(rNames.size() * Width, '\0');
    for (std::size_t i = 0; i < rNames.size(); ++i) {
        const auto& r_name = rNames[i];
        KRATOS_ERROR_IF(r_name.size() >= Width)
            << "Name \"" << r_name << "\" exceeds the fixed-width limit of " << Width << " characters!" << std::endl;
        std::copy(r_name.begin(), r_name.end(), buffer.begin() + i * Width);
    }

    const hsize_t dims[1] = { rNames.size() };
    const hid_t space_id = H5Screate_simple(1, dims, nullptr);

    const hid_t lcpl_id = H5Pcreate(H5P_LINK_CREATE);
    KRATOS_ERROR_IF(H5Pset_create_intermediate_group(lcpl_id, 1) < 0) << "H5Pset_create_intermediate_group failed!" << std::endl;

    const hid_t dset_id = H5Dcreate2(FileId, rPath.c_str(), array_type_id, space_id, lcpl_id, H5P_DEFAULT, H5P_DEFAULT);
    KRATOS_ERROR_IF(dset_id < 0) << "Could not create dataset \"" << rPath << "\"!" << std::endl;

    if (!rNames.empty()) {
        KRATOS_ERROR_IF(H5Dwrite(dset_id, array_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data()) < 0)
            << "Could not write dataset \"" << rPath << "\"!" << std::endl;
    }

    KRATOS_ERROR_IF(H5Dclose(dset_id) < 0) << "H5Dclose failed!" << std::endl;
    KRATOS_ERROR_IF(H5Sclose(space_id) < 0) << "H5Sclose failed!" << std::endl;
    KRATOS_ERROR_IF(H5Pclose(lcpl_id) < 0) << "H5Pclose failed!" << std::endl;
    KRATOS_ERROR_IF(H5Tclose(array_type_id) < 0) << "H5Tclose failed!" << std::endl;

    KRATOS_CATCH("")
}

} // namespace MedHdf5Utilities
} // namespace Kratos
