//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  license: HDF5Application/license.txt
//
//  Main author:    Michael Andre, https://github.com/msandre
//                  Suneth Warnakulasuriya
//

// System includes
#include <algorithm>
#include <array>
#include <numeric>
#include <regex>
#include <sstream>
#include <utility>

// External includes

// Project includes
#include "includes/kratos_parameters.h"
#include "input_output/logger.h"
#include "utilities/builtin_timer.h"
#include "utilities/data_type_traits.h"
#include "utilities/string_utilities.h"

// Application includes
#include "custom_utilities/h5_data_type_traits.h"

// Include base h
#include "hdf5_file.h"

namespace Kratos
{
namespace HDF5
{
namespace Internals
{

bool IsPath(const std::string& rPath)
{
    return regex_match(rPath, std::regex("(/[\\w\\(\\)]+)+"));
}

std::vector<hsize_t> GetDataDimensions(
    const File& rFile,
    const std::string& rPath)
{
    const std::vector<unsigned> dims = rFile.GetDataDimensions(rPath);
    std::vector<hsize_t> h5_dims;
    h5_dims.reserve(dims.size());
    std::for_each(dims.begin(), dims.end(),
                  [&h5_dims](unsigned d) { h5_dims.push_back(d); });
    return h5_dims;
}
} // namespace Internals.

File::File(
    const DataCommunicator& rDataCommunicator,
    Parameters Settings)
    : mpDataCommunicator(&rDataCommunicator)
{
    KRATOS_TRY;

    Parameters default_params(R"(
            {
                "file_name" : "PLEASE_SPECIFY_HDF5_FILENAME",
                "file_access_mode": "exclusive",
                "file_driver": "sec2",
                "echo_level" : 0
            })");

    Settings.RecursivelyAddMissingParameters(default_params);

    m_file_name = Settings["file_name"].GetString();
    KRATOS_ERROR_IF(m_file_name == "PLEASE_SPECIFY_HDF5_FILENAME") << "Invalid file name: " << m_file_name << std::endl;

    hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    std::string file_driver = Settings["file_driver"].GetString();
    SetFileDriver(file_driver, fapl_id);

    std::string file_access_mode = Settings["file_access_mode"].GetString();
    if (file_access_mode == "exclusive") {
        m_file_id = H5Fcreate(m_file_name.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, fapl_id);
    } else if (file_access_mode == "truncate") {
        m_file_id = H5Fcreate(m_file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
    } else {
        // Save old error handler
        herr_t (*old_func)(hid_t, void*);
        void *old_client_data;

        H5Eget_auto(H5E_DEFAULT, &old_func, &old_client_data);

        // Turn off error handling
        H5Eset_auto(H5E_DEFAULT, NULL, NULL);

        // follwoing will fail if the hdf5 file is not found, and will print a failure msg to the output.
        // if "file_access_mode" is "read_write", then it is ok to fail the following call because,
        // if it is failed then the file will be created in the subsequent section.
        htri_t is_hdf5 = H5Fis_hdf5(m_file_name.c_str());

        // Restore previous error handler
        H5Eset_auto(H5E_DEFAULT, old_func, old_client_data);

        if (file_access_mode == "read_only") {
            KRATOS_ERROR_IF(is_hdf5 <= 0 && file_driver != "core") << "Invalid HDF5 file: " << m_file_name << std::endl;
            m_file_id = H5Fopen(m_file_name.c_str(), H5F_ACC_RDONLY, fapl_id);
        } else if (file_access_mode == "read_write") {
            if (is_hdf5 <= 0) {
                // creates the hdf5 file if the file is not found
                m_file_id = H5Fcreate(m_file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
            } else {
                // open the existing hdf5 file if file is found
                m_file_id = H5Fopen(m_file_name.c_str(), H5F_ACC_RDWR, fapl_id);
            }
        } else {
            KRATOS_ERROR << "Invalid \"file_access_mode\": " << file_access_mode
                         << std::endl;
        }
    }

    KRATOS_ERROR_IF(m_file_id < 0) << "Failed to open file: " << m_file_name << std::endl;

    KRATOS_ERROR_IF(H5Pclose(fapl_id) < 0) << "H5Pclose failed." << std::endl;

    m_echo_level = Settings["echo_level"].GetInt();

    KRATOS_CATCH("");
}

File::File(File&& rOther)
{
    mpDataCommunicator = rOther.mpDataCommunicator;
    m_file_name = std::move(rOther.m_file_name);
    m_file_id = rOther.m_file_id;
    rOther.m_file_id = -1;
    m_echo_level = rOther.m_echo_level;
    rOther.m_echo_level = 0;
}

File& File::operator=(File&& rOther)
{
    mpDataCommunicator = rOther.mpDataCommunicator;
    m_file_name = std::move(rOther.m_file_name);
    m_file_id = rOther.m_file_id;
    rOther.m_file_id = -1;
    m_echo_level = rOther.m_echo_level;
    rOther.m_echo_level = 0;
    return *this;
}

File::~File()
{
    if (0 <= m_file_id) {
        H5Fclose(m_file_id);
    }
}

bool File::HasPath(const std::string& rPath) const
{
    KRATOS_TRY;

    // Expects a valid path.
    KRATOS_ERROR_IF_NOT(Internals::IsPath(rPath))
        << "Invalid path: \"" << rPath
        << "\". Path should start with \"/\" and should only have characters A-Z, a-z, 0-9, \"/\", and \"_\"."
        << std::endl;

    std::vector<std::string> splitted_path = StringUtilities::SplitStringByDelimiter(rPath, '/');
    splitted_path.erase(std::remove_if(splitted_path.begin(), splitted_path.end(), [](const std::string& s) {return (s.size() == 0);}));
    std::string sub_path;

    for (const auto& r_link: splitted_path) {
        sub_path += '/' + r_link;

        htri_t link_found = H5Lexists(m_file_id, sub_path.c_str(), H5P_DEFAULT);
        KRATOS_ERROR_IF(link_found < 0) << "H5Lexists failed." << std::endl;
        if (!link_found) {
            return false;
        }

        htri_t object_found = H5Oexists_by_name(m_file_id, sub_path.c_str(), H5P_DEFAULT);
        KRATOS_ERROR_IF(object_found < 0) << "H5Oexists_by_name failed." << std::endl;
        if (!object_found) {
            return false;
        }
    }

    return true;

    KRATOS_CATCH("");
}

bool File::IsGroup(const std::string& rPath) const
{
    KRATOS_TRY;

    if (HasPath(rPath) == false) {// Expects a valid path.
        return false;
    }

    H5O_info_t object_info;
    KRATOS_ERROR_IF(H5Oget_info_by_name(m_file_id, rPath.c_str(), &object_info, H5P_DEFAULT) < 0)
        << "H5Oget_info_by_name failed." << std::endl;

    return (object_info.type == H5O_TYPE_GROUP);

    KRATOS_CATCH("");
}

bool File::IsDataSet(const std::string& rPath) const
{
    KRATOS_TRY;

    if (HasPath(rPath) == false) {// Expects a valid path.
        return false;
    }

    H5O_info_t object_info;
    KRATOS_ERROR_IF(H5Oget_info_by_name(m_file_id, rPath.c_str(), &object_info, H5P_DEFAULT) < 0)
        << "H5Oget_info_by_name failed." << std::endl;

    return (object_info.type == H5O_TYPE_DATASET);

    KRATOS_CATCH("");
}

bool File::HasAttribute(
    const std::string& rObjectPath,
    const std::string& rName) const
{
    KRATOS_TRY;

    htri_t status = H5Aexists_by_name(m_file_id, rObjectPath.c_str(), rName.c_str(), H5P_DEFAULT);
    KRATOS_ERROR_IF(status < 0) << "H5Aexists_by_name failed." << std::endl;
    return (status > 0);

    KRATOS_CATCH("");
}

void File::DeleteAttribute(
    const std::string& rObjectPath,
    const std::string& rName)
{
    KRATOS_TRY;

    KRATOS_ERROR_IF(H5Adelete_by_name(m_file_id, rObjectPath.c_str(), rName.c_str(), H5P_DEFAULT) < 0)
        << "H5Adelete_by_name failed.";

    KRATOS_CATCH("");
}

std::vector<std::string> File::GetAttributeNames(const std::string& rObjectPath) const
{
    KRATOS_TRY;

    constexpr unsigned max_ssize = 100;
    char buffer[max_ssize];
    // Get number of attributes.
    hid_t object_id = H5Oopen(m_file_id, rObjectPath.c_str(), H5P_DEFAULT);
    KRATOS_ERROR_IF(object_id < 0) << "H5Oopen failed." << std::endl;
    H5O_info_t object_info;
    KRATOS_ERROR_IF(H5Oget_info(object_id, &object_info) < 0)
        << "H5Oget_info failed." << std::endl;
    hsize_t num_attrs = object_info.num_attrs;
    std::vector<std::string> names(num_attrs);

    for (hsize_t i = 0; i < num_attrs; ++i) {
        // Get size of name.
        ssize_t ssize;
        ssize = H5Aget_name_by_idx(m_file_id, rObjectPath.c_str(), H5_INDEX_CRT_ORDER, H5_ITER_INC, i, buffer, max_ssize, H5P_DEFAULT);
        KRATOS_ERROR_IF(ssize < 0) << "H5Aget_name_by_idx failed." << std::endl;
        KRATOS_ERROR_IF(ssize > max_ssize)
            << "Attribute name size exceeds "
            << max_ssize << std::endl;

        names[i].resize(ssize);
        std::copy_n(buffer, ssize, names[i].begin());
    }

    KRATOS_ERROR_IF(H5Oclose(object_id) < 0) << "H5Oclose failed." << std::endl;
    return names;

    KRATOS_CATCH("");
}

void File::CreateGroup(const std::string& rPath)
{
    KRATOS_TRY;

    hid_t group_id = H5Gcreate(m_file_id,rPath.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    KRATOS_ERROR_IF(group_id < 0) << "H5Gcreate failed." << std::endl;
    KRATOS_ERROR_IF(H5Gclose(group_id) < 0) << "H5Gclose failed." << std::endl;

    KRATOS_CATCH("");
}

std::vector<std::string> File::GetLinkNames(const std::string& rGroupPath) const
{
    KRATOS_TRY;

    constexpr unsigned max_ssize = 100;
    char buffer[max_ssize];
    // Get number of links.
    hid_t group_id = H5Gopen(m_file_id, rGroupPath.c_str(), H5P_DEFAULT);
    KRATOS_ERROR_IF(group_id < 0) << "H5Gopen failed." << std::endl;

    H5G_info_t group_info;
    KRATOS_ERROR_IF(H5Gget_info(group_id, &group_info) < 0)
        << "H5Gget_info failed." << std::endl;
    hsize_t num_links = group_info.nlinks;
    std::vector<std::string> names(num_links);

    for (hsize_t i=0; i < num_links; ++i) {
        // Get size of name.
        ssize_t ssize;
        ssize = H5Lget_name_by_idx(m_file_id, rGroupPath.c_str(), H5_INDEX_NAME,
                                    H5_ITER_INC, i, buffer, max_ssize, H5P_DEFAULT);
        KRATOS_ERROR_IF(ssize < 0) << "H5Lget_name_by_idx failed." << std::endl;
        KRATOS_ERROR_IF(ssize > max_ssize) << "Link name size exceeds "
                                           << max_ssize << std::endl;
        names[i].resize(ssize);
        std::copy_n(buffer, ssize, names[i].begin());
    }
    KRATOS_ERROR_IF(H5Gclose(group_id) < 0) << "H5Gclose failed." << std::endl;
    return names;

    KRATOS_CATCH("");
}

std::vector<std::string> File::GetGroupNames(const std::string& rGroupPath) const
{
    KRATOS_TRY;

    std::vector<std::string> names;
    std::vector<std::string> link_names = GetLinkNames(rGroupPath);
    names.reserve(link_names.size());
    for (const auto& r_name : link_names) {
        if (IsGroup(rGroupPath + '/' + r_name)) {
            names.push_back(r_name);
        }
    }

    return names;

    KRATOS_CATCH("");
}

std::vector<std::string> File::GetDataSetNames(const std::string& rGroupPath) const
{
    KRATOS_TRY;

    std::vector<std::string> names;
    std::vector<std::string> link_names = GetLinkNames(rGroupPath);
    names.reserve(link_names.size());
    for (const auto& r_name : link_names) {
        if (IsDataSet(rGroupPath + '/' + r_name)) {
            names.push_back(r_name);
        }
    }
    return names;

    KRATOS_CATCH("");
}

void File::AddPath(const std::string& rPath)
{
    KRATOS_ERROR_IF_NOT(Internals::IsPath(rPath)) << "Invalid path: \"" << rPath << "\". Path should start with \"/\" and should only have characters A-Z, a-z, 0-9, \"/\", and \"_\"." << std::endl;

    std::vector<std::string> splitted_path = StringUtilities::SplitStringByDelimiter(rPath, '/');
    splitted_path.erase(std::remove_if(splitted_path.begin(), splitted_path.end(), [](const std::string& s) {return (s.size() == 0);}));
    std::string sub_path;
    for (const auto& r_link: splitted_path) {
        sub_path += '/' + r_link;
        if (!HasPath(sub_path)) {
            CreateGroup(sub_path); // Add missing link.
        } else {
            KRATOS_ERROR_IF_NOT(IsGroup(sub_path))
                << "Path exists and is not a group: " << sub_path << std::endl;
        }
    }
}

template<class TDataType>
void File::WriteAttribute(
    const std::string& rObjectPath,
    const std::string& rName,
    const TDataType& rValue)
{
    KRATOS_TRY;

    using type_traits = DataTypeTraits<TDataType>;

    static_assert(type_traits::IsContiguous, "Attributes needs to have data contiguous in memory.");

    constexpr auto local_dimension = type_traits::Dimension;

    BuiltinTimer timer;
    hid_t type_id, space_id, attr_id;

    if (HasAttribute(rObjectPath, rName)) {
        DeleteAttribute(rObjectPath, rName);
    }

    if constexpr(local_dimension == 0) {
        space_id = H5Screate(H5S_SCALAR);
    } else {
        std::vector<hsize_t> shape(local_dimension);
        type_traits::Shape(rValue, shape.data(), shape.data() + local_dimension);
        space_id = H5Screate_simple(local_dimension, shape.data(), nullptr);
    }

    KRATOS_ERROR_IF(space_id < 0) << "H5Screate failed." << std::endl;

    type_id = Internals::GetPrimitiveH5Type<TDataType>();
    attr_id = H5Acreate_by_name(m_file_id, rObjectPath.c_str(), rName.c_str(), type_id, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    KRATOS_ERROR_IF(attr_id < 0) << "H5Acreate_by_name failed." << std::endl;
    KRATOS_ERROR_IF(H5Awrite(attr_id, type_id, type_traits::GetContiguousData(rValue)) < 0) << "H5Awrite failed." << std::endl;
    KRATOS_ERROR_IF(H5Sclose(space_id) < 0) << "H5Sclose failed." << std::endl;
    KRATOS_ERROR_IF(H5Aclose(attr_id) < 0) << "H5Aclose failed." << std::endl;

    KRATOS_INFO_IF("HDF5Application", GetEchoLevel() == 2)
        << "Write time \"" << rObjectPath << '/' << rName
        << "\": " << timer.ElapsedSeconds() << std::endl;

    KRATOS_CATCH("Path: \"" + rObjectPath + '/' + rName + "\".");
}

template<class TDataType>
void File::WriteDataSet(
    const std::string& rPath,
    const TDataType& rData,
    WriteInfo& rInfo)
{
    WriteDataSetImpl<TDataType, DataTransferMode::Collective>(rPath, rData, rInfo);
}

template<class TDataType>
void File::WriteDataSetIndependent(
    const std::string& rPath,
    const TDataType& rData,
    WriteInfo& rInfo)
{
    WriteDataSetImpl<TDataType, DataTransferMode::Independent>(rPath, rData, rInfo);
}

std::vector<unsigned> File::GetDataDimensions(const std::string& rPath) const
{
    KRATOS_TRY;

    int ndims;
    hid_t dset_id, dspace_id;
    KRATOS_ERROR_IF((dset_id = H5Dopen(m_file_id, rPath.c_str(), H5P_DEFAULT)) < 0)
        << "H5Dopen failed." << std::endl;
    KRATOS_ERROR_IF((dspace_id = H5Dget_space(dset_id)) < 0)
        << "H5Dget_space failed." << std::endl;
    KRATOS_ERROR_IF((ndims = H5Sget_simple_extent_ndims(dspace_id)) < 0)
        << "H5Sget_simple_extent_ndims failed." << std::endl;

    std::vector<hsize_t> dims(ndims);
    KRATOS_ERROR_IF(H5Sget_simple_extent_dims(dspace_id, dims.data(), nullptr) < 0)
        << "H5Sget_simple_extent_dims failed" << std::endl;
    KRATOS_ERROR_IF(H5Sclose(dspace_id) < 0) << "H5Sclose failed." << std::endl;
    KRATOS_ERROR_IF(H5Dclose(dset_id) < 0) << "H5Dclose failed." << std::endl;

    return std::vector<unsigned>(dims.begin(), dims.end());

    KRATOS_CATCH("");
}

bool File::HasIntDataType(const std::string& rPath) const
{
    return HasDataType<int>(rPath);
}

bool File::HasFloatDataType(const std::string& rPath) const
{
    return HasDataType<double>(rPath);
}

void File::Flush()
{
    KRATOS_ERROR_IF(H5Fflush(m_file_id, H5F_SCOPE_GLOBAL) < 0)
        << "H5Fflush failed." << std::endl;
}

void File::Close()
{
    if (0 <= m_file_id) {
        const auto close_result = H5Fclose(m_file_id);
        KRATOS_ERROR_IF(close_result < 0) << "Failed to close " << m_file_name << " with error code " << close_result;
        m_file_id = -1;
    } else {
        KRATOS_WARNING("Invalid file handle") << "Attempt to close an invalid file" << std::endl;
    }
}

unsigned File::GetFileSize() const
{
    hsize_t size;
    KRATOS_ERROR_IF(H5Fget_filesize(m_file_id, &size) < 0)
        << "H5Fget_filesize failed." << std::endl;

    return size;
}

std::string File::GetFileName() const
{
    return m_file_name;
}

int File::GetEchoLevel() const
{
    return m_echo_level;
}

void File::SetEchoLevel(int Level)
{
    m_echo_level = Level;
}

const DataCommunicator& File::GetDataCommunicator() const
{
    return *mpDataCommunicator;
}

unsigned File::GetPID() const
{
    return mpDataCommunicator->Rank();
}

unsigned File::GetTotalProcesses() const
{
    return mpDataCommunicator->Size();
}

template<class TDataType>
void File::ReadAttribute(
    const std::string& rObjectPath,
    const std::string& rName,
    TDataType& rValue)
{
    KRATOS_TRY;

    using type_traits = DataTypeTraits<TDataType>;

    static_assert(type_traits::IsContiguous, "Attribute data should be contiguous in memory.");

    constexpr auto local_dimension = type_traits::Dimension;

    BuiltinTimer timer;

    hid_t mem_type_id, attr_type_id, space_id, attr_id;
    std::vector<hsize_t> shape(local_dimension);

    int ndims;

    mem_type_id = Internals::GetPrimitiveH5Type<TDataType>();

    attr_id = H5Aopen_by_name(m_file_id, rObjectPath.c_str(), rName.c_str(), H5P_DEFAULT, H5P_DEFAULT);
    KRATOS_ERROR_IF(attr_id < 0) << "H5Aopen_by_name failed." << std::endl;

    // Check data type.
    attr_type_id = H5Aget_type(attr_id);
    KRATOS_ERROR_IF(attr_type_id < 0) << "H5Aget_type failed." << std::endl;
    htri_t is_valid_type = H5Tequal(mem_type_id, attr_type_id);
    KRATOS_ERROR_IF(H5Tclose(attr_type_id) < 0) << "H5Tclose failed." << std::endl;
    KRATOS_ERROR_IF(is_valid_type < 0) << "H5Tequal failed." << std::endl;
    KRATOS_ERROR_IF(is_valid_type == 0) << "Memory and file data types are different." << std::endl;

    // Check dimensions.
    space_id = H5Aget_space(attr_id);
    KRATOS_ERROR_IF(space_id < 0)
        << "H5Aget_space failed." << std::endl;
    KRATOS_ERROR_IF((ndims = H5Sget_simple_extent_ndims(space_id)) < 0)
        << "H5Sget_simple_extent_ndims failed." << std::endl;
    KRATOS_ERROR_IF(ndims != local_dimension)
        << "Attribute \"" << rName << "\" has dimension mismatch [ memory dimension = "
        << local_dimension << ", file dimension = " << ndims  << " ].\n";
    KRATOS_ERROR_IF(H5Sget_simple_extent_dims(space_id, shape.data(), nullptr) < 0)
        << "H5Sget_simple_extent_dims failed" << std::endl;
    KRATOS_ERROR_IF(H5Sclose(space_id) < 0) << "H5Sclose failed." << std::endl;
    type_traits::Reshape(rValue, shape);

    // Read attribute.
    KRATOS_ERROR_IF(H5Aread(attr_id, mem_type_id, type_traits::GetContiguousData(rValue)) < 0)
        << "H5Aread failed." << std::endl;
    KRATOS_ERROR_IF(H5Aclose(attr_id) < 0)
        << "H5Aclose failed." << std::endl;

    KRATOS_INFO_IF("HDF5Application", GetEchoLevel() == 2)
        << "Read time \"" << rObjectPath << '/' << rName
        << "\": " << timer.ElapsedSeconds() << std::endl;

    KRATOS_CATCH("Path: \"" + rObjectPath + '/' + rName + "\".");
}

template<class TDataType>
void File::ReadDataSet(
    const std::string& rPath,
    TDataType& rData,
    const unsigned StartIndex,
    const unsigned BlockSize)
{
    ReadDataSetImpl<TDataType, DataTransferMode::Collective>(rPath, rData, StartIndex, BlockSize);
}

template<class TDataType>
void File::ReadDataSetIndependent(
    const std::string& rPath,
    TDataType& rData,
    const unsigned StartIndex,
    const unsigned BlockSize)
{
    ReadDataSetImpl<TDataType, DataTransferMode::Independent>(rPath, rData, StartIndex, BlockSize);
}

unsigned File::GetOpenObjectsCount() const
{
    KRATOS_TRY;
    ssize_t num_open_objects = H5Fget_obj_count(m_file_id, H5F_OBJ_ALL);
    KRATOS_ERROR_IF(num_open_objects < 0) << "H5Fget_obj_count failed." << std::endl;
    return num_open_objects;
    KRATOS_CATCH("");
}

hid_t File::GetFileId() const
{
    return m_file_id;
}

template<class TDataType>
bool File::HasDataType(const std::string& rPath) const
{
    KRATOS_TRY;

    hid_t dset_id, dtype_id;

    KRATOS_ERROR_IF((dset_id = H5Dopen(GetFileId(), rPath.c_str(), H5P_DEFAULT)) < 0)
        << "H5Dopen failed." << std::endl;
    KRATOS_ERROR_IF((dtype_id = H5Dget_type(dset_id)) < 0)
        << "H5Dget_type failed." << std::endl;
    H5T_class_t type = H5Tget_class(dtype_id);

    KRATOS_ERROR_IF(type == H5T_NO_CLASS) << "Invalid data type." << std::endl;
    KRATOS_ERROR_IF(H5Tclose(dtype_id) < 0) << "H5Tclose failed." << std::endl;
    KRATOS_ERROR_IF(H5Dclose(dset_id) < 0) << "H5Dclose failed." << std::endl;

    if constexpr(std::is_same_v<TDataType, int>) {
        return (type == H5T_INTEGER);
    } else if constexpr(std::is_same_v<TDataType, double>) {
        return (type == H5T_FLOAT);
    } else {
        static_assert(!std::is_same_v<TDataType, TDataType>, "Unsupported data type.");
    }

    KRATOS_CATCH("");
}

template<class TDataType>
void File::GetDataSet(
    hid_t& rDataSetId,
    hid_t& rDataSpaceId,
    const std::vector<hsize_t>& rDims,
    const std::string& rPath)
{
    KRATOS_TRY

    if (!HasPath(rPath)) {
        CreateNewDataSet(rDataSetId, rDataSpaceId, Internals::GetPrimitiveH5Type<TDataType>(), rDims, rPath);
    } else {
        KRATOS_ERROR_IF_NOT(HasDataType<TDataType>(rPath))
            << "Wrong scalar data type: " << rPath << std::endl;
        KRATOS_ERROR_IF(Internals::GetDataDimensions(*this, rPath) != rDims)
            << "Wrong dimensions: " << rPath << std::endl;

        rDataSetId = OpenExistingDataSet(rPath);
        KRATOS_ERROR_IF(rDataSetId < 0) << "H5Dopen failed." << std::endl;
        rDataSpaceId = H5Dget_space(rDataSetId);
    }

    KRATOS_CATCH("");
}

void File::CreateNewDataSet(
    hid_t& rDataSetId,
    hid_t& rDataSpaceId,
    const hid_t DataTypeId,
    const std::vector<hsize_t>& rDims,
    const std::string& rPath)
{
    rDataSpaceId = H5Screate_simple(rDims.size(), rDims.data(), nullptr);
    rDataSetId = H5Dcreate(GetFileId(), rPath.c_str(), DataTypeId, rDataSpaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    KRATOS_ERROR_IF(rDataSetId < 0) << "H5Dcreate failed." << std::endl;
}

hid_t File::OpenExistingDataSet(const std::string& rPath)
{
    const hid_t dset_id = H5Dopen(GetFileId(), rPath.c_str(), H5P_DEFAULT);
    KRATOS_ERROR_IF(dset_id < 0) << "H5Dopen failed." << std::endl;
    return dset_id;
}

void File::SetFileDriver(const std::string& rDriver, hid_t FaplId) const
{
    KRATOS_TRY;
#if (defined(_WIN32) || defined(_WIN64))
    KRATOS_ERROR_IF(rDriver != "windows")
        << "Unsupported (Windows) \"file_driver\": " << rDriver << std::endl;
    KRATOS_ERROR_IF(H5Pset_fapl_windows(FaplId) < 0)
        << "H5Pset_fapl_windows failed." << std::endl;
#else
    if (rDriver == "sec2")
    {
        KRATOS_ERROR_IF(H5Pset_fapl_sec2(FaplId) < 0)
            << "H5Pset_fapl_sec2 failed." << std::endl;
    }
    else if (rDriver == "stdio")
    {
        KRATOS_ERROR_IF(H5Pset_fapl_stdio(FaplId) < 0)
            << "H5Pset_fapl_stdio failed." << std::endl;
    }
    else if (rDriver == "core")
    {
        KRATOS_ERROR_IF(H5Pset_fapl_core(FaplId, 1000000, 0) < 0)
            << "H5Pset_fapl_core failed." << std::endl;
    }
    else if (rDriver == "mpio")
    {
#if defined(KRATOS_USING_MPI)
        KRATOS_ERROR_IF(H5Pset_fapl_mpio(FaplId, MPI_COMM_WORLD, MPI_INFO_NULL) < 0)
            << "H5Pset_fapl_mpio failed." << std::endl;
#else
        KRATOS_ERROR
            << "Kratos must be built with MPI for \"file_driver\"=\"mpio\"."
            << std::endl;
#endif
    }
    else
        KRATOS_ERROR << "Unsupported \"file_driver\": " << rDriver << std::endl;
#endif
    KRATOS_CATCH("");
}

template<class TDataType, File::DataTransferMode TDataTransferMode>
void File::WriteDataSetImpl(
    const std::string& rPath,
    const TDataType& rData,
    WriteInfo& rInfo)
{
    KRATOS_TRY

    using type_trait = DataTypeTraits<TDataType>;

    static_assert(type_trait::IsContiguous, "HDF5File can only write contiguous data sets.");

    constexpr auto local_dimension = type_trait::Dimension;

    constexpr auto global_dimension = (local_dimension == 1 ? 1 : 2);

    BuiltinTimer timer;

    // Create any missing subpaths.
    auto pos = rPath.find_last_of('/');
    if (pos != 0) {// Skip if last '/' is root.
        std::string sub_path = rPath.substr(0, pos);
        AddPath(sub_path);
    }

    // Initialize data space dimensions.
    std::vector<hsize_t> local_shape(local_dimension);
    type_trait::Shape(rData, local_shape.data(), local_shape.data() + local_dimension);

    const hsize_t number_of_local_primitive_data_values = type_trait::Size(rData);

    const auto& r_data_communicator = GetDataCommunicator();
    // get the maximized dimensions of the underlying data. Max is taken because,
    // there can be empty ranks which will give wrong sizes in the case of dynamic
    // data types.
    if constexpr(local_dimension >= 2) {
        if constexpr(type_trait::template IsDimensionDynamic<1>()) {
            // this is the matrix version. Hence it is better to get the max size
            // from all ranks.
            const auto max_size = r_data_communicator.MaxAll(local_shape[1]);

            // now check every non-empty ranks have the same sizes since this dimension
            // is a dynamic dimension.
            KRATOS_ERROR_IF(number_of_local_primitive_data_values > 0 && max_size != local_shape[1])
                << "Mismatching shapes found in different ranks. All ranks should have the same shapes in data sets.";

            local_shape[1] = max_size;
        }
    }

    // local_reduced_shape holds the max 2d flattened shape if the local_shape dimensions
    // are higher than 2.
    std::vector<hsize_t> global_shape(global_dimension, 0), local_reduced_shape(global_dimension, 0), local_shape_start(global_dimension, 0);

    // get total number of items to be written in the data set to the first dimension.
    global_shape[0] = r_data_communicator.SumAll(local_shape[0]);
    local_reduced_shape[0] = local_shape[0];

    if constexpr(global_dimension > 1) {
        // flattens higher dimensions into one since we write matrices which is the highest dimension
        // supported by paraview for visualization
        global_shape[1] = std::accumulate(local_shape.begin() + 1, local_shape.end(), hsize_t{1}, std::multiplies<hsize_t>());
        local_reduced_shape[1] = global_shape[1];
    }

    // Set the data type.
    hid_t dtype_id = Internals::GetPrimitiveH5Type<TDataType>();

    // Create and write the data set.
    hid_t dset_id{}, fspace_id{};
    GetDataSet<typename type_trait::PrimitiveType>(dset_id, fspace_id, global_shape, rPath);

    // here onwards the procedure differs shared memory and distributed memeory runs.
    // The same steps need to be applied if the HDF5Application is compiled with serial hdf5lib
    // and if the HDF5Application is compiled with mpi hdf5lib but running in shared memory parallelization
    // only.
    // Different steps has to be taken if it is run in a distributed memory environment.

    if (!r_data_communicator.IsDistributed()) {
        if (number_of_local_primitive_data_values > 0) {
            KRATOS_ERROR_IF(H5Dwrite(dset_id, dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, type_trait::GetContiguousData(rData)) < 0)
                << "H5Dwrite failed." << std::endl;
        } else {
            KRATOS_ERROR_IF(H5Dwrite(dset_id, dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, nullptr) < 0)
                << "H5Dwrite failed. Please ensure global data set is non-empty."
                << std::endl;
        }
    } else {
        #ifdef KRATOS_USING_MPI
            if constexpr(TDataTransferMode == DataTransferMode::Collective) {
                local_shape_start[0] = r_data_communicator.ScanSum(local_reduced_shape[0]) - local_reduced_shape[0];
            }

            if (TDataTransferMode == DataTransferMode::Collective || number_of_local_primitive_data_values > 0) {
                hid_t dxpl_id = H5Pcreate(H5P_DATASET_XFER);
                if constexpr(TDataTransferMode == DataTransferMode::Collective) {
                    H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE);
                } else {
                    H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_INDEPENDENT);
                }

                // select the local hyperslab
                H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, local_shape_start.data(), nullptr, local_reduced_shape.data(), nullptr);
                hid_t mspace_id = H5Screate_simple(global_dimension, local_reduced_shape.data(), nullptr);
                if (number_of_local_primitive_data_values > 0) {
                    KRATOS_ERROR_IF(H5Dwrite(dset_id, dtype_id, mspace_id, fspace_id, dxpl_id, type_trait::GetContiguousData(rData)) < 0)
                        << "H5Dwrite failed." << std::endl;
                } else {
                    KRATOS_ERROR_IF(H5Dwrite(dset_id, dtype_id, mspace_id, fspace_id, dxpl_id, nullptr) < 0)
                        << "H5Dwrite failed. Please ensure global data set is non-empty."
                        << std::endl;
                }

                KRATOS_ERROR_IF(H5Pclose(dxpl_id) < 0) << "H5Pclose failed." << std::endl;
                KRATOS_ERROR_IF(H5Sclose(mspace_id) < 0) << "H5Sclose failed." << std::endl;
            }
        #else
            KRATOS_ERROR << "HDFApplication is not compiled with MPI enabled";
        #endif
    }

    KRATOS_ERROR_IF(H5Sclose(fspace_id) < 0) << "H5Sclose failed." << std::endl;
    KRATOS_ERROR_IF(H5Dclose(dset_id) < 0) << "H5Dclose failed." << std::endl;

    // Set the write info.
    rInfo.StartIndex = local_shape_start[0];
    rInfo.BlockSize = local_reduced_shape[0];
    rInfo.TotalSize = global_shape[0];

    KRATOS_INFO_IF("HDF5Application", GetEchoLevel() == 2)
        << "Write time \"" << rPath << "\": " << timer.ElapsedSeconds() << std::endl;

    KRATOS_CATCH("Path: \"" + rPath + "\".");
}

template<class TDataType, File::DataTransferMode TDataTransferMode>
void File::ReadDataSetImpl(
    const std::string& rPath,
    TDataType& rData,
    const unsigned StartIndex,
    const unsigned BlockSize)
{
    KRATOS_TRY;

    using type_trait = DataTypeTraits<TDataType>;

    static_assert(type_trait::IsContiguous, "HDF5File can only write contiguous data sets.");

    constexpr auto local_dimension = type_trait::Dimension;

    constexpr auto global_dimension = (local_dimension == 1 ? 1 : 2);

    BuiltinTimer timer;
    // Check that full path exists.
    KRATOS_ERROR_IF_NOT(IsDataSet(rPath))
        << "Path is not a data set: " << rPath << std::endl;

    const auto& file_space_dims = GetDataDimensions(rPath);

    // Check consistency of file's data set dimensions.
    KRATOS_ERROR_IF(file_space_dims.size() != global_dimension)
        << "Invalid data set dimension." << std::endl;
    KRATOS_ERROR_IF(StartIndex + BlockSize > file_space_dims[0])
        << "StartIndex (" << StartIndex << ") + BlockSize (" << BlockSize
        << ") > size of data set (" << file_space_dims[0] << ")." << std::endl;

    std::vector<hsize_t> memory_space_dims(local_dimension);
    type_trait::Shape(rData, memory_space_dims.data(), memory_space_dims.data() + local_dimension);

    if constexpr(local_dimension >= 2) {
        if constexpr(type_trait::template IsDimensionDynamic<1>()) {
            const auto v = std::accumulate(memory_space_dims.begin() + 2, memory_space_dims.end(), hsize_t{1}, std::multiplies<hsize_t>());
            KRATOS_ERROR_IF_NOT(file_space_dims[1] % v == 0) << "Size mismatch with memory space and file space.";
            memory_space_dims[1] = file_space_dims[1] / v;
        }
    }
    memory_space_dims[0] = BlockSize;

    // now reshape the memory space data
    type_trait::Reshape(rData, memory_space_dims.data(), memory_space_dims.data() + local_dimension);

    std::vector<hsize_t> local_reduced_space_dims(file_space_dims.begin(), file_space_dims.end()), local_space_start(global_dimension, 0);
    local_reduced_space_dims[0] = BlockSize;
    local_space_start[0] = StartIndex;

    KRATOS_ERROR_IF_NOT(HasDataType<typename type_trait::PrimitiveType>(rPath))
        << "Data type mismatch at " << rPath << std::endl;

    // Set the data type.
    hid_t dtype_id = Internals::GetPrimitiveH5Type<TDataType>();

    hid_t file_id = GetFileId();

    hid_t dset_id = H5Dopen(file_id, rPath.c_str(), H5P_DEFAULT);
    KRATOS_ERROR_IF(dset_id < 0) << "H5Dopen failed." << std::endl;
    hid_t file_space_id = H5Dget_space(dset_id);
    hid_t mem_space_id = H5Screate_simple(global_dimension, local_reduced_space_dims.data(), nullptr);
    KRATOS_ERROR_IF(H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, local_space_start.data(), nullptr, local_reduced_space_dims.data(), nullptr) < 0)
        << "H5Sselect_hyperslab failed." << std::endl;

    if (!GetDataCommunicator().IsDistributed()) {
        if (type_trait::Size(rData) > 0) {
            KRATOS_ERROR_IF(H5Dread(dset_id, dtype_id, mem_space_id, file_space_id, H5P_DEFAULT, type_trait::GetContiguousData(rData)) < 0)
                << "H5Dread failed." << std::endl;
        } else {
            KRATOS_ERROR_IF(H5Dread(dset_id, dtype_id, mem_space_id, file_space_id, H5P_DEFAULT, nullptr) < 0)
                << "H5Dread failed." << std::endl;
        }
    } else {
        #ifdef KRATOS_USING_MPI
            hid_t dxpl_id = H5Pcreate(H5P_DATASET_XFER);
            if constexpr(TDataTransferMode == DataTransferMode::Collective) {
                H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE);
            } else {
                H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_INDEPENDENT);
            }
            if (type_trait::Size(rData) > 0) {
                KRATOS_ERROR_IF(H5Dread(dset_id, dtype_id, mem_space_id, file_space_id, dxpl_id, type_trait::GetContiguousData(rData)) < 0)
                    << "H5Dread failed." << std::endl;
            } else {
                KRATOS_ERROR_IF(H5Dread(dset_id, dtype_id, mem_space_id, file_space_id, dxpl_id, nullptr) < 0)
                    << "H5Dread failed." << std::endl;
            }
            KRATOS_ERROR_IF(H5Pclose(dxpl_id) < 0) << "H5Pclose failed." << std::endl;
        #else
            KRATOS_ERROR << "HDF5Application is not compiled with MPI.";
        #endif
    }

    KRATOS_ERROR_IF(H5Dclose(dset_id) < 0) << "H5Dclose failed." << std::endl;
    KRATOS_ERROR_IF(H5Sclose(file_space_id) < 0) << "H5Sclose failed." << std::endl;
    KRATOS_ERROR_IF(H5Sclose(mem_space_id) < 0) << "H5Sclose failed." << std::endl;

    KRATOS_INFO_IF("HDF5Application", GetEchoLevel() == 2)
        << "Read time \"" << rPath << "\": " << timer.ElapsedSeconds() << std::endl;

    KRATOS_CATCH("Path: \"" + rPath + "\".");
}

// template instantiations

template void File::GetDataSet<int>(hid_t&, hid_t&, const std::vector<hsize_t>&, const std::string&);
template void File::GetDataSet<double>(hid_t&, hid_t&, const std::vector<hsize_t>&, const std::string&);

#ifndef KRATOS_HDF5_FILE_DATA_SET_METHOD_INSTANTIATION
#define KRATOS_HDF5_FILE_DATA_SET_METHOD_INSTANTIATION(...)                                                                                                                              \
    template KRATOS_API(HDF5_APPLICATION) void File::WriteDataSetImpl<__VA_ARGS__, File::DataTransferMode::Collective>(const std::string&, const __VA_ARGS__&, WriteInfo&);              \
    template KRATOS_API(HDF5_APPLICATION) void File::ReadDataSetImpl<__VA_ARGS__, File::DataTransferMode::Collective>(const std::string&, __VA_ARGS__&, const unsigned, const unsigned); \
    template KRATOS_API(HDF5_APPLICATION) void File::WriteDataSetImpl<__VA_ARGS__, File::DataTransferMode::Independent>(const std::string&, const __VA_ARGS__&, WriteInfo&);             \
    template KRATOS_API(HDF5_APPLICATION) void File::ReadDataSetImpl<__VA_ARGS__, File::DataTransferMode::Independent>(const std::string&, __VA_ARGS__&, const unsigned, const unsigned);\
    template KRATOS_API(HDF5_APPLICATION) void File::WriteDataSet(const std::string&, const __VA_ARGS__&, WriteInfo&);                                                                   \
    template KRATOS_API(HDF5_APPLICATION) void File::WriteDataSetIndependent(const std::string&, const __VA_ARGS__&, WriteInfo&);                                                        \
    template KRATOS_API(HDF5_APPLICATION) void File::ReadDataSet(const std::string&, __VA_ARGS__&, const unsigned, const unsigned);                                                      \
    template KRATOS_API(HDF5_APPLICATION) void File::ReadDataSetIndependent(const std::string&, __VA_ARGS__&, const unsigned, const unsigned);                                           \

#endif

#ifndef KRATOS_HDF5_FILE_ATTRIBUTE_METHOD_INSTANTIATION
#define KRATOS_HDF5_FILE_ATTRIBUTE_METHOD_INSTANTIATION(...)                                                                        \
    template KRATOS_API(HDF5_APPLICATION) void File::WriteAttribute(const std::string&, const std::string&, const __VA_ARGS__&);    \
    template KRATOS_API(HDF5_APPLICATION) void File::ReadAttribute(const std::string&, const std::string&, __VA_ARGS__&);           \

#endif

KRATOS_HDF5_FILE_ATTRIBUTE_METHOD_INSTANTIATION(int);
KRATOS_HDF5_FILE_ATTRIBUTE_METHOD_INSTANTIATION(double);
KRATOS_HDF5_FILE_ATTRIBUTE_METHOD_INSTANTIATION(std::string);
KRATOS_HDF5_FILE_ATTRIBUTE_METHOD_INSTANTIATION(Vector<int>);
KRATOS_HDF5_FILE_ATTRIBUTE_METHOD_INSTANTIATION(Vector<double>);
KRATOS_HDF5_FILE_ATTRIBUTE_METHOD_INSTANTIATION(Matrix<int>);
KRATOS_HDF5_FILE_ATTRIBUTE_METHOD_INSTANTIATION(Matrix<double>);
KRATOS_HDF5_FILE_ATTRIBUTE_METHOD_INSTANTIATION(array_1d<double, 3>);
KRATOS_HDF5_FILE_ATTRIBUTE_METHOD_INSTANTIATION(array_1d<double, 4>);
KRATOS_HDF5_FILE_ATTRIBUTE_METHOD_INSTANTIATION(array_1d<double, 6>);
KRATOS_HDF5_FILE_ATTRIBUTE_METHOD_INSTANTIATION(array_1d<double, 9>);

KRATOS_HDF5_FILE_DATA_SET_METHOD_INSTANTIATION(Vector<int>);
KRATOS_HDF5_FILE_DATA_SET_METHOD_INSTANTIATION(Vector<double>);
KRATOS_HDF5_FILE_DATA_SET_METHOD_INSTANTIATION(Vector<array_1d<double, 3>>);
KRATOS_HDF5_FILE_DATA_SET_METHOD_INSTANTIATION(Vector<array_1d<double, 4>>);
KRATOS_HDF5_FILE_DATA_SET_METHOD_INSTANTIATION(Vector<array_1d<double, 6>>);
KRATOS_HDF5_FILE_DATA_SET_METHOD_INSTANTIATION(Vector<array_1d<double, 9>>);
KRATOS_HDF5_FILE_DATA_SET_METHOD_INSTANTIATION(Matrix<int>);
KRATOS_HDF5_FILE_DATA_SET_METHOD_INSTANTIATION(Matrix<double>);

#undef KRATOS_HDF5_FILE_DATA_SET_METHOD_INSTANTIATION
#undef KRATOS_HDF5_FILE_ATTRIBUTE_METHOD_INSTANTIATION

} // namespace HDF5.
} // namespace Kratos.
