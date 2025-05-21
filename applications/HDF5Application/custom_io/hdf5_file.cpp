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
#ifdef KRATOS_USING_MPI
#include "mpi/includes/mpi_data_communicator.h"
#endif

// Application includes
#include "custom_utilities/data_type_utilities.h"

// Include base h
#include "hdf5_file.h"

#ifndef KRATOS_CHECK_FOR_HDF5_ERROR
#define KRATOS_CHECK_FOR_HDF5_ERROR(METHOD, /*int*/ RETURN_VALUE)                                           \
    if (RETURN_VALUE < 0) {                                                                                 \
        std::stringstream message;                                                                          \
        message << "call to " << #METHOD << " returned error code " << RETURN_VALUE << " with message:\n";  \
        auto callback = [](unsigned depth,                                                                  \
                        const H5E_error_t* p_error_descriptor,                                              \
                        void* p_stream) -> herr_t {                                                         \
            std::stringstream& r_stream = *static_cast<std::stringstream*>(p_stream);                       \
            if (not depth) r_stream << p_error_descriptor->desc << "\n";                                    \
            r_stream << "in file " << p_error_descriptor->file_name                                         \
                    << ":" << p_error_descriptor->line                                                      \
                    << " in function " << p_error_descriptor->func_name << "\n";                            \
            return 0;                                                                                       \
        };                                                                                                  \
        H5Ewalk(H5E_DEFAULT, H5E_WALK_UPWARD, callback, static_cast<void*>(&message));                      \
        KRATOS_ERROR << message.str();                                                                      \
    }
#endif

#ifndef KRATOS_HDF5_CALL_WITH_RETURN
#define KRATOS_HDF5_CALL_WITH_RETURN(RETURN_VALUE, METHOD, ...) \
    RETURN_VALUE = METHOD(__VA_ARGS__);                         \
    KRATOS_CHECK_FOR_HDF5_ERROR(METHOD, RETURN_VALUE)           \

#endif

#ifndef KRATOS_HDF5_CALL
#define KRATOS_HDF5_CALL(METHOD, ...)              \
    {                                              \
        const auto error = METHOD(__VA_ARGS__);    \
        KRATOS_CHECK_FOR_HDF5_ERROR(METHOD, error) \
    }

#endif

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

    // Unset error handling callbacks set by dependencies.
    // Some dependencies (such as h5py) may set callbacks via the legacy HDF5 API,
    // which breaks Kratos HDF5 calls. We don't do any kind of error handling
    // apart from reporting, so there's no behavior change on our end, but if
    // dependencies issue HDF5 calls after Kratos, they might run into trouble.
    KRATOS_HDF5_CALL(H5Eset_auto, H5E_DEFAULT, NULL, NULL)

    Parameters default_params(R"(
            {
                "file_name" : "PLEASE_SPECIFY_HDF5_FILENAME",
                "file_access_mode": "exclusive",
                "file_driver": "sec2",
                "echo_level" : 0
            })");

    Settings.RecursivelyAddMissingParameters(default_params);

    mFileName = Settings["file_name"].GetString();
    KRATOS_ERROR_IF(mFileName == "PLEASE_SPECIFY_HDF5_FILENAME") << "Invalid file name: " << mFileName << std::endl;

    hid_t fapl_id;
    KRATOS_HDF5_CALL_WITH_RETURN(fapl_id, H5Pcreate, H5P_FILE_ACCESS)

    std::string file_driver = Settings["file_driver"].GetString();
    SetFileDriver(file_driver, fapl_id);

    std::string file_access_mode = Settings["file_access_mode"].GetString();
    if (file_access_mode == "exclusive") {
        mFileId = H5Fcreate(mFileName.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, fapl_id);
    } else if (file_access_mode == "truncate") {
        mFileId = H5Fcreate(mFileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
    } else {
        // Save old error handler
        herr_t (*old_func)(hid_t, void*);
        void *old_client_data;

        KRATOS_HDF5_CALL(H5Eget_auto, H5E_DEFAULT, &old_func, &old_client_data)

        // Turn off error handling
        KRATOS_HDF5_CALL(H5Eset_auto, H5E_DEFAULT, NULL, NULL)

        // follwoing will fail if the hdf5 file is not found, and will print a failure msg to the output.
        // if "file_access_mode" is "read_write", then it is ok to fail the following call because,
        // if it is failed then the file will be created in the subsequent section.
        htri_t is_hdf5 = H5Fis_hdf5(mFileName.c_str());

        // Restore previous error handler
        KRATOS_HDF5_CALL(H5Eset_auto, H5E_DEFAULT, old_func, old_client_data)

        if (file_access_mode == "read_only") {
            KRATOS_ERROR_IF(is_hdf5 <= 0 && file_driver != "core") << "Invalid HDF5 file: " << mFileName << std::endl;
            mFileId = H5Fopen(mFileName.c_str(), H5F_ACC_RDONLY, fapl_id);
        } else if (file_access_mode == "read_write") {
            if (is_hdf5 <= 0) {
                // creates the hdf5 file if the file is not found
                mFileId = H5Fcreate(mFileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
            } else {
                // open the existing hdf5 file if file is found
                mFileId = H5Fopen(mFileName.c_str(), H5F_ACC_RDWR, fapl_id);
            }
        } else {
            KRATOS_ERROR << "Invalid \"file_access_mode\": " << file_access_mode
                         << std::endl;
        }
    }

    KRATOS_ERROR_IF(mFileId < 0) << "Failed to open file: " << mFileName << " [ hdf5 error code = " << mFileId << " ].\n";

    KRATOS_HDF5_CALL(H5Pclose, fapl_id)

    mEchoLevel = Settings["echo_level"].GetInt();

    KRATOS_CATCH("");
}

File::File(File&& rOther)
{
    mpDataCommunicator = rOther.mpDataCommunicator;
    mFileName = std::move(rOther.mFileName);
    mFileId = rOther.mFileId;
    rOther.mFileId = -1;
    mEchoLevel = rOther.mEchoLevel;
    rOther.mEchoLevel = 0;
}

File& File::operator=(File&& rOther)
{
    mpDataCommunicator = rOther.mpDataCommunicator;
    mFileName = std::move(rOther.mFileName);
    mFileId = rOther.mFileId;
    rOther.mFileId = -1;
    mEchoLevel = rOther.mEchoLevel;
    rOther.mEchoLevel = 0;
    return *this;
}

File::~File()
{
    if (0 <= mFileId) {
        H5Fclose(mFileId);
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

    std::vector<std::string> split_path = StringUtilities::SplitStringByDelimiter(rPath, '/');
    split_path.erase(std::remove_if(split_path.begin(), split_path.end(), [](const std::string& s) {return (s.size() == 0);}));
    std::string sub_path;

    for (const auto& r_link: split_path) {
        sub_path += '/' + r_link;

        htri_t link_found;
        KRATOS_HDF5_CALL_WITH_RETURN(link_found, H5Lexists, mFileId, sub_path.c_str(), H5P_DEFAULT)
        if (!link_found) {
            return false;
        }

        htri_t object_found;
        KRATOS_HDF5_CALL_WITH_RETURN(object_found, H5Oexists_by_name, mFileId, sub_path.c_str(), H5P_DEFAULT)
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

    if (!HasPath(rPath)) {// Expects a valid path.
        return false;
    }

    H5O_info_t object_info;
    KRATOS_HDF5_CALL(H5Oget_info_by_name, mFileId, rPath.c_str(), &object_info, H5P_DEFAULT)
    return object_info.type == H5O_TYPE_GROUP;

    KRATOS_CATCH("");
}

bool File::IsDataSet(const std::string& rPath) const
{
    KRATOS_TRY;

    if (!HasPath(rPath)) {// Expects a valid path.
        return false;
    }

    H5O_info_t object_info;
    KRATOS_HDF5_CALL(H5Oget_info_by_name, mFileId, rPath.c_str(), &object_info, H5P_DEFAULT)
    return object_info.type == H5O_TYPE_DATASET;

    KRATOS_CATCH("");
}

bool File::HasAttribute(
    const std::string& rObjectPath,
    const std::string& rName) const
{
    KRATOS_TRY;

    htri_t status;
    KRATOS_HDF5_CALL_WITH_RETURN(status, H5Aexists_by_name, mFileId, rObjectPath.c_str(), rName.c_str(), H5P_DEFAULT)
    return status > 0;

    KRATOS_CATCH("");
}

template<class TDataType>
bool File::HasAttributeType(
    const std::string& rObjectPath,
    const std::string& rName) const
{
    KRATOS_TRY

    auto mem_type_id = Internals::GetPrimitiveH5Type<TDataType>();

    hid_t attr_id;
    KRATOS_HDF5_CALL_WITH_RETURN(attr_id, H5Aopen_by_name, mFileId, rObjectPath.c_str(), rName.c_str(), H5P_DEFAULT, H5P_DEFAULT)

    hid_t attr_type_id;
    KRATOS_HDF5_CALL_WITH_RETURN(attr_type_id, H5Aget_type, attr_id)

    htri_t is_valid_type;
    KRATOS_HDF5_CALL_WITH_RETURN(is_valid_type, H5Tequal, mem_type_id, attr_type_id)

    KRATOS_HDF5_CALL(H5Tclose, attr_type_id)
    KRATOS_HDF5_CALL(H5Aclose, attr_id)

    return is_valid_type != 0;

    KRATOS_CATCH("");
}

std::vector<hsize_t> File::GetAttributeDimensions(
    const std::string& rObjectPath,
    const std::string& rName) const
{
    KRATOS_TRY

    hid_t attr_id;
    KRATOS_HDF5_CALL_WITH_RETURN(attr_id, H5Aopen_by_name, mFileId, rObjectPath.c_str(), rName.c_str(), H5P_DEFAULT, H5P_DEFAULT)

    hid_t space_id;
    KRATOS_HDF5_CALL_WITH_RETURN(space_id, H5Aget_space, attr_id)
    int ndims;
    KRATOS_HDF5_CALL_WITH_RETURN(ndims, H5Sget_simple_extent_ndims, space_id)
    std::vector<hsize_t> shape(ndims);
    KRATOS_HDF5_CALL(H5Sget_simple_extent_dims, space_id, shape.data(), nullptr)
    KRATOS_HDF5_CALL(H5Sclose, space_id)
    KRATOS_HDF5_CALL(H5Aclose, attr_id)

    return shape;

    KRATOS_CATCH("");
}

void File::DeleteAttribute(
    const std::string& rObjectPath,
    const std::string& rName)
{
    KRATOS_TRY;

    KRATOS_HDF5_CALL(H5Adelete_by_name, mFileId, rObjectPath.c_str(), rName.c_str(), H5P_DEFAULT)

    KRATOS_CATCH("");
}

std::vector<std::string> File::GetAttributeNames(const std::string& rObjectPath) const
{
    KRATOS_TRY;

    constexpr unsigned max_ssize = 100;
    std::array<char, max_ssize> buffer;

    // Get number of attributes.
    hid_t object_id;
    KRATOS_HDF5_CALL_WITH_RETURN(object_id, H5Oopen, mFileId, rObjectPath.c_str(), H5P_DEFAULT)

    H5O_info_t object_info;
    KRATOS_HDF5_CALL(H5Oget_info, object_id, &object_info)

    hsize_t num_attrs = object_info.num_attrs;
    std::vector<std::string> names(num_attrs);

    for (hsize_t i = 0; i < num_attrs; ++i) {
        // Get size of name.
        ssize_t ssize;
        KRATOS_HDF5_CALL_WITH_RETURN(ssize, H5Aget_name_by_idx, mFileId, rObjectPath.c_str(), H5_INDEX_CRT_ORDER, H5_ITER_INC, i, buffer.data(), max_ssize, H5P_DEFAULT)
        KRATOS_ERROR_IF(ssize > max_ssize)
            << "Attribute name size exceeds "
            << max_ssize << std::endl;

        names[i].resize(ssize);
        std::copy_n(buffer.begin(), ssize, names[i].begin());
    }

    KRATOS_HDF5_CALL(H5Oclose, object_id)

    return names;

    KRATOS_CATCH("");
}

void File::CreateGroup(const std::string& rPath)
{
    KRATOS_TRY;

    hid_t group_id;
    KRATOS_HDF5_CALL_WITH_RETURN(group_id, H5Gcreate, mFileId,rPath.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)
    KRATOS_HDF5_CALL(H5Gclose, group_id)

    KRATOS_CATCH("");
}

std::vector<std::string> File::GetLinkNames(const std::string& rGroupPath) const
{
    KRATOS_TRY;

    constexpr unsigned max_ssize = 100;
    char buffer[max_ssize];
    // Get number of links.
    hid_t group_id;
    KRATOS_HDF5_CALL_WITH_RETURN(group_id, H5Gopen, mFileId, rGroupPath.c_str(), H5P_DEFAULT)

    H5G_info_t group_info;
    KRATOS_HDF5_CALL(H5Gget_info, group_id, &group_info)
    hsize_t num_links = group_info.nlinks;
    std::vector<std::string> names(num_links);

    for (hsize_t i=0; i < num_links; ++i) {
        // Get size of name.
        ssize_t ssize;
        KRATOS_HDF5_CALL_WITH_RETURN(ssize, H5Lget_name_by_idx, mFileId, rGroupPath.c_str(), H5_INDEX_NAME,
                                    H5_ITER_INC, i, buffer, max_ssize, H5P_DEFAULT)
        KRATOS_ERROR_IF(ssize > max_ssize) << "Link name size exceeds "
                                           << max_ssize << std::endl;
        names[i].resize(ssize);
        std::copy_n(buffer, ssize, names[i].begin());
    }
    KRATOS_HDF5_CALL(H5Gclose, group_id)
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
            names.emplace_back(std::move(r_name));
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

    auto split_path = StringUtilities::SplitStringByDelimiter(rPath, '/');
    split_path.erase(std::remove_if(split_path.begin(), split_path.end(), [](const std::string& s) {return (s.size() == 0);}));
    std::string sub_path;
    for (const auto& r_link: split_path) {
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

    using TypeTraits = DataTypeTraits<TDataType>;

    static_assert(TypeTraits::IsContiguous, "Attributes needs to have data contiguous in memory.");

    constexpr auto local_dimension = TypeTraits::Dimension;

    BuiltinTimer timer;
    hid_t type_id, space_id, attr_id;

    if (HasAttribute(rObjectPath, rName)) {
        DeleteAttribute(rObjectPath, rName);
    }

    if constexpr(local_dimension == 0) {
        KRATOS_HDF5_CALL_WITH_RETURN(space_id, H5Screate, H5S_SCALAR)
    } else {
        std::vector<hsize_t> shape(local_dimension);
        TypeTraits::Shape(rValue, shape.data(), shape.data() + local_dimension);
        KRATOS_HDF5_CALL_WITH_RETURN(space_id, H5Screate_simple, local_dimension, shape.data(), nullptr)
    }

    type_id = Internals::GetPrimitiveH5Type<TDataType>();
    KRATOS_HDF5_CALL_WITH_RETURN(attr_id, H5Acreate_by_name, mFileId, rObjectPath.c_str(), rName.c_str(), type_id, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)
    KRATOS_HDF5_CALL(H5Awrite, attr_id, type_id, TypeTraits::GetContiguousData(rValue))
    KRATOS_HDF5_CALL(H5Sclose, space_id)
    KRATOS_HDF5_CALL(H5Aclose, attr_id)

    KRATOS_INFO_IF("HDF5Application", GetEchoLevel() == 2)
        << "Write time \"" << rObjectPath << '/' << rName
        << "\": " << timer.ElapsedSeconds() << std::endl;

    KRATOS_CATCH("Path: \"" + rObjectPath + '/' + rName + "\".");
}

void File::WriteAttribute(
    const std::string& rObjectPath,
    const Parameters Attributes)
{
    KRATOS_TRY

    for (auto itr = Attributes.begin(); itr != Attributes.end(); ++itr) {
        const auto attrib_name = itr.name();

        if (itr->IsInt()) {
            WriteAttribute(rObjectPath, attrib_name, itr->GetInt());
        } else if (itr->IsDouble()) {
            WriteAttribute(rObjectPath, attrib_name, itr->GetDouble());
        } else if (itr->IsString()) {
            WriteAttribute(rObjectPath, attrib_name, itr->GetString());
        } else if (itr->IsMatrix()) {
            WriteAttribute(rObjectPath, attrib_name, itr->GetMatrix());
        } else if (itr->IsVector()) {
            WriteAttribute(rObjectPath, attrib_name, itr->GetVector());
            bool is_int_array = true;
            bool is_double_array = true;

            for (const auto& v : *itr) {
                is_int_array = is_int_array && v.IsInt();
                is_double_array = is_double_array && v.IsDouble();
            }

            if (is_int_array) {
                Vector<int> values(itr->size());
                for (unsigned int i = 0; i < values.size(); ++i) {
                    values[i] = itr->GetArrayItem(i).GetInt();
                }
                WriteAttribute(rObjectPath, attrib_name, values);
            } else if (is_double_array) {
                WriteAttribute(rObjectPath, attrib_name, itr->GetVector());
            } else {
                KRATOS_ERROR << "Only supports int or double arrays.";
            }
        } else {
            KRATOS_ERROR << "Only supports following data types:\n\tint\n\tdouble\n\tstring\n\tVector\n\tMatrix\n\tIntArray";
        }
    }

    KRATOS_CATCH("");
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
    KRATOS_HDF5_CALL_WITH_RETURN(dset_id, H5Dopen, mFileId, rPath.c_str(), H5P_DEFAULT)
    KRATOS_HDF5_CALL_WITH_RETURN(dspace_id, H5Dget_space, dset_id)
    KRATOS_HDF5_CALL_WITH_RETURN(ndims, H5Sget_simple_extent_ndims, dspace_id)

    std::vector<hsize_t> dims(ndims);
    KRATOS_HDF5_CALL(H5Sget_simple_extent_dims, dspace_id, dims.data(), nullptr)
    KRATOS_HDF5_CALL(H5Sclose, dspace_id)
    KRATOS_HDF5_CALL(H5Dclose, dset_id)

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
    KRATOS_HDF5_CALL(H5Fflush, mFileId, H5F_SCOPE_GLOBAL)
}

void File::Close()
{
    if (0 <= mFileId) {
        KRATOS_HDF5_CALL(H5Fclose, mFileId)
        mFileId = -1;
    } else {
        KRATOS_WARNING("Invalid file handle") << "Attempt to close an invalid file" << std::endl;
    }
}

unsigned File::GetFileSize() const
{
    hsize_t size;
    KRATOS_HDF5_CALL(H5Fget_filesize, mFileId, &size)
    return size;
}

std::string File::GetFileName() const
{
    return mFileName;
}

int File::GetEchoLevel() const
{
    return mEchoLevel;
}

void File::SetEchoLevel(int Level)
{
    mEchoLevel = Level;
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
    TDataType& rValue) const
{
    KRATOS_TRY;

    using TypeTraits = DataTypeTraits<TDataType>;

    static_assert(TypeTraits::IsContiguous, "Attribute data should be contiguous in memory.");

    constexpr auto local_dimension = TypeTraits::Dimension;

    BuiltinTimer timer;

    KRATOS_ERROR_IF_NOT(HasAttributeType<TDataType>(rObjectPath, rName))
        << "Memory and file data types are different." << std::endl;

    const auto& shape = GetAttributeDimensions(rObjectPath, rName);
    KRATOS_ERROR_IF(shape.size() != local_dimension)
        << "Attribute \"" << rName << "\" has dimension mismatch [ memory dimension = "
        << local_dimension << ", file dimension = " << shape.size()  << " ].\n";

    auto mem_type_id = Internals::GetPrimitiveH5Type<TDataType>();

    hid_t attr_id;
    KRATOS_HDF5_CALL_WITH_RETURN(attr_id, H5Aopen_by_name, mFileId, rObjectPath.c_str(), rName.c_str(), H5P_DEFAULT, H5P_DEFAULT)

    TypeTraits::Reshape(rValue, shape);

    // Read attribute.
    KRATOS_HDF5_CALL(H5Aread, attr_id, mem_type_id, TypeTraits::GetContiguousData(rValue))
    KRATOS_HDF5_CALL(H5Aclose, attr_id)

    KRATOS_INFO_IF("HDF5Application", GetEchoLevel() == 2)
        << "Read time \"" << rObjectPath << '/' << rName
        << "\": " << timer.ElapsedSeconds() << std::endl;

    KRATOS_CATCH("Path: \"" + rObjectPath + '/' + rName + "\".");
}

Parameters File::ReadAttribute(const std::string& rObjectPath) const
{
    KRATOS_TRY

    Parameters result;

    const auto attribute_names = GetAttributeNames(rObjectPath);

    for (const auto& attribute_name : attribute_names) {
        const auto shape = GetAttributeDimensions(rObjectPath, attribute_name);

        if (HasAttributeType<int>(rObjectPath, attribute_name)) {
            if (shape.empty()) {
                int value;
                ReadAttribute(rObjectPath, attribute_name, value);
                result.AddInt(attribute_name, value);
            } else if (shape.size() == 1) {
                Vector<int> values;
                ReadAttribute(rObjectPath, attribute_name, values);
                result.AddEmptyArray(attribute_name);
                for (const auto v : values) {
                    result[attribute_name].Append(v);
                }
            } else {
                KRATOS_ERROR << "Unsupported dimension for int data type.";
            }
        } else if (HasAttributeType<double>(rObjectPath, attribute_name)) {
            if (shape.empty()) {
                double value;
                ReadAttribute(rObjectPath, attribute_name, value);
                result.AddDouble(attribute_name, value);
            } else if (shape.size() == 1) {
                Vector<double> values;
                ReadAttribute(rObjectPath, attribute_name, values);
                result.AddVector(attribute_name, values);
            } else if (shape.size() == 2) {
                Matrix<double> values;
                ReadAttribute(rObjectPath, attribute_name, values);
                result.AddMatrix(attribute_name, values);
            } else {
                KRATOS_ERROR << "Unsupported dimension for double data type.";
            }
        } else if (HasAttributeType<char>(rObjectPath, attribute_name)) {
            if (shape.size() != 0) {
                std::string values;
                ReadAttribute(rObjectPath, attribute_name, values);
                result.AddString(attribute_name, values);
            } else {
                KRATOS_ERROR << "Unsupported std::string type.";
            }
        } else {
            KRATOS_ERROR << "Unsupported data type.";
        }
    }

    return result;

    KRATOS_CATCH("");
}

template<class TDataType>
void File::ReadDataSet(
    const std::string& rPath,
    TDataType& rData,
    const unsigned StartIndex,
    const unsigned BlockSize) const
{
    ReadDataSetImpl<TDataType, DataTransferMode::Collective>(rPath, rData, StartIndex, BlockSize);
}

template<class TDataType>
void File::ReadDataSetIndependent(
    const std::string& rPath,
    TDataType& rData,
    const unsigned StartIndex,
    const unsigned BlockSize) const
{
    ReadDataSetImpl<TDataType, DataTransferMode::Independent>(rPath, rData, StartIndex, BlockSize);
}

unsigned File::GetOpenObjectsCount() const
{
    KRATOS_TRY;
    ssize_t num_open_objects;
    KRATOS_HDF5_CALL_WITH_RETURN(num_open_objects, H5Fget_obj_count, mFileId, H5F_OBJ_ALL)
    return num_open_objects;
    KRATOS_CATCH("");
}

hid_t File::GetFileId() const
{
    return mFileId;
}

template<class TDataType>
bool File::HasDataType(const std::string& rPath) const
{
    KRATOS_TRY;

    hid_t dset_id, dtype_id;

    KRATOS_HDF5_CALL_WITH_RETURN(dset_id, H5Dopen, GetFileId(), rPath.c_str(), H5P_DEFAULT)
    KRATOS_HDF5_CALL_WITH_RETURN(dtype_id, H5Dget_type, dset_id);
    H5T_class_t type = H5Tget_class(dtype_id);

    KRATOS_ERROR_IF(type == H5T_NO_CLASS) << "Invalid data type." << std::endl;
    KRATOS_HDF5_CALL(H5Tclose, dtype_id)
    KRATOS_HDF5_CALL(H5Dclose, dset_id)

    if constexpr(std::is_same_v<TDataType, int>) {
        return type == H5T_INTEGER;
    } else if constexpr(std::is_same_v<TDataType, double>) {
        return type == H5T_FLOAT;
    } else if constexpr(std::is_same_v<TDataType, char>) {
        return type == H5T_INTEGER;
    } else {
        static_assert(!std::is_same_v<TDataType, TDataType>, "Unsupported data type.");
    }

    KRATOS_CATCH("");
}

hid_t File::OpenExistingDataSet(const std::string& rPath)
{
    hid_t dset_id;
    KRATOS_HDF5_CALL_WITH_RETURN(dset_id, H5Dopen, GetFileId(), rPath.c_str(), H5P_DEFAULT)
    return dset_id;
}

void File::SetFileDriver(const std::string& rDriver, hid_t FaplId) const
{
    KRATOS_TRY;
#if (defined(_WIN32) || defined(_WIN64))
    KRATOS_ERROR_IF(rDriver != "windows")
        << "Unsupported (Windows) \"file_driver\": " << rDriver << std::endl;
    KRATOS_HDF5_CALL(H5Pset_fapl_windows, FaplId)
#else
    if (rDriver == "sec2") {
        KRATOS_HDF5_CALL(H5Pset_fapl_sec2, FaplId)
    } else if (rDriver == "stdio") {
        KRATOS_HDF5_CALL(H5Pset_fapl_stdio, FaplId)
    } else if (rDriver == "core") {
        KRATOS_HDF5_CALL(H5Pset_fapl_core, FaplId, 1000000, 0)
    } else if (rDriver == "mpio") {
#if defined(KRATOS_USING_MPI)
        KRATOS_HDF5_CALL(H5Pset_fapl_mpio, FaplId, MPIDataCommunicator::GetMPICommunicator(GetDataCommunicator()), MPI_INFO_NULL)
#else
        KRATOS_ERROR
            << "Kratos must be built with MPI for \"file_driver\"=\"mpio\"."
            << std::endl;
#endif
    } else {
        KRATOS_ERROR << "Unsupported \"file_driver\": " << rDriver << std::endl;
    }
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

    using TypeTraits = DataTypeTraits<TDataType>;

    static_assert(TypeTraits::IsContiguous, "HDF5File can only write contiguous data sets.");

    constexpr auto local_dimension = TypeTraits::Dimension;

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
    TypeTraits::Shape(rData, local_shape.data(), local_shape.data() + local_dimension);

    const hsize_t number_of_local_primitive_data_values = TypeTraits::Size(rData);

    const auto& r_data_communicator = GetDataCommunicator();
    // get the maximized dimensions of the underlying data. Max is taken because,
    // there can be empty ranks which will give wrong sizes in the case of dynamic
    // data types.
    if constexpr(local_dimension >= 2) {
        if constexpr(TypeTraits::template IsDimensionDynamic<1>()) {
            // this is the matrix version. Hence it is better to get the max size
            // from all ranks.
            const auto max_size = r_data_communicator.MaxAll(static_cast<unsigned int>(local_shape[1]));

            // now check every non-empty ranks have the same sizes since this dimension
            // is a dynamic dimension.
            KRATOS_ERROR_IF(number_of_local_primitive_data_values > 0 && max_size != local_shape[1])
                << "Mismatching shapes found in different ranks. All ranks should have the same shapes in data sets.";

            local_shape[1] = max_size;
        }
    }

    // local_reduced_shape holds the max 2d flattened shape if the local_shape dimensions
    // are higher than 2.
    std::array<hsize_t, global_dimension> local_reduced_shape{}, local_shape_start{}, global_shape{};

    // get total number of items to be written in the data set to the first dimension.
    global_shape[0] = r_data_communicator.SumAll(static_cast<unsigned int>(local_shape[0]));
    local_reduced_shape[0] = local_shape[0];

    if constexpr(global_dimension > 1) {
        // flattens higher dimensions into one since we write matrices which is the highest dimension
        // supported by paraview for visualization
        global_shape[1] = std::accumulate(local_shape.begin() + 1, local_shape.end(), hsize_t{1}, std::multiplies<hsize_t>());
        local_reduced_shape[1] = global_shape[1];
    }

    // Set the data type.
    hid_t dtype_id = Internals::GetPrimitiveH5Type<TDataType>();

    hid_t dset_id{}, fspace_id{};
    if (!HasPath(rPath)) {
        // Create and write the data set.
        KRATOS_HDF5_CALL_WITH_RETURN(fspace_id, H5Screate_simple, global_shape.size(), global_shape.data(), nullptr)
        KRATOS_HDF5_CALL_WITH_RETURN(dset_id, H5Dcreate, GetFileId(), rPath.c_str(), dtype_id, fspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)
    } else {
        // if it is existing, get it.
        KRATOS_ERROR_IF_NOT(HasDataType<typename TypeTraits::PrimitiveType>(rPath))
            << "Wrong scalar data type: " << rPath << std::endl;
        const auto& data_dimensions = Internals::GetDataDimensions(*this, rPath);
        KRATOS_ERROR_IF_NOT(data_dimensions.size() == global_dimension)
            << "Wrong number of dimensions [ file number of data dimensions = "
            << data_dimensions.size() << ", memory data number of dimensions = "
            << global_dimension << " ].\n";
        for (IndexType i = 0; i < global_dimension; ++i) {
            KRATOS_ERROR_IF_NOT(data_dimensions[i] == global_shape[i])
                << "Wrong dimensional value at index = " << i
                << " [ data dimension in file = " << data_dimensions[i]
                << ", memory dimension = " << global_shape[i] << " ].\n";
        }

        KRATOS_HDF5_CALL_WITH_RETURN(dset_id, OpenExistingDataSet, rPath)
        KRATOS_HDF5_CALL_WITH_RETURN(fspace_id, H5Dget_space, dset_id)
    }

    // here onwards the procedure differs shared memory and distributed memeory runs.
    // The same steps need to be applied if the HDF5Application is compiled with serial hdf5lib
    // and if the HDF5Application is compiled with mpi hdf5lib but running in shared memory parallelization
    // only.
    // Different steps has to be taken if it is run in a distributed memory environment.

    if (!r_data_communicator.IsDistributed()) {
        if (number_of_local_primitive_data_values > 0) {
            KRATOS_HDF5_CALL(H5Dwrite, dset_id, dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, TypeTraits::GetContiguousData(rData))
        } else {
            KRATOS_HDF5_CALL(H5Dwrite, dset_id, dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, nullptr)
        }
    } else {
        #ifdef KRATOS_USING_MPI
            if constexpr(TDataTransferMode == DataTransferMode::Collective) {
                local_shape_start[0] = r_data_communicator.ScanSum(static_cast<unsigned int>(local_reduced_shape[0])) - local_reduced_shape[0];
            }

            bool write_data = TDataTransferMode == DataTransferMode::Collective || number_of_local_primitive_data_values > 0;
            #if H5_VERS_MAJOR < 2 && ((H5_VERS_MINOR == 14 && H5_VERS_RELEASE < 2) ||  H5_VERS_MINOR < 14)
                /**
                 *  Until hdf5 1.14.2, if someone tries to read/write empty containers from every rank, it throws an error. This is
                 *  fixed in the later versions.
                 *
                 * TODO: Remove this pre-compiler directive once we move to compatible versions.
                 */
                write_data &= r_data_communicator.SumAll(static_cast<std::size_t>(number_of_local_primitive_data_values)) > 0;
            #endif

            if (write_data) {
                hid_t dxpl_id, mspace_id;
                KRATOS_HDF5_CALL_WITH_RETURN(dxpl_id, H5Pcreate, H5P_DATASET_XFER)
                if constexpr(TDataTransferMode == DataTransferMode::Collective) {
                    KRATOS_HDF5_CALL(H5Pset_dxpl_mpio, dxpl_id, H5FD_MPIO_COLLECTIVE)
                } else {
                    KRATOS_HDF5_CALL(H5Pset_dxpl_mpio, dxpl_id, H5FD_MPIO_INDEPENDENT)
                }

                // select the local hyperslab
                KRATOS_HDF5_CALL(H5Sselect_hyperslab, fspace_id, H5S_SELECT_SET, local_shape_start.data(), nullptr, local_reduced_shape.data(), nullptr)
                KRATOS_HDF5_CALL_WITH_RETURN(mspace_id, H5Screate_simple, global_dimension, local_reduced_shape.data(), nullptr)
                if (number_of_local_primitive_data_values > 0) {
                    KRATOS_HDF5_CALL(H5Dwrite, dset_id, dtype_id, mspace_id, fspace_id, dxpl_id, TypeTraits::GetContiguousData(rData))
                } else {
                    KRATOS_HDF5_CALL(H5Dwrite, dset_id, dtype_id, mspace_id, fspace_id, dxpl_id, nullptr)
                }

                KRATOS_HDF5_CALL(H5Pclose, dxpl_id)
                KRATOS_HDF5_CALL(H5Sclose, mspace_id)
            }
        #else
            KRATOS_ERROR << "HDFApplication is not compiled with MPI enabled";
        #endif
    }

    KRATOS_HDF5_CALL(H5Sclose, fspace_id)
    KRATOS_HDF5_CALL(H5Dclose, dset_id)

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
    const unsigned BlockSize) const
{
    KRATOS_TRY;

    using TypeTraits = DataTypeTraits<TDataType>;

    static_assert(TypeTraits::IsContiguous, "HDF5File can only write contiguous data sets.");

    constexpr auto local_dimension = TypeTraits::Dimension;

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
    TypeTraits::Shape(rData, memory_space_dims.data(), memory_space_dims.data() + local_dimension);

    if constexpr(local_dimension >= 2) {
        if constexpr(TypeTraits::template IsDimensionDynamic<1>()) {
            const auto v = std::accumulate(memory_space_dims.begin() + 2, memory_space_dims.end(), hsize_t{1}, std::multiplies<hsize_t>());
            KRATOS_ERROR_IF_NOT(file_space_dims[1] % v == 0) << "Size mismatch with memory space and file space.";
            memory_space_dims[1] = file_space_dims[1] / v;
        }
    }
    memory_space_dims[0] = BlockSize;

    // now reshape the memory space data
    TypeTraits::Reshape(rData, memory_space_dims.data(), memory_space_dims.data() + local_dimension);

    std::vector<hsize_t> local_reduced_space_dims(file_space_dims.begin(), file_space_dims.end()), local_space_start(global_dimension, 0);
    local_reduced_space_dims[0] = BlockSize;
    local_space_start[0] = StartIndex;

    KRATOS_ERROR_IF_NOT(HasDataType<typename TypeTraits::PrimitiveType>(rPath))
        << "Data type mismatch at " << rPath << std::endl;

    // Set the data type.
    hid_t dtype_id = Internals::GetPrimitiveH5Type<TDataType>();

    hid_t file_id = GetFileId();

    hid_t dset_id, file_space_id, mem_space_id;
    KRATOS_HDF5_CALL_WITH_RETURN(dset_id, H5Dopen, file_id, rPath.c_str(), H5P_DEFAULT)
    KRATOS_HDF5_CALL_WITH_RETURN(file_space_id, H5Dget_space, dset_id)
    KRATOS_HDF5_CALL_WITH_RETURN(mem_space_id, H5Screate_simple, global_dimension, local_reduced_space_dims.data(), nullptr)
    KRATOS_HDF5_CALL(H5Sselect_hyperslab, file_space_id, H5S_SELECT_SET, local_space_start.data(), nullptr, local_reduced_space_dims.data(), nullptr)

    if (!GetDataCommunicator().IsDistributed()) {
        if (TypeTraits::Size(rData) > 0) {
            KRATOS_HDF5_CALL(H5Dread, dset_id, dtype_id, mem_space_id, file_space_id, H5P_DEFAULT, TypeTraits::GetContiguousData(rData))
        } else {
            KRATOS_HDF5_CALL(H5Dread, dset_id, dtype_id, mem_space_id, file_space_id, H5P_DEFAULT, nullptr)
        }
    } else {
        #ifdef KRATOS_USING_MPI
            bool read_data  = true;
            #if H5_VERS_MAJOR < 2 && ((H5_VERS_MINOR == 14 && H5_VERS_RELEASE < 2) ||  H5_VERS_MINOR < 14)
                /**
                 *  Until hdf5 1.14.2, if someone tries to read/write empty containers from every rank, it throws an error. This is
                 *  fixed in the later versions.
                 *
                 * TODO: Remove this pre-compiler directive once we move to compatible versions.
                 */
                read_data &= mpDataCommunicator->SumAll(TypeTraits::Size(rData)) > 0;
            #endif

            if (read_data) {
                hid_t dxpl_id;
                KRATOS_HDF5_CALL_WITH_RETURN(dxpl_id, H5Pcreate, H5P_DATASET_XFER)
                if constexpr(TDataTransferMode == DataTransferMode::Collective) {
                    KRATOS_HDF5_CALL(H5Pset_dxpl_mpio, dxpl_id, H5FD_MPIO_COLLECTIVE)
                } else {
                    KRATOS_HDF5_CALL(H5Pset_dxpl_mpio, dxpl_id, H5FD_MPIO_INDEPENDENT)
                }
                if (TypeTraits::Size(rData) > 0) {
                    KRATOS_HDF5_CALL(H5Dread, dset_id, dtype_id, mem_space_id, file_space_id, dxpl_id, TypeTraits::GetContiguousData(rData))
                } else {
                    KRATOS_HDF5_CALL(H5Dread, dset_id, dtype_id, mem_space_id, file_space_id, dxpl_id, nullptr)
                }
                KRATOS_HDF5_CALL(H5Pclose, dxpl_id)
            }
        #else
            KRATOS_ERROR << "HDF5Application is not compiled with MPI.";
        #endif
    }

    KRATOS_HDF5_CALL(H5Dclose, dset_id)
    KRATOS_HDF5_CALL(H5Sclose, file_space_id)
    KRATOS_HDF5_CALL(H5Sclose, mem_space_id)

    KRATOS_INFO_IF("HDF5Application", GetEchoLevel() == 2)
        << "Read time \"" << rPath << "\": " << timer.ElapsedSeconds() << std::endl;

    KRATOS_CATCH("Path: \"" + rPath + "\".");
}

// template instantiations
template KRATOS_API(HDF5_APPLICATION) bool File::HasDataType<int>(const std::string&) const;
template KRATOS_API(HDF5_APPLICATION) bool File::HasDataType<double>(const std::string&) const;

#ifndef KRATOS_HDF5_FILE_DATA_SET_METHOD_INSTANTIATION
#define KRATOS_HDF5_FILE_DATA_SET_METHOD_INSTANTIATION(...)                                                                                                                                     \
    template KRATOS_API(HDF5_APPLICATION) void File::WriteDataSetImpl<__VA_ARGS__, File::DataTransferMode::Collective>(const std::string&, const __VA_ARGS__&, WriteInfo&);                     \
    template KRATOS_API(HDF5_APPLICATION) void File::ReadDataSetImpl<__VA_ARGS__, File::DataTransferMode::Collective>(const std::string&, __VA_ARGS__&, const unsigned, const unsigned) const;  \
    template KRATOS_API(HDF5_APPLICATION) void File::WriteDataSetImpl<__VA_ARGS__, File::DataTransferMode::Independent>(const std::string&, const __VA_ARGS__&, WriteInfo&);                    \
    template KRATOS_API(HDF5_APPLICATION) void File::ReadDataSetImpl<__VA_ARGS__, File::DataTransferMode::Independent>(const std::string&, __VA_ARGS__&, const unsigned, const unsigned) const; \
    template KRATOS_API(HDF5_APPLICATION) void File::WriteDataSet(const std::string&, const __VA_ARGS__&, WriteInfo&);                                                                          \
    template KRATOS_API(HDF5_APPLICATION) void File::WriteDataSetIndependent(const std::string&, const __VA_ARGS__&, WriteInfo&);                                                               \
    template KRATOS_API(HDF5_APPLICATION) void File::ReadDataSet(const std::string&, __VA_ARGS__&, const unsigned, const unsigned) const;                                                       \
    template KRATOS_API(HDF5_APPLICATION) void File::ReadDataSetIndependent(const std::string&, __VA_ARGS__&, const unsigned, const unsigned) const;                                            \

#endif

#ifndef KRATOS_HDF5_FILE_ATTRIBUTE_METHOD_INSTANTIATION
#define KRATOS_HDF5_FILE_ATTRIBUTE_METHOD_INSTANTIATION(...)                                                                        \
    template KRATOS_API(HDF5_APPLICATION) bool File::HasAttributeType<__VA_ARGS__>(const std::string&, const std::string&) const;   \
    template KRATOS_API(HDF5_APPLICATION) void File::WriteAttribute(const std::string&, const std::string&, const __VA_ARGS__&);    \
    template KRATOS_API(HDF5_APPLICATION) void File::ReadAttribute(const std::string&, const std::string&, __VA_ARGS__&) const;     \

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

KRATOS_HDF5_FILE_DATA_SET_METHOD_INSTANTIATION(Vector<char>);
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
#undef KRATOS_HDF5_CALL
#undef KRATOS_HDF5_CALL_WITH_RETURN

} // namespace HDF5.
} // namespace Kratos.
