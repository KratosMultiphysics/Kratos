#include "hdf5_file.h"

#include <algorithm>
#include <sstream>
#include <regex>
#include <utility>
#include "includes/kratos_parameters.h"

namespace Kratos
{
namespace HDF5
{
File::File(Parameters Settings)
{
    KRATOS_TRY;

    Parameters default_params(R"(
            {
                "file_name" : "PLEASE_SPECIFY_HDF5_FILENAME",
                "file_access_mode": "exclusive",
                "file_driver": "sec2",
                "echo_level" : 0
            })");

    Settings.RecursivelyValidateAndAssignDefaults(default_params);

    m_file_name = Settings["file_name"].GetString();
    KRATOS_ERROR_IF(m_file_name == "PLEASE_SPECIFY_HDF5_FILENAME") << "Invalid file name: " << m_file_name << std::endl;

    hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    std::string file_driver = Settings["file_driver"].GetString();
    SetFileDriver(file_driver, fapl_id);

    std::string file_access_mode = Settings["file_access_mode"].GetString();
    if (file_access_mode == "exclusive")
        m_file_id = H5Fcreate(m_file_name.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, fapl_id);
    else if (file_access_mode == "truncate")
        m_file_id = H5Fcreate(m_file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
    else
    {
        if (file_driver != "core")
            KRATOS_ERROR_IF(H5Fis_hdf5(m_file_name.c_str()) <= 0)
                << "Invalid HDF5 file: " << m_file_name << std::endl;
        if (file_access_mode == "read_only")
        {
            m_file_id = H5Fopen(m_file_name.c_str(), H5F_ACC_RDONLY, fapl_id);
        }
        else if (file_access_mode == "read_write")
            m_file_id = H5Fopen(m_file_name.c_str(), H5F_ACC_RDWR, fapl_id);
        else
            KRATOS_ERROR << "Invalid \"file_access_mode\": " << file_access_mode
                         << std::endl;
    }

    KRATOS_ERROR_IF(m_file_id < 0) << "Failed to open file: " << m_file_name << std::endl;

    KRATOS_ERROR_IF(H5Pclose(fapl_id) < 0) << "H5Pclose failed." << std::endl;

    m_echo_level = Settings["echo_level"].GetInt();

    KRATOS_CATCH("");
}

File::File(File&& rOther)
{
    m_file_name = std::move(rOther.m_file_name);
    m_file_id = rOther.m_file_id;
    rOther.m_file_id = -1;
    m_echo_level = rOther.m_echo_level;
    rOther.m_echo_level = 0;
}

File& File::operator=(File&& rOther)
{
    m_file_name = std::move(rOther.m_file_name);
    m_file_id = rOther.m_file_id;
    rOther.m_file_id = -1;
    m_echo_level = rOther.m_echo_level;
    rOther.m_echo_level = 0;
    return *this;
}

File::~File()
{
    H5Fclose(m_file_id);
}

bool File::HasPath(const std::string& rPath) const
{
    KRATOS_TRY;
    // Expects a valid path.
    KRATOS_ERROR_IF_NOT(Internals::IsPath(rPath)) << "Invalid path: \"" << rPath << '"' << std::endl;

    std::vector<std::string> splitted_path = Internals::Split(rPath, '/');
    std::string sub_path;
    for (const auto& r_link: splitted_path)
    {
        sub_path += '/' + r_link;

        htri_t link_found = H5Lexists(m_file_id, sub_path.c_str(), H5P_DEFAULT);
        KRATOS_ERROR_IF(link_found < 0) << "H5Lexists failed." << std::endl;
        if (!link_found)
            return false;

        htri_t object_found = H5Oexists_by_name(m_file_id, sub_path.c_str(), H5P_DEFAULT);
        KRATOS_ERROR_IF(object_found < 0) << "H5Oexists_by_name failed." << std::endl;
        if (!object_found)
            return false;
    }

    return true;
    KRATOS_CATCH("");
}

bool File::IsGroup(const std::string& rPath) const
{
    KRATOS_TRY;
    if (HasPath(rPath) == false) // Expects a valid path.
        return false;

    H5O_info_t object_info;
    KRATOS_ERROR_IF(H5Oget_info_by_name(m_file_id, rPath.c_str(), &object_info, H5P_DEFAULT) < 0)
        << "H5Oget_info_by_name failed." << std::endl;

    return (object_info.type == H5O_TYPE_GROUP);
    KRATOS_CATCH("");
}

bool File::IsDataSet(const std::string& rPath) const
{
    KRATOS_TRY;
    if (HasPath(rPath) == false) // Expects a valid path.
        return false;

    H5O_info_t object_info;
    KRATOS_ERROR_IF(H5Oget_info_by_name(m_file_id, rPath.c_str(), &object_info, H5P_DEFAULT) < 0)
        << "H5Oget_info_by_name failed." << std::endl;

    return (object_info.type == H5O_TYPE_DATASET);
    KRATOS_CATCH("");
}

bool File::HasAttribute(const std::string& rObjectPath, const std::string& rName) const
{
    KRATOS_TRY;
    htri_t status =
        H5Aexists_by_name(m_file_id, rObjectPath.c_str(), rName.c_str(), H5P_DEFAULT);
    KRATOS_ERROR_IF(status < 0) << "H5Aexists_by_name failed" << std::endl;
    return (status > 0);
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

    for (hsize_t i = 0; i < num_attrs; ++i)
    {
        // Get size of name.
        ssize_t ssize;
        ssize = H5Aget_name_by_idx(m_file_id, rObjectPath.c_str(), H5_INDEX_CRT_ORDER,
                                   H5_ITER_INC, i, buffer, max_ssize, H5P_DEFAULT);
        KRATOS_ERROR_IF(ssize < 0) << "H5Aget_name_by_idx failed." << std::endl;
        KRATOS_ERROR_IF(ssize > max_ssize) << "Attribute name size exceeds "
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
    hid_t group_id =
        H5Gcreate(m_file_id,rPath.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
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

    for (hsize_t i=0; i < num_links; ++i)
    {
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
    for (const auto& r_name : link_names)
        if (IsGroup(rGroupPath + '/' + r_name))
            names.push_back(r_name);
    return names;
    KRATOS_CATCH("");
}

void File::AddPath(const std::string& rPath)
{
    KRATOS_ERROR_IF_NOT(Internals::IsPath(rPath)) << "Invalid path: " << rPath << std::endl;

    std::vector<std::string> splitted_path = Internals::Split(rPath, '/');
    std::string sub_path;
    for (const auto& r_link: splitted_path)
    {
        sub_path += '/' + r_link;
        if (HasPath(sub_path) == false)
            CreateGroup(sub_path); // Add missing link.
        else
            KRATOS_ERROR_IF_NOT(IsGroup(sub_path))
                << "Path exists and is not a group: " << sub_path << std::endl;
    }
}

void File::WriteAttribute(const std::string& rObjectPath, const std::string& rName, const std::string& rValue)
{
    KRATOS_TRY;
    BuiltinTimer timer;
    hid_t type_id, space_id, attr_id;

    type_id = H5T_NATIVE_CHAR;
    const unsigned ndims = 1;
    hsize_t dims[ndims];
    dims[0] = rValue.size();
    space_id = H5Screate_simple(ndims, dims, nullptr);
    KRATOS_ERROR_IF(space_id < 0) << "H5Screate failed." << std::endl;
    attr_id = H5Acreate_by_name(m_file_id, rObjectPath.c_str(), rName.c_str(), type_id,
                                space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    KRATOS_ERROR_IF(attr_id < 0) << "H5Acreate_by_name failed." << std::endl;
    KRATOS_ERROR_IF(H5Awrite(attr_id, type_id, rValue.c_str()) < 0) << "H5Awrite failed." << std::endl;
    KRATOS_ERROR_IF(H5Sclose(space_id) < 0) << "H5Sclose failed." << std::endl;
    KRATOS_ERROR_IF(H5Aclose(attr_id) < 0) << "H5Aclose failed." << std::endl;
    if (GetEchoLevel() == 2 && GetPID() == 0)
        std::cout << "Write time \"" << rObjectPath << '/' << rName << "\": " << timer.ElapsedSeconds() << std::endl;
    KRATOS_CATCH("Path: \"" + rObjectPath + '/' + rName + "\".");
}

void File::WriteAttribute(const std::string& rObjectPath, const std::string& rName, const array_1d<double, 3>& rValue)
{
    KRATOS_TRY;
    Vector<double> vector_value = rValue;
    WriteAttribute(rObjectPath, rName, vector_value);
    KRATOS_CATCH("Path: \"" + rObjectPath + '/' + rName + "\".");
}

void File::WriteDataSet(const std::string& rPath, const Vector<int>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the base class method. Please override in the derived class." << std::endl;
    KRATOS_CATCH("");
}

void File::WriteDataSet(const std::string& rPath, const Vector<double>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the base class method. Please override in the derived class." << std::endl;
    KRATOS_CATCH("");
}

void File::WriteDataSet(const std::string& rPath, const Vector<array_1d<double, 3>>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the base class method. Please override in the derived class." << std::endl;
    KRATOS_CATCH("");
}

void File::WriteDataSet(const std::string& rPath, const Matrix<int>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the base class method. Please override in the derived class." << std::endl;
    KRATOS_CATCH("");
}

void File::WriteDataSet(const std::string& rPath, const Matrix<double>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the base class method. Please override in the derived class." << std::endl;
    KRATOS_CATCH("");
}

void File::WriteDataSetIndependent(const std::string& rPath, const Vector<int>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the base class method. Please override in the derived class." << std::endl;
    KRATOS_CATCH("");
}

void File::WriteDataSetIndependent(const std::string& rPath, const Vector<double>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the base class method. Please override in the derived class." << std::endl;
    KRATOS_CATCH("");
}

void File::WriteDataSetIndependent(const std::string& rPath, const Vector<array_1d<double, 3>>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the base class method. Please override in the derived class." << std::endl;
    KRATOS_CATCH("");
}

void File::WriteDataSetIndependent(const std::string& rPath, const Matrix<int>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the base class method. Please override in the derived class." << std::endl;
    KRATOS_CATCH("");
}

void File::WriteDataSetIndependent(const std::string& rPath, const Matrix<double>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the base class method. Please override in the derived class." << std::endl;
    KRATOS_CATCH("");
}

std::vector<unsigned> File::GetDataDimensions(const std::string& rPath) const
{
    KRATOS_TRY;
    constexpr int max_ndims = 5;
    int ndims;
    hsize_t dims[max_ndims];
    hid_t dset_id, dspace_id;
    KRATOS_ERROR_IF((dset_id = H5Dopen(m_file_id, rPath.c_str(), H5P_DEFAULT)) < 0)
        << "H5Dopen failed." << std::endl;
    KRATOS_ERROR_IF((dspace_id = H5Dget_space(dset_id)) < 0)
        << "H5Dget_space failed." << std::endl;
    KRATOS_ERROR_IF((ndims = H5Sget_simple_extent_ndims(dspace_id)) < 0)
        << "H5Sget_simple_extent_ndims failed." << std::endl;
    KRATOS_ERROR_IF(max_ndims < ndims) << "Maximum dimension: " << max_ndims
                                       << ", dimension: " << ndims << std::endl;
    KRATOS_ERROR_IF(H5Sget_simple_extent_dims(dspace_id, dims, nullptr) < 0)
        << "H5Sget_simple_extent_dims failed" << std::endl;
    KRATOS_ERROR_IF(H5Sclose(dspace_id) < 0) << "H5Sclose failed." << std::endl;
    KRATOS_ERROR_IF(H5Dclose(dset_id) < 0) << "H5Dclose failed." << std::endl;

    return std::vector<unsigned>(dims, dims + ndims);
    KRATOS_CATCH("");
}

bool File::HasIntDataType(const std::string& rPath) const
{
    KRATOS_TRY;
    hid_t dset_id, dtype_id;
    KRATOS_ERROR_IF((dset_id = H5Dopen(m_file_id, rPath.c_str(), H5P_DEFAULT)) < 0)
        << "H5Dopen failed." << std::endl;
    KRATOS_ERROR_IF((dtype_id = H5Dget_type(dset_id)) < 0)
        << "H5Dget_type failed." << std::endl;
    H5T_class_t type = H5Tget_class(dtype_id);
    KRATOS_ERROR_IF(type == H5T_NO_CLASS) << "Invalid data type." << std::endl;
    KRATOS_ERROR_IF(H5Tclose(dtype_id) < 0) << "H5Tclose failed." << std::endl;
    KRATOS_ERROR_IF(H5Dclose(dset_id) < 0) << "H5Dclose failed." << std::endl;

    return (type == H5T_INTEGER);
    KRATOS_CATCH("");
}

bool File::HasFloatDataType(const std::string& rPath) const
{
    KRATOS_TRY;
    hid_t dset_id, dtype_id;
    KRATOS_ERROR_IF((dset_id = H5Dopen(m_file_id, rPath.c_str(), H5P_DEFAULT)) < 0)
        << "H5Dopen failed." << std::endl;
    KRATOS_ERROR_IF((dtype_id = H5Dget_type(dset_id)) < 0)
        << "H5Dget_type failed." << std::endl;
    H5T_class_t type = H5Tget_class(dtype_id);
    KRATOS_ERROR_IF(type == H5T_NO_CLASS) << "Invalid data type." << std::endl;
    KRATOS_ERROR_IF(H5Tclose(dtype_id) < 0) << "H5Tclose failed." << std::endl;
    KRATOS_ERROR_IF(H5Dclose(dset_id) < 0) << "H5Dclose failed." << std::endl;

    return (type == H5T_FLOAT);
    KRATOS_CATCH("");
}

void File::Flush()
{
    KRATOS_ERROR_IF(H5Fflush(m_file_id, H5F_SCOPE_GLOBAL) < 0)
        << "H5Fflush failed." << std::endl;
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

unsigned File::GetPID() const
{
    return 0;
}

unsigned File::GetTotalProcesses() const
{
    return 1;
}

void File::ReadAttribute(const std::string& rObjectPath, const std::string& rName, std::string& rValue)
{
    KRATOS_TRY;
    BuiltinTimer timer;
    hid_t mem_type_id, attr_type_id, space_id, attr_id;
    int ndims;
    hsize_t dims[2];
    const unsigned max_ssize = 100;
    char buffer[max_ssize];

    mem_type_id = H5T_NATIVE_CHAR;
    attr_id = H5Aopen_by_name(m_file_id, rObjectPath.c_str(), rName.c_str(),
                                    H5P_DEFAULT, H5P_DEFAULT);
    KRATOS_ERROR_IF(attr_id < 0) << "H5Aopen_by_name failed." << std::endl;

    // Check data type.
    attr_type_id = H5Aget_type(attr_id);
    KRATOS_ERROR_IF(attr_type_id < 0) << "H5Aget_type failed." << std::endl;
    htri_t is_valid_type = H5Tequal(mem_type_id, attr_type_id);
    KRATOS_ERROR_IF(H5Tclose(attr_type_id) < 0) << "H5Tclose failed." << std::endl;
    KRATOS_ERROR_IF(is_valid_type < 0) << "H5Tequal failed." << std::endl;
    KRATOS_ERROR_IF(is_valid_type == 0) << "Attribute \"" << rName << "\" is not a string." << std::endl;

    // Check dimensions.
    space_id = H5Aget_space(attr_id);
    KRATOS_ERROR_IF(space_id < 0) << "H5Aget_space failed." << std::endl;
    KRATOS_ERROR_IF((ndims = H5Sget_simple_extent_ndims(space_id)) < 0)
        << "H5Sget_simple_extent_ndims failed." << std::endl;
    KRATOS_ERROR_IF(ndims != 1) << "Attribute \"" << rName << "\" is not string." << std::endl;
    KRATOS_ERROR_IF(H5Sget_simple_extent_dims(space_id, dims, nullptr) < 0)
        << "H5Sget_simple_extent_dims failed" << std::endl;
    KRATOS_ERROR_IF(H5Sclose(space_id) < 0) << "H5Sclose failed." << std::endl;
    KRATOS_ERROR_IF(max_ssize < dims[0]) << "String size is greater than " << max_ssize << '.' << std::endl;
    // Read attribute.
    KRATOS_ERROR_IF(H5Aread(attr_id, mem_type_id, buffer) < 0) << "H5Aread failed." << std::endl;
    KRATOS_ERROR_IF(H5Aclose(attr_id) < 0) << "H5Aclose failed." << std::endl;
    rValue = std::string(buffer, dims[0]);
    if (GetEchoLevel() == 2 && GetPID() == 0)
        std::cout << "Read time \"" << rObjectPath << '/' << rName << "\": " << timer.ElapsedSeconds() << std::endl;
    KRATOS_CATCH("Path: \"" + rObjectPath + '/' + rName + "\".");
}

void File::ReadAttribute(const std::string& rObjectPath, const std::string& rName, array_1d<double, 3>& rValue)
{
    KRATOS_TRY;
    Vector<double> vector_value;
    ReadAttribute(rObjectPath, rName, vector_value);
    KRATOS_ERROR_IF(vector_value.size() > 3)
        << "Invalid size (" << vector_value.size() << ") for array_1d!" << std::endl;
    rValue = vector_value;
    KRATOS_CATCH("Path: \"" + rObjectPath + '/' + rName + "\".");
}

void File::ReadDataSet(const std::string& rPath, Vector<int>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the base class method. Please override in the derived class." << std::endl;
    KRATOS_CATCH("");
}

void File::ReadDataSet(const std::string& rPath, Vector<double>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the base class method. Please override in the derived class." << std::endl;
    KRATOS_CATCH("");
}

void File::ReadDataSet(const std::string& rPath, Vector<array_1d<double, 3>>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the base class method. Please override in the derived class." << std::endl;
    KRATOS_CATCH("");
}

void File::ReadDataSet(const std::string& rPath, Matrix<int>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the base class method. Please override in the derived class." << std::endl;
    KRATOS_CATCH("");
}

void File::ReadDataSet(const std::string& rPath, Matrix<double>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the base class method. Please override in the derived class." << std::endl;
    KRATOS_CATCH("");
}

void File::ReadDataSetIndependent(const std::string& rPath, Vector<int>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the base class method. Please override in the derived class." << std::endl;
    KRATOS_CATCH("");
}

void File::ReadDataSetIndependent(const std::string& rPath, Vector<double>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the base class method. Please override in the derived class." << std::endl;
    KRATOS_CATCH("");
}

void File::ReadDataSetIndependent(const std::string& rPath, Vector<array_1d<double, 3>>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the base class method. Please override in the derived class." << std::endl;
    KRATOS_CATCH("");
}

void File::ReadDataSetIndependent(const std::string& rPath, Matrix<int>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the base class method. Please override in the derived class." << std::endl;
    KRATOS_CATCH("");
}

void File::ReadDataSetIndependent(const std::string& rPath, Matrix<double>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the base class method. Please override in the derived class." << std::endl;
    KRATOS_CATCH("");
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

namespace Internals
{
bool IsPath(const std::string& rPath)
{
#if defined(__GNUC__) && !defined(__clang__) && (__GNUC__ < 4 || (__GNUC__ == 4 && (__GNUC_MINOR__ < 9)))
    KRATOS_ERROR << "This method is not compiled well. You should use a GCC 4.9 or higher" << std::endl;
#else
    return regex_match(rPath, std::regex("(/[\\w\\(\\)]+)+"));
#endif
}

std::vector<std::string> Split(const std::string& rPath, char Delimiter)
{
    std::vector<std::string> splitted;
    splitted.reserve(10);
    std::stringstream ss(rPath);
    std::string sub_string;
    while (std::getline(ss, sub_string, Delimiter))
        if (sub_string.size() > 0)
            splitted.push_back(sub_string);
    return splitted;
}
} // namespace Internals.

} // namespace HDF5.
} // namespace Kratos.
