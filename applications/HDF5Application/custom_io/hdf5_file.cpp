#include "hdf5_file.h"

#include <algorithm>
#include <sstream>
#include <regex>

namespace Kratos
{
namespace HDF5
{
File::File(Parameters& rParams)
{
    KRATOS_TRY;

    Parameters default_params(R"(
            {
                "file_name" : "PLEASE_SPECIFY_HDF5_FILENAME",
                "file_access_mode": "exclusive",
                "file_driver": "sec2",
                "echo_level" : 0
            })");

    rParams.RecursivelyValidateAndAssignDefaults(default_params);

    m_file_name = rParams["file_name"].GetString();

    hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    std::string file_driver = rParams["file_driver"].GetString();
    SetFileDriver(file_driver, fapl_id);

    std::string file_access_mode = rParams["file_access_mode"].GetString();
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

    m_echo_level = rParams["echo_level"].GetInt();

    KRATOS_CATCH("");
}

File::~File()
{
    H5Fclose(m_file_id);
}

bool File::HasPath(std::string Path) const
{
    KRATOS_TRY;
    // Expects a valid path.
    KRATOS_ERROR_IF_NOT(Internals::IsPath(Path)) << "Invalid path: \"" << Path << '"' << std::endl;

    std::vector<std::string> splitted_path = Internals::Split(Path, '/');
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

bool File::IsGroup(std::string Path) const
{
    KRATOS_TRY;
    if (HasPath(Path) == false) // Expects a valid path.
        return false;

    H5O_info_t object_info;
    KRATOS_ERROR_IF(H5Oget_info_by_name(m_file_id, Path.c_str(), &object_info, H5P_DEFAULT) < 0)
        << "H5Oget_info_by_name failed." << std::endl;

    return (object_info.type == H5O_TYPE_GROUP);
    KRATOS_CATCH("");
}

bool File::IsDataSet(std::string Path) const
{
    KRATOS_TRY;
    if (HasPath(Path) == false) // Expects a valid path.
        return false;

    H5O_info_t object_info;
    KRATOS_ERROR_IF(H5Oget_info_by_name(m_file_id, Path.c_str(), &object_info, H5P_DEFAULT) < 0)
        << "H5Oget_info_by_name failed." << std::endl;

    return (object_info.type == H5O_TYPE_DATASET);
    KRATOS_CATCH("");
}

bool File::HasAttribute(std::string ObjectPath, std::string Name) const
{
    KRATOS_TRY;
    htri_t status =
        H5Aexists_by_name(m_file_id, ObjectPath.c_str(), Name.c_str(), H5P_DEFAULT);
    KRATOS_ERROR_IF(status < 0) << "H5Aexists_by_name failed" << std::endl;
    return (status > 0);
    KRATOS_CATCH("");
}

void File::GetAttributeNames(std::string ObjectPath, std::vector<std::string>& rNames) const
{
    KRATOS_TRY;
    constexpr unsigned max_ssize = 100;
    char buffer[max_ssize];
    // Get number of attributes.
    hid_t object_id = H5Oopen(m_file_id, ObjectPath.c_str(), H5P_DEFAULT);
    KRATOS_ERROR_IF(object_id < 0) << "H5Oopen failed." << std::endl;
    H5O_info_t object_info;
    KRATOS_ERROR_IF(H5Oget_info(object_id, &object_info) < 0)
        << "H5Oget_info failed." << std::endl;
    hsize_t num_attrs = object_info.num_attrs;
    rNames.resize(num_attrs);

    for (hsize_t i = 0; i < num_attrs; ++i)
    {
        // Get size of name.
        ssize_t ssize;
        ssize = H5Aget_name_by_idx(m_file_id, ObjectPath.c_str(), H5_INDEX_CRT_ORDER,
                                   H5_ITER_INC, i, buffer, max_ssize, H5P_DEFAULT);
        KRATOS_ERROR_IF(ssize < 0) << "H5Aget_name_by_idx failed." << std::endl;
        KRATOS_ERROR_IF(ssize > max_ssize) << "Attribute name size exceeds "
                                           << max_ssize << std::endl;
        rNames[i].resize(ssize);
        std::copy_n(buffer, ssize, rNames[i].begin());
    }
    KRATOS_ERROR_IF(H5Oclose(object_id) < 0) << "H5Oclose failed." << std::endl;
    KRATOS_CATCH("");
}

void File::CreateGroup(std::string Path)
{
    KRATOS_TRY;
    hid_t group_id =
        H5Gcreate(m_file_id, Path.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    KRATOS_ERROR_IF(group_id < 0) << "H5Gcreate failed." << std::endl;
    KRATOS_ERROR_IF(H5Gclose(group_id) < 0) << "H5Gclose failed." << std::endl;
    KRATOS_CATCH("");
}

void File::AddPath(std::string Path)
{
    KRATOS_ERROR_IF(Internals::IsPath(Path) == false) << "Invalid path: " << Path << std::endl;

    std::vector<std::string> splitted_path = Internals::Split(Path, '/');
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

void File::WriteDataSet(std::string Path, const Vector<int>& rData)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the base class method. Please override in the derived class." << std::endl;
    KRATOS_CATCH("");
}

void File::WriteDataSet(std::string Path, const Vector<double>& rData)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the base class method. Please override in the derived class." << std::endl;
    KRATOS_CATCH("");
}

void File::WriteDataSet(std::string Path, const Vector<array_1d<double, 3>>& rData)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the base class method. Please override in the derived class." << std::endl;
    KRATOS_CATCH("");
}

void File::WriteDataSet(std::string Path, const Matrix<int>& rData)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the base class method. Please override in the derived class." << std::endl;
    KRATOS_CATCH("");
}

void File::WriteDataSet(std::string Path, const Matrix<double>& rData)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the base class method. Please override in the derived class." << std::endl;
    KRATOS_CATCH("");
}

void File::WriteDataPartition(std::string Path, const Vector<int>& rData)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the base class method. Please override in the derived class." << std::endl;
    KRATOS_CATCH("");
}

void File::WriteDataPartition(std::string Path, const Vector<double>& rData)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the base class method. Please override in the derived class." << std::endl;
    KRATOS_CATCH("");
}

void File::WriteDataPartition(std::string Path, const Vector<array_1d<double,3>>& rData)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the base class method. Please override in the derived class." << std::endl;
    KRATOS_CATCH("");
}

void File::WriteDataPartition(std::string Path, const Matrix<int>& rData)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the base class method. Please override in the derived class." << std::endl;
    KRATOS_CATCH("");
}

void File::WriteDataPartition(std::string Path, const Matrix<double>& rData)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the base class method. Please override in the derived class." << std::endl;
    KRATOS_CATCH("");
}

void File::WriteDataSetIndependent(std::string Path, const Vector<int>& rData)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the base class method. Please override in the derived class." << std::endl;
    KRATOS_CATCH("");
}

void File::WriteDataSetIndependent(std::string Path, const Vector<double>& rData)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the base class method. Please override in the derived class." << std::endl;
    KRATOS_CATCH("");
}

void File::WriteDataSetIndependent(std::string Path, const Vector<array_1d<double, 3>>& rData)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the base class method. Please override in the derived class." << std::endl;
    KRATOS_CATCH("");
}

void File::WriteDataSetIndependent(std::string Path, const Matrix<int>& rData)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the base class method. Please override in the derived class." << std::endl;
    KRATOS_CATCH("");
}

void File::WriteDataSetIndependent(std::string Path, const Matrix<double>& rData)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the base class method. Please override in the derived class." << std::endl;
    KRATOS_CATCH("");
}

std::vector<unsigned> File::GetDataDimensions(std::string Path) const
{
    KRATOS_TRY;
    constexpr int max_ndims = 5;
    int ndims;
    hsize_t dims[max_ndims];
    hid_t dset_id, dspace_id;
    KRATOS_ERROR_IF((dset_id = H5Dopen(m_file_id, Path.c_str(), H5P_DEFAULT)) < 0)
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

bool File::HasIntDataType(std::string Path) const
{
    KRATOS_TRY;
    hid_t dset_id, dtype_id;
    KRATOS_ERROR_IF((dset_id = H5Dopen(m_file_id, Path.c_str(), H5P_DEFAULT)) < 0)
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

bool File::HasFloatDataType(std::string Path) const
{
    KRATOS_TRY;
    hid_t dset_id, dtype_id;
    KRATOS_ERROR_IF((dset_id = H5Dopen(m_file_id, Path.c_str(), H5P_DEFAULT)) < 0)
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

void File::ReadDataSet(std::string Path, Vector<int>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the base class method. Please override in the derived class." << std::endl;
    KRATOS_CATCH("");
}

void File::ReadDataSet(std::string Path, Vector<double>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the base class method. Please override in the derived class." << std::endl;
    KRATOS_CATCH("");
}

void File::ReadDataSet(std::string Path, Vector<array_1d<double, 3>>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the base class method. Please override in the derived class." << std::endl;
    KRATOS_CATCH("");
}

void File::ReadDataSet(std::string Path, Matrix<int>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the base class method. Please override in the derived class." << std::endl;
    KRATOS_CATCH("");
}

void File::ReadDataSet(std::string Path, Matrix<double>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the base class method. Please override in the derived class." << std::endl;
    KRATOS_CATCH("");
}

void File::ReadDataSetIndependent(std::string Path, Vector<int>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the base class method. Please override in the derived class." << std::endl;
    KRATOS_CATCH("");
}

void File::ReadDataSetIndependent(std::string Path, Vector<double>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the base class method. Please override in the derived class." << std::endl;
    KRATOS_CATCH("");
}

void File::ReadDataSetIndependent(std::string Path, Vector<array_1d<double, 3>>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the base class method. Please override in the derived class." << std::endl;
    KRATOS_CATCH("");
}

void File::ReadDataSetIndependent(std::string Path, Matrix<int>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the base class method. Please override in the derived class." << std::endl;
    KRATOS_CATCH("");
}

void File::ReadDataSetIndependent(std::string Path, Matrix<double>& rData, unsigned StartIndex, unsigned BlockSize)
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
bool IsPath(std::string Path)
{
    return regex_match(Path, std::regex("(/[\\w\\(\\)]+)+"));
}

std::vector<std::string> Split(std::string Path, char Delimiter)
{
    std::vector<std::string> result;
    result.reserve(10);
    std::stringstream ss(Path);
    std::string sub_string;
    while (std::getline(ss, sub_string, Delimiter))
        if (sub_string.size() > 0)
            result.push_back(sub_string);

    return result;
}
} // namespace Internals.

} // namespace HDF5.
} // namespace Kratos.