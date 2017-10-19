#include "hdf5_file.h"

#include <cassert>

namespace Kratos
{
HDF5File::HDF5File(Parameters& rParams)
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

    std::string file_access_mode = rParams["file_access_mode"].GetString();
    hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);

#if (defined(_WIN32) || defined(_WIN64))
    KRATOS_ERROR_IF(H5Pset_fapl_windows(fapl_id) < 0)
        << "H5Pset_fapl_windows failed." << std::endl;
#else
    std::string file_driver = rParams["file_driver"].GetString();
    if (file_driver == "sec2")
    {
        KRATOS_ERROR_IF(H5Pset_fapl_sec2(fapl_id) < 0)
            << "H5Pset_fapl_sec2 failed." << std::endl;
    }
    else if (file_driver == "stdio")
    {
        KRATOS_ERROR_IF(H5Pset_fapl_stdio(fapl_id) < 0)
            << "H5Pset_fapl_stdio failed." << std::endl;
    }
    else if (file_driver == "core")
    {
        KRATOS_ERROR_IF(H5Pset_fapl_core(fapl_id, 1000000, 0) < 0)
            << "H5Pset_fapl_core failed." << std::endl;
    }
    else
        KRATOS_ERROR << "Unsupported \"file_driver\": " << file_driver << std::endl;
#endif

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

HDF5File::~HDF5File()
{
#if !defined(NDEBUG)
    // Check for open data sets, groups, datatypes, attributes before
    // closing the file.
    ssize_t num_open_objects = H5Fget_obj_count(m_file_id, H5F_OBJ_ALL);
    assert(num_open_objects == 1); // 1 file object only
#endif
    H5Fclose(m_file_id);
}

bool HDF5File::HasPath(std::string Path) const
{
    KRATOS_TRY;

    KRATOS_ERROR_IF_NOT(HDF5Utils::IsPath(Path)) << "Invalid path: " << Path << std::endl;

    std::string sub_path;
    decltype(Path.size()) pos = 0; /* Make sure pos has same size as
        std::string::npos. */
    while (pos < Path.size()) // Check each link in the path.
    {
        pos = Path.find('/', ++pos);
        if (pos != std::string::npos)
            sub_path = Path.substr(0, pos); // Check current subpath.
        else
        {
            sub_path = Path;   // Check the complete path.
            pos = Path.size(); // Exit on loop completion.
        }

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

bool HDF5File::IsGroup(std::string Path) const
{
    KRATOS_TRY;
    if (HasPath(Path) == false)
        return false;

    H5O_info_t object_info;
    KRATOS_ERROR_IF(H5Oget_info_by_name(m_file_id, Path.c_str(), &object_info, H5P_DEFAULT) < 0)
        << "H5Oget_info_by_name failed." << std::endl;

    return (object_info.type == H5O_TYPE_GROUP);
    KRATOS_CATCH("");
}

bool HDF5File::IsDataSet(std::string Path) const
{
    KRATOS_TRY;
    if (HasPath(Path) == false)
        return false;

    H5O_info_t object_info;
    KRATOS_ERROR_IF(H5Oget_info_by_name(m_file_id, Path.c_str(), &object_info, H5P_DEFAULT) < 0)
        << "H5Oget_info_by_name failed." << std::endl;

    return (object_info.type == H5O_TYPE_DATASET);
    KRATOS_CATCH("");
}

void HDF5File::CreateGroup(std::string Path)
{
    KRATOS_TRY;
    hid_t group_id =
        H5Gcreate(m_file_id, Path.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    KRATOS_ERROR_IF(group_id < 0) << "H5Gcreate failed." << std::endl;
    KRATOS_ERROR_IF(H5Gclose(group_id) < 0) << "H5Gclose failed." << std::endl;
    KRATOS_CATCH("");
}

void HDF5File::AddPath(std::string Path)
{
    KRATOS_ERROR_IF(HDF5Utils::IsPath(Path) == false) << "Invalid path: " << Path << std::endl;

    std::vector<std::string> splitted_path = HDF5Utils::Split(Path, '/');
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

void HDF5File::WriteDataSet(std::string Path, const std::vector<int>& rData)
{
    KRATOS_TRY;
    WriteDataSetImpl(Path, rData);
    KRATOS_CATCH("");
}

void HDF5File::WriteDataSet(std::string Path, const std::vector<double>& rData)
{
    KRATOS_TRY;
    WriteDataSetImpl(Path, rData);
    KRATOS_CATCH("");
}

void HDF5File::WriteDataSet(std::string Path, const std::vector<array_1d<double, 3>>& rData)
{
    KRATOS_TRY;
    WriteDataSetImpl(Path, rData);
    KRATOS_CATCH("");
}

std::vector<unsigned> HDF5File::GetDataDimensions(std::string Path) const
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

bool HDF5File::HasIntDataType(std::string Path) const
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

bool HDF5File::HasFloatDataType(std::string Path) const
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

void HDF5File::Flush()
{
    KRATOS_ERROR_IF(H5Fflush(m_file_id, H5F_SCOPE_GLOBAL) < 0)
        << "H5Fflush failed." << std::endl;
}

unsigned HDF5File::GetFileSize() const
{
    hsize_t size;
    KRATOS_ERROR_IF(H5Fget_filesize(m_file_id, &size) < 0)
        << "H5Fget_filesize failed." << std::endl;

    return size;
}

std::string HDF5File::GetFileName() const
{
    return m_file_name;
}

void HDF5File::ReadDataSet(std::string Path, std::vector<int>& rData, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetImpl(Path, rData, BlockSize);
    KRATOS_CATCH("");
}

void HDF5File::ReadDataSet(std::string Path, std::vector<double>& rData, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetImpl(Path, rData, BlockSize);
    KRATOS_CATCH("");
}

void HDF5File::ReadDataSet(std::string Path, std::vector<array_1d<double, 3>>& rData, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetImpl(Path, rData, BlockSize);
    KRATOS_CATCH("");
}

} // // namespace Kratos.