#include "hdf5_file.h"

#include <regex>
#include <cassert>

namespace Kratos
{

    HDF5File::HDF5File(Parameters& rParams)
    {
        KRATOS_TRY;

        herr_t status;

        Parameters default_params(R"(
            {
                "file_name" : "PLEASE_SPECIFY_FILENAME",
                "file_access_mode": "exclusive",
                "file_driver": "sec2",
                "echo_level" : 0
            })");

        rParams.RecursivelyValidateAndAssignDefaults(default_params);

        m_file_name = rParams["file_name"].GetString();

        std::string file_access_mode = rParams["file_access_mode"].GetString();
        hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);

#if (defined(_WIN32) || defined(_WIN64))
        status = H5Pset_fapl_windows(fapl_id);
#else
        std::string file_driver = rParams["file_driver"].GetString();
        if (file_driver == "sec2")
            status = H5Pset_fapl_sec2(fapl_id); // default posix sec 2
        else if (file_driver == "stdio")
            status = H5Pset_fapl_stdio(fapl_id);
        else if (file_driver == "core")
            status = H5Pset_fapl_core(fapl_id, 1000000, 0); // Use for testing so file is not written to disk.
        else
            KRATOS_ERROR << "Unsupported \"file_driver\": " << file_driver << std::endl;
#endif
        KRATOS_ERROR_IF(status < 0) << "Failed to set file driver." << std::endl;

        if (file_access_mode == "exclusive")
            m_file_id = H5Fcreate(m_file_name.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, fapl_id);
        else if (file_access_mode == "truncate")
            m_file_id = H5Fcreate(m_file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
        else
        {
            if (file_driver != "core")
                KRATOS_ERROR_IF(H5Fis_hdf5(m_file_name.c_str()) <= 0) << "Invalid HDF5 file: " << m_file_name << std::endl;
            if (file_access_mode == "read_only")
            {
                m_file_id = H5Fopen(m_file_name.c_str(), H5F_ACC_RDONLY, fapl_id);
            }
            else if (file_access_mode == "read_write")
                m_file_id = H5Fopen(m_file_name.c_str(), H5F_ACC_RDWR, fapl_id);
            else
                KRATOS_ERROR << "Invalid \"file_access_mode\": " << file_access_mode << std::endl;
        }
        
        KRATOS_ERROR_IF(m_file_id < 0) << "Failed to open file: " << m_file_name << std::endl;

        status = H5Pclose(fapl_id);

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

    /// Check if string is a valid path.
    /**
     * Valid paths are similar to linux file system with alphanumeric names
     * and possible underscores separated by '/'. All paths are absolute.
     */
    bool HDF5File::IsPath(std::string rPath) const
    {
        std::regex pattern("(/\\w+)+"); 
        return regex_match(rPath, pattern);
    }

    /// Check if path exists in HDF5 file.
    bool HDF5File::HasPath(std::string rPath) const
    {
        KRATOS_ERROR_IF_NOT(IsPath(rPath)) << "Invalid path: " << rPath << std::endl;
        
        std::string sub_path;
        decltype(rPath.size()) pos = 0; // Make sure pos is large enough to store std::string::npos.
        // Check each link in the path.
        while (pos < rPath.size())
        {
            pos = rPath.find('/', ++pos);
            if (pos != std::string::npos)
                sub_path = rPath.substr(0, pos); // Check current subpath.
            else
            {
                sub_path = rPath; // Check the complete path.
                pos = rPath.size(); // Exit on loop completion.
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
    }

    bool HDF5File::IsGroup(std::string rPath) const
    {
        if (HasPath(rPath) == false)
            return false;

        H5O_info_t object_info;
        herr_t status = H5Oget_info_by_name(m_file_id, rPath.c_str(), &object_info, H5P_DEFAULT);
        KRATOS_ERROR_IF(status < 0) << "H5Oget_info_by_name failed." << std::endl;

        return (object_info.type == H5O_TYPE_GROUP);
    }

    bool HDF5File::IsDataSet(std::string rPath) const
    {
        if (HasPath(rPath) == false)
            return false;

        H5O_info_t object_info;
        herr_t status = H5Oget_info_by_name(m_file_id, rPath.c_str(), &object_info, H5P_DEFAULT);
        KRATOS_ERROR_IF(status < 0) << "H5Oget_info_by_name failed." << std::endl;

        return (object_info.type == H5O_TYPE_DATASET);
    }

    void HDF5File::CreateGroup(std::string rPath)
    {
        hid_t group_id = H5Gcreate(m_file_id, rPath.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        KRATOS_ERROR_IF(group_id < 0) << "H5Gcreate failed." << std::endl;
        herr_t status = H5Gclose(group_id);
        KRATOS_ERROR_IF(status < 0) << "H5Gclose failed." << std::endl;
    }

    void HDF5File::WriteDataSet(std::string rPath, const std::vector<int>& rData)
    {
        WriteDataSetImpl(rPath, rData);
    }
    
    void HDF5File::WriteDataSet(std::string rPath, const std::vector<double>& rData)
    {
        WriteDataSetImpl(rPath, rData);        
    }
    
    void HDF5File::WriteDataSet(std::string rPath, const std::vector<array_1d<double,3>>& rData)
    {
        WriteDataSetImpl(rPath, rData);        
    }

    std::vector<unsigned> HDF5File::GetDataDimensions(std::string rPath)
    {
        KRATOS_ERROR_IF_NOT(IsDataSet(rPath)) << "Invalid path: " << rPath << std::endl;
        
        hsize_t* dims = nullptr;
        hsize_t* maxdims = nullptr;
        herr_t status;
        hid_t dset_id = H5Dopen(m_file_id, rPath.c_str(), H5P_DEFAULT);
        KRATOS_ERROR_IF(dset_id < 0) << "H5Dopen failed." << std::endl;
        hid_t dspace_id = H5Dget_space(dset_id);
        KRATOS_ERROR_IF(dspace_id < 0) << "H5Dget_space failed." << std::endl;
        int ndims = H5Sget_simple_extent_dims(dspace_id, dims, maxdims);
        KRATOS_ERROR_IF(ndims < 0) << "H5Sget_simple_extent_dims failed." << std::endl;
        status = H5Sclose(dspace_id);
        KRATOS_ERROR_IF(status < 0) << "H5Sclose failed." << std::endl;
        status = H5Dclose(dset_id);
        KRATOS_ERROR_IF(status < 0) << "H5Dclose failed." << std::endl;
        
        return std::vector<unsigned>(dims, dims + ndims);
    }
    
    bool HDF5File::HasIntDataType(std::string rPath)
    {
        KRATOS_ERROR_IF_NOT(IsDataSet(rPath)) << "Invalid path: " << rPath << std::endl;
        
        herr_t status;
        hid_t dset_id = H5Dopen(m_file_id, rPath.c_str(), H5P_DEFAULT);
        KRATOS_ERROR_IF(dset_id < 0) << "H5Dopen failed." << std::endl;
        hid_t dtype_id = H5Dget_type(dset_id);
        KRATOS_ERROR_IF(dtype_id < 0) << "H5Dget_type failed." << std::endl;
        H5T_class_t type = H5Tget_class(dtype_id);
        KRATOS_ERROR_IF(type == H5T_NO_CLASS) << "Invalid data type." << std::endl;
        status = H5Tclose(dtype_id);
        KRATOS_ERROR_IF(status < 0) << "H5Tclose failed." << std::endl;
        status = H5Dclose(dset_id);
        KRATOS_ERROR_IF(status < 0) << "H5Dclose failed." << std::endl;

        return (type == H5T_INTEGER);
    }
    
    bool HDF5File::HasFloatDataType(std::string rPath)
    {
        KRATOS_ERROR_IF_NOT(IsDataSet(rPath)) << "Invalid path: " << rPath << std::endl;
        
        herr_t status;
        hid_t dset_id = H5Dopen(m_file_id, rPath.c_str(), H5P_DEFAULT);
        KRATOS_ERROR_IF(dset_id < 0) << "H5Dopen failed." << std::endl;
        hid_t dtype_id = H5Dget_type(dset_id);
        KRATOS_ERROR_IF(dtype_id < 0) << "H5Dget_type failed." << std::endl;
        H5T_class_t type = H5Tget_class(dtype_id);
        KRATOS_ERROR_IF(type == H5T_NO_CLASS) << "Invalid data type." << std::endl;
        status = H5Tclose(dtype_id);
        KRATOS_ERROR_IF(status < 0) << "H5Tclose failed." << std::endl;
        status = H5Dclose(dset_id);
        KRATOS_ERROR_IF(status < 0) << "H5Dclose failed." << std::endl;

        return (type == H5T_FLOAT);
    }

    void HDF5File::Flush()
    {
        herr_t status = H5Fflush(m_file_id, H5F_SCOPE_GLOBAL);
        KRATOS_ERROR_IF(status < 0) << "H5Fflush failed." << std::endl;
    }

    unsigned HDF5File::GetFileSize() const
    {
        hsize_t size;
        herr_t status = H5Fget_filesize(m_file_id, &size);
        KRATOS_ERROR_IF(status < 0) << "H5Fget_filesize failed." << std::endl;
        
        return size;
    }

    std::string HDF5File::GetFileName() const
    {
        return m_file_name;
    }

    void HDF5File::ReadDataSet(std::string rPath, std::vector<int>& rData, unsigned BlockSize)
    {
        ReadDataSetImpl(rPath, rData, BlockSize);
    }
    
    void HDF5File::ReadDataSet(std::string rPath, std::vector<double>& rData, unsigned BlockSize)
    {
        ReadDataSetImpl(rPath, rData, BlockSize);        
    }
    
    void HDF5File::ReadDataSet(std::string rPath, std::vector<array_1d<double,3>>& rData, unsigned BlockSize)
    {
        ReadDataSetImpl(rPath, rData, BlockSize);        
    }
    

} // // namespace Kratos.