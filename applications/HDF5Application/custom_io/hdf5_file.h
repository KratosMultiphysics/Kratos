//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main author:    Michael Andre, https://github.com/msandre
//

#if !defined(KRATOS_HDF5_FILE_H_INCLUDED)
#define KRATOS_HDF5_FILE_H_INCLUDED

// System includes
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <type_traits>
#include <regex>

// External includes
extern "C" {
#include "hdf5.h"
}

// Project includes
#include "includes/define.h"
#include "containers/array_1d.h"
#include "includes/kratos_parameters.h"

namespace Kratos
{
///@addtogroup HDF5Application
///@{

/// Provides helper functions that don't depend on the representation of a
/// particular class.
struct HDF5Utils
{
    /// Check if string is a valid path.
    /**
     * Valid paths are similar to linux file system with alphanumeric names
     * and possible underscores separated by '/'. All paths are absolute.
     */
    static bool IsPath(std::string Path)
    { 
        return regex_match(Path, std::regex("(/\\w+)+"));
    }

    // Return vector of non-empty substrings separated by a delimiter.
    static std::vector<std::string> Split(std::string Path, char Delimiter)
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
};

class HDF5File
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(HDF5File);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    HDF5File() {}

    explicit HDF5File(Parameters& rParams);

    // Copy constructor.
    HDF5File(const HDF5File& rOther) = delete;

    /// Destructor.
    virtual ~HDF5File();

    // Assignment operator.
    HDF5File& operator=(const HDF5File& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    /// Check if path exists in HDF5 file.
    virtual bool HasPath(std::string Path) const;

    virtual bool IsGroup(std::string Path) const;

    virtual bool IsDataSet(std::string Path) const;

    virtual void CreateGroup(std::string Path);

    virtual void AddPath(std::string Path);

    /// Write a data set to the HDF5 file.
    /**
     * Performs independent write in MPI. Must be called collectively with only
     * one process having non-empty data.
     */
    virtual void WriteDataSet(std::string Path, const std::vector<int>& rData);

    virtual void WriteDataSet(std::string Path, const std::vector<double>& rData);

    virtual void WriteDataSet(std::string Path,
                              const std::vector<array_1d<double, 3>>& rData);
    
    /// Collectively write data set to the HDF5 file.
    /** 
     *  Performs collective write in MPI. The data is written blockwise according to
     *  processor rank. In serial, it calls the default write.
     */
    virtual void WriteDataSetCollective(std::string Path, const std::vector<int>& rData);
                              
    virtual void WriteDataSetCollective(std::string Path, const std::vector<double>& rData);

    virtual void WriteDataSetCollective(std::string Path,
                                        const std::vector<array_1d<double, 3>>& rData);

    virtual std::vector<unsigned> GetDataDimensions(std::string Path) const;

    virtual bool HasIntDataType(std::string Path) const;

    virtual bool HasFloatDataType(std::string Path) const;

    virtual void Flush();

    virtual unsigned GetFileSize() const;

    virtual std::string GetFileName() const;

    /// Read a data set from the HDF5 file.
    /**
     * Performs independent read in MPI.
     */
    virtual void ReadDataSet(std::string Path,
                             std::vector<int>& rData,
                             unsigned StartIndex,
                             unsigned BlockSize);

    virtual void ReadDataSet(std::string Path,
                             std::vector<double>& rData,
                             unsigned StartIndex,
                             unsigned BlockSize);

    virtual void ReadDataSet(std::string Path,
                             std::vector<array_1d<double, 3>>& rData,
                             unsigned StartIndex,
                             unsigned BlockSize);

    // Collectively read data set from the HDF5 file.
    /**
     *  Performs collective read in MPI. In serial, it calls the default read.
     */
    virtual void ReadDataSetCollective(std::string Path,
                                       std::vector<int>& rData,
                                       unsigned StartIndex,
                                       unsigned BlockSize);

    virtual void ReadDataSetCollective(std::string Path,
                                       std::vector<double>& rData,
                                       unsigned StartIndex,
                                       unsigned BlockSize);

    virtual void ReadDataSetCollective(std::string Path,
                                       std::vector<array_1d<double, 3>>& rData,
                                       unsigned StartIndex,
                                       unsigned BlockSize);

    ///@}

protected:
    ///@name Member Variables
    ///@{
    std::string m_file_name;
    hid_t m_file_id = -1; // Default invalid file id.
    int m_echo_level;
    ///@}

private:
    ///@name Private Operations
    ///@{
    template <class T>
    void WriteDataSetImpl(std::string Path, const std::vector<T>& rData)
    {
        // Check that full path does not exist before trying to write data.
        KRATOS_ERROR_IF(HasPath(Path)) << "Path already exists: " << Path << std::endl;

        // Create any missing subpaths.
        auto pos = Path.find_last_of('/');
        std::string sub_path = Path.substr(0, pos);
        AddPath(sub_path);

        // Initialize data space dimensions.
        constexpr bool is_int_type = std::is_same<int, T>::value;
        constexpr bool is_double_type = std::is_same<double, T>::value;
        constexpr bool is_array_1d_type = std::is_same<array_1d<double, 3>, T>::value;
        constexpr unsigned ndims = (!is_array_1d_type) ? 1 : 2;
        hsize_t dims[ndims] = {0};
        dims[0] = rData.size(); // Set first data space dimension.
        if (is_array_1d_type)
            dims[1] = 3; // Set second data space dimension (array_1d<double,3>).
        
        // Set the data type.
        hid_t dtype_id;
        if (is_int_type)
            dtype_id = H5T_NATIVE_INT;
        else if (is_double_type)
            dtype_id = H5T_NATIVE_DOUBLE;
        else if (is_array_1d_type)
            dtype_id = H5T_NATIVE_DOUBLE;
        else
            static_assert(is_int_type || is_double_type || is_array_1d_type,
                          "Unsupported data type.");

        // Create and write the data set.
        hid_t dspace_id = H5Screate_simple(ndims, dims, nullptr);
        hid_t dset_id = H5Dcreate(m_file_id, Path.c_str(), dtype_id, dspace_id,
                                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        KRATOS_ERROR_IF(dset_id < 0) << "H5Dcreate failed." << std::endl;
        KRATOS_ERROR_IF(H5Dwrite(dset_id, dtype_id, H5S_ALL, H5S_ALL,
                                 H5P_DEFAULT, rData.data()) < 0)
            << "H5Dwrite failed." << std::endl;
        KRATOS_ERROR_IF(H5Dclose(dset_id) < 0) << "H5Dclose failed." << std::endl;
        KRATOS_ERROR_IF(H5Sclose(dspace_id) < 0) << "H5Sclose failed." << std::endl;
    }

    template <class T>
    void ReadDataSetImpl(std::string Path, std::vector<T>& rData, unsigned StartIndex, unsigned BlockSize)
    {
        // Check that full path exists.
        KRATOS_ERROR_IF_NOT(IsDataSet(Path))
            << "Path is not a data set: " << Path << std::endl;

        constexpr bool is_int_type = std::is_same<int, T>::value;
        constexpr bool is_double_type = std::is_same<double, T>::value;
        constexpr bool is_array_1d_type = std::is_same<array_1d<double, 3>, T>::value;
        constexpr unsigned ndims = (!is_array_1d_type) ? 1 : 2;

        // Check consistency of file's data set dimensions.
        std::vector<unsigned> file_space_dims = GetDataDimensions(Path);
        KRATOS_ERROR_IF(file_space_dims.size() != ndims)
            << "Invalid data set dimension." << std::endl;
        KRATOS_ERROR_IF(StartIndex + BlockSize > file_space_dims[0])
            << "StartIndex (" << StartIndex << ") + BlockSize (" << BlockSize
            << ") > size of data set (" << file_space_dims[0] << ")." << std::endl;
        if (is_array_1d_type)
            KRATOS_ERROR_IF(file_space_dims[1] != 3)
                << "Invalid data set dimension." << std::endl;

        if (rData.size() != BlockSize)
            rData.resize(BlockSize);

        // Initialize memory space dimensions.
        hsize_t start[ndims] = {0}, mem_dims[ndims];
        start[0] = StartIndex;
        mem_dims[0] = BlockSize; // Set first dimension.
        if (is_array_1d_type)
            mem_dims[1] = 3; // Set second dimension.

        // Set the data type.
        hid_t dtype_id;
        if (is_int_type)
        {
            KRATOS_ERROR_IF_NOT(HasIntDataType(Path))
                << "Data type is not int: " << Path << std::endl;
            dtype_id = H5T_NATIVE_INT;
        }
        else if (is_double_type)
        {
            KRATOS_ERROR_IF_NOT(HasFloatDataType(Path))
                << "Data type is not float: " << Path << std::endl;
            dtype_id = H5T_NATIVE_DOUBLE;
        }
        else if (is_array_1d_type)
        {
            KRATOS_ERROR_IF_NOT(HasFloatDataType(Path))
                << "Data type is not float: " << Path << std::endl;
            dtype_id = H5T_NATIVE_DOUBLE;
        }
        else
            static_assert(is_int_type || is_double_type || is_array_1d_type,
                          "Unsupported data type.");

        hid_t dset_id = H5Dopen(m_file_id, Path.c_str(), H5P_DEFAULT);
        KRATOS_ERROR_IF(dset_id < 0) << "H5Dopen failed." << std::endl;
        hid_t file_space_id = H5Dget_space(dset_id);
        hid_t mem_space_id =
            H5Screate_simple(ndims, mem_dims, nullptr);
        KRATOS_ERROR_IF(H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, start, nullptr, mem_dims, nullptr) < 0)
            << "H5Sselect_hyperslab failed." << std::endl;
        KRATOS_ERROR_IF(H5Dread(dset_id, dtype_id, mem_space_id, file_space_id, H5P_DEFAULT, rData.data()) < 0)
            << "H5Dread failed." << std::endl;
        KRATOS_ERROR_IF(H5Dclose(dset_id) < 0) << "H5Dclose failed." << std::endl;
        KRATOS_ERROR_IF(H5Sclose(file_space_id) < 0) << "H5Sclose failed." << std::endl;
        KRATOS_ERROR_IF(H5Sclose(mem_space_id) < 0) << "H5Sclose failed." << std::endl;
    }
    ///@}
};

///@} addtogroup
} // namespace Kratos.

#endif // KRATOS_HDF5_FILE_H_INCLUDED defined