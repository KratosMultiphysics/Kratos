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
#include <type_traits>

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
    HDF5File(Parameters& rParams);

    // Copy constructor.
    HDF5File(const HDF5File& rOther) = delete;

    /// Destructor.
    virtual ~HDF5File();

    // Assignment operator.
    HDF5File& operator=(const HDF5File& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    virtual bool IsPath(std::string rPath) const;

    virtual bool HasPath(std::string rPath) const;

    virtual bool IsGroup(std::string rPath) const;

    virtual bool IsDataSet(std::string rPath) const;

    virtual void CreateGroup(std::string rPath);

    /// Write a data set to the HDF5 file.
    /** This function must always be called colletively in MPI. For more
     *  than one process, data blocks are ordered by rank.
     */
    virtual void WriteDataSet(std::string rPath, const std::vector<int>& rData);

    virtual void WriteDataSet(std::string rPath, const std::vector<double>& rData);

    virtual void WriteDataSet(std::string rPath,
                              const std::vector<array_1d<double, 3>>& rData);

    virtual std::vector<unsigned> GetDataDimensions(std::string rPath);

    virtual bool HasIntDataType(std::string rPath);

    virtual bool HasFloatDataType(std::string rPath);

    virtual void Flush();

    virtual unsigned GetFileSize() const;

    virtual std::string GetFileName() const;

    /// Read a data set from the HDF5 file.
    /** This function must always be called colletively in MPI. For more
     *  than one process, data blocks are ordered by rank. For example,
     *  for 10 processes each with a BlockSize of 1000, the process
     *  with rank 2 will read 1000 elements beginning at index 2000.
     */
    virtual void ReadDataSet(std::string rPath, std::vector<int>& rData, unsigned BlockSize);

    virtual void ReadDataSet(std::string rPath, std::vector<double>& rData, unsigned BlockSize);

    virtual void ReadDataSet(std::string rPath,
                             std::vector<array_1d<double, 3>>& rData,
                             unsigned BlockSize);

    ///@}

protected:
private:
    ///@name Member Variables
    ///@{
    std::string m_file_name;
    hid_t m_file_id;
    int m_echo_level;
    ///@}

    ///@name Private Operations
    ///@{
    template <class T>
    void WriteDataSetImpl(std::string rPath, const std::vector<T>& rData)
    {
        KRATOS_ERROR_IF_NOT(IsPath(rPath)) << "Invalid path: " << rPath << std::endl;

        // Create any missing subpaths.
        std::string sub_path;
        decltype(rPath.size()) pos = 0; /* Make sure pos has same size as
           std::string::npos. */
        while (pos < rPath.size()) // Check each link in the path.
        {
            pos = rPath.find('/', ++pos);
            if (pos != std::string::npos)
                sub_path = rPath.substr(0, pos); // Check current subpath.
            else
                break; // Exit loop after all links are present.

            if (HasPath(sub_path) == false)
                CreateGroup(sub_path); // Create missing link.
            else
                KRATOS_ERROR_IF_NOT(IsGroup(sub_path))
                    << "Path exists and is not a group: " << sub_path << std::endl;
        }

        // Check that full path does not exist before trying to write data.
        KRATOS_ERROR_IF(HasPath(rPath)) << "Path already exists: " << rPath << std::endl;

        // Create and write data set.
        int ndims = 1; // Default rank is 1 (scalar data type).
        hid_t dtype_id;
        hsize_t dims[2] = {0};
        dims[0] = rData.size(); // Set first dataspace dimension.
        constexpr bool is_int_type = std::is_same<int, T>::value;
        constexpr bool is_double_type = std::is_same<double, T>::value;
        constexpr bool is_array_1d_type = std::is_same<array_1d<double, 3>, T>::value;
        if (is_int_type)
            dtype_id = H5T_NATIVE_INT;
        else if (is_double_type)
            dtype_id = H5T_NATIVE_DOUBLE;
        else if (is_array_1d_type)
        {
            dtype_id = H5T_NATIVE_DOUBLE;
            ndims = 2; // Increase data space rank to 2 for array_1d data type.
            dims[1] = 3; // Second dimension is length of array_1d<double, 3>.
        }
        else
            static_assert(is_int_type || is_double_type || is_array_1d_type,
                          "Unsupported data type.");

        hid_t dspace_id = H5Screate_simple(ndims, dims, nullptr);
        hid_t dset_id = H5Dcreate(m_file_id, rPath.c_str(), dtype_id, dspace_id,
                                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        KRATOS_ERROR_IF(dset_id < 0) << "H5Dcreate failed." << std::endl;
        KRATOS_ERROR_IF(H5Dwrite(dset_id, dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, rData.data()) < 0) << "H5Dwrite failed." << std::endl;
        KRATOS_ERROR_IF(H5Dclose(dset_id) < 0) << "H5Dclose failed." << std::endl;
        KRATOS_ERROR_IF(H5Sclose(dspace_id) < 0) << "H5Sclose failed." << std::endl;
    }

    template <class T>
    void ReadDataSetImpl(std::string rPath, std::vector<T>& rData, unsigned BlockSize)
    {
        // Check that full path exists.
        KRATOS_ERROR_IF_NOT(IsDataSet(rPath))
            << "Path is not a data set: " << rPath << std::endl;

        // Check data types and dimensions.
        hid_t dtype_id;
        std::vector<unsigned> file_space_dims = GetDataDimensions(rPath);
        KRATOS_ERROR_IF(file_space_dims.size() == 0) << "Invalid data set." << std::endl;
        KRATOS_ERROR_IF(BlockSize > file_space_dims[0])
            << "BlockSize exceeds data set dimension by: "
            << BlockSize - file_space_dims[0] << std::endl;
        if (rData.size() != BlockSize)
            rData.resize(BlockSize);

        // Define the hyperslab.
        hsize_t start[] = {0, 0}, stride[] = {1, 1}, count[] = {0, 0};
        count[0] = BlockSize;
        std::vector<hsize_t> mem_space_dims(file_space_dims.size());
        mem_space_dims[0] = BlockSize; // First dimension is BlockSize.
        constexpr bool is_int_type = std::is_same<int, T>::value;
        constexpr bool is_double_type = std::is_same<double, T>::value;
        constexpr bool is_array_1d_type = std::is_same<array_1d<double, 3>, T>::value;
        if (is_int_type)
        {
            KRATOS_ERROR_IF_NOT(HasIntDataType(rPath))
                << "Data type is not int: " << rPath << std::endl;
            KRATOS_ERROR_IF(file_space_dims.size() != 1)
                << "Invalid data set dimension." << std::endl;
            dtype_id = H5T_NATIVE_INT;
        }
        else if (is_double_type)
        {
            KRATOS_ERROR_IF_NOT(HasFloatDataType(rPath))
                << "Data type is not float: " << rPath << std::endl;
            KRATOS_ERROR_IF(file_space_dims.size() != 1)
                << "Invalid data set dimension." << std::endl;
            dtype_id = H5T_NATIVE_DOUBLE;
        }
        else if (is_array_1d_type)
        {
            KRATOS_ERROR_IF_NOT(HasFloatDataType(rPath))
                << "Data type is not float: " << rPath << std::endl;
            KRATOS_ERROR_IF(file_space_dims.size() != 2 || file_space_dims[1] != 3)
                << "Invalid data set dimension." << std::endl;
            dtype_id = H5T_NATIVE_DOUBLE;
            count[1] = 3; // Second dimension is length of array_1d<double, 3>.
            mem_space_dims[1] = 3; 
        }
        else
            static_assert(is_int_type || is_double_type || is_array_1d_type,
                          "Unsupported data type.");

        hid_t dset_id = H5Dopen(m_file_id, rPath.c_str(), H5P_DEFAULT);
        KRATOS_ERROR_IF(dset_id < 0) << "H5Dopen failed." << std::endl;
        hid_t file_space_id = H5Dget_space(dset_id);
        hid_t mem_space_id =
            H5Screate_simple(mem_space_dims.size(), mem_space_dims.data(), nullptr);
        KRATOS_ERROR_IF(H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, start, stride, count, nullptr) < 0)
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