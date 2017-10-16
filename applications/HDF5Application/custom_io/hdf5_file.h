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

    HDF5File(const HDF5File& rOther) = delete;

    /// Destructor.
    virtual ~HDF5File();

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

    virtual std::vector<unsigned int> GetDataDimensions(std::string rPath);

    virtual bool HasIntDataType(std::string rPath);

    virtual bool HasFloatDataType(std::string rPath);

    virtual void Flush();

    virtual unsigned GetFileSize() const;

    virtual std::string GetFileName() const;

    /// Read a data set from the HDF5 file.
    /** This function must always be called colletively in MPI. For more
     *  than one process, data blocks are ordered by rank. For example,
     *  for 10 processes each with a block_size of 1000, the process
     *  with rank 2 will read 1000 elements beginning at index 2000.
     */
    virtual void ReadDataSet(std::string rPath, std::vector<int>& rData, unsigned int block_size);

    virtual void ReadDataSet(std::string rPath, std::vector<double>& rData, unsigned int block_size);

    virtual void ReadDataSet(std::string rPath,
                             std::vector<array_1d<double, 3>>& rData,
                             unsigned int block_size);

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
        static_assert(std::is_same<int, T>::value || std::is_same<double, T>::value ||
                      std::is_same<array_1d<double, 3>, T>::value, "Unsupported data type.");

        // Create any missing path links.
        std::string sub_path;
        unsigned int pos = 1;
        while (pos < rPath.size())
        {
            pos = rPath.find('/', pos);
            if (pos != std::string::npos)
                sub_path = rPath.substr(0, pos); // Check current subpath.
            else
                break; // Exit loop after all links are present.

            if (HasPath(sub_path) == false)
                CreateGroup(sub_path);
            else
                KRATOS_ERROR_IF_NOT(IsGroup(sub_path))
                    << "Path exists and is not a group: " << sub_path << std::endl;
        }

        // Check that full path does not exist.
        KRATOS_ERROR_IF(HasPath(rPath)) << "Path already exists: " << rPath << std::endl;

        // Create and write data set.
        herr_t status;
        int ndims = 1;
        hsize_t dims[2];
        dims[0] = rData.size();
        hid_t dtype_id;
        if (std::is_same<int, T>::value)
            dtype_id = H5T_NATIVE_INT;
        else if (std::is_same<double, T>::value)
            dtype_id = H5T_NATIVE_DOUBLE;
        else if (std::is_same<array_1d<double, 3>, T>::value)
        {
            dtype_id = H5T_NATIVE_DOUBLE;
            ndims = 2;
            dims[1] = 3;
        }
        else
            KRATOS_ERROR << "Unsupported data type." << std::endl;

        hid_t dspace_id = H5Screate_simple(ndims, dims, nullptr);
        KRATOS_ERROR_IF(dspace_id < 0) << "H5Screate_simple failed." << std::endl;
        hid_t dset_id = H5Dcreate(m_file_id, rPath.c_str(), dtype_id, dspace_id,
                                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        KRATOS_ERROR_IF(dset_id < 0) << "H5Dcreate failed." << std::endl;
        status =
            H5Dwrite(dset_id, dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, rData.data());
        KRATOS_ERROR_IF(status < 0) << "H5Dwrite failed." << std::endl;
        status = H5Dclose(dset_id);
        KRATOS_ERROR_IF(status < 0) << "H5Dclose failed." << std::endl;
        status = H5Sclose(dspace_id);
        KRATOS_ERROR_IF(status < 0) << "H5Sclose failed." << std::endl;
    }

    template <class T>
    void ReadDataSetImpl(std::string rPath, std::vector<T>& rData, unsigned int block_size)
    {
        static_assert(std::is_same<int, T>::value || std::is_same<double, T>::value ||
                      std::is_same<array_1d<double, 3>, T>::value, "Unsupported data type.");

        // Check that full path exists.
        KRATOS_ERROR_IF_NOT(IsDataSet(rPath))
            << "Path is not a data set: " << rPath << std::endl;

        // Check data types and dimensions.
        hid_t dtype_id;
        std::vector<unsigned int> file_space_dims = GetDataDimensions(rPath);
        KRATOS_ERROR_IF(file_space_dims.size() == 0) << "Empty data set." << std::endl;
        KRATOS_ERROR_IF(block_size > file_space_dims[1])
            << "block_size exceeds data set dimension by: "
            << block_size - file_space_dims[1] << std::endl;
        if (rData.size() != block_size)
            rData.resize(block_size);

        hsize_t start[2] = {0}, stride[2] = {1}, count[2];
        count[0] = block_size;
        std::vector<hsize_t> mem_space_dims(file_space_dims.size());
        mem_space_dims[0] = block_size;
        if (std::is_same<int, T>::value)
        {
            KRATOS_ERROR_IF_NOT(HasIntDataType(rPath))
                << "Data type is not int: " << rPath << std::endl;
            KRATOS_ERROR_IF(file_space_dims.size() != 1)
                << "Invalid data set dimension." << std::endl;
            dtype_id = H5T_NATIVE_INT;
        }
        else if (std::is_same<double, T>::value)
        {
            KRATOS_ERROR_IF_NOT(HasFloatDataType(rPath))
                << "Data type is not float: " << rPath << std::endl;
            KRATOS_ERROR_IF(file_space_dims.size() != 1)
                << "Invalid data set dimension." << std::endl;
            dtype_id = H5T_NATIVE_DOUBLE;
        }
        else if (std::is_same<array_1d<double, 3>, T>::value)
        {
            KRATOS_ERROR_IF_NOT(HasFloatDataType(rPath))
                << "Data type is not float: " << rPath << std::endl;
            KRATOS_ERROR_IF(file_space_dims.size() != 2 || file_space_dims[1] != 3)
                << "Invalid data set dimension." << std::endl;
            dtype_id = H5T_NATIVE_DOUBLE;
            count[1] = 3;
            mem_space_dims[1] = 3;
        }
        else
            KRATOS_ERROR << "Unsupported data type." << std::endl;

        herr_t status;
        hid_t dset_id = H5Dopen(m_file_id, rPath.c_str(), H5P_DEFAULT);
        KRATOS_ERROR_IF(dset_id < 0) << "H5Dopen failed." << std::endl;
        hid_t file_space_id = H5Dget_space(dset_id);
        KRATOS_ERROR_IF(file_space_id < 0) << "H5Dget_space failed." << std::endl;
        hid_t mem_space_id = H5Screate_simple(
            mem_space_dims.size(), mem_space_dims.data(), nullptr);
        KRATOS_ERROR_IF(mem_space_id < 0) << "H5Screate_simple failed." << std::endl;
        status = H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, start,
                                     stride, count, nullptr);
        KRATOS_ERROR_IF(status < 0) << "H5Sselect_hyperslab failed." << std::endl;
        status = H5Dread(dset_id, dtype_id, mem_space_id, file_space_id,
                         H5P_DEFAULT, rData.data());
        KRATOS_ERROR_IF(status < 0) << "H5Dread failed." << std::endl;
        status = H5Dclose(dset_id);
        KRATOS_ERROR_IF(status < 0) << "H5Dclose failed." << std::endl;
        status = H5Sclose(file_space_id);
        KRATOS_ERROR_IF(status < 0) << "H5Dclose failed." << std::endl;
        status = H5Sclose(mem_space_id);
        KRATOS_ERROR_IF(status < 0) << "H5Dclose failed." << std::endl;
    }
    ///@}
};

///@} addtogroup
} // namespace Kratos.

#endif // KRATOS_HDF5_FILE_H_INCLUDED defined