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

#if !defined(KRATOS_HDF5_FILE_MPI_H_INCLUDED)
#define KRATOS_HDF5_FILE_MPI_H_INCLUDED

// System includes

// External includes

// Project includes

// Application includes
#include "custom_io/hdf5_file.h"

namespace Kratos
{
///@addtogroup HDF5Application
///@{

class HDF5FileMPI : public HDF5File
{
    enum class DataTransferMode { independent, collective };
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(HDF5FileMPI);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    HDF5FileMPI(Parameters& rParams);

    // Copy constructor.
    HDF5FileMPI(const HDF5FileMPI& rOther) = delete;

    /// Destructor.
    ~HDF5FileMPI() override;

    /// Assignment operator.
    HDF5FileMPI& operator=(const HDF5FileMPI& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    void WriteDataSet(std::string Path, const std::vector<int>& rData) override;

    void WriteDataSet(std::string Path, const std::vector<double>& rData) override;

    void WriteDataSet(std::string Path, const std::vector<array_1d<double, 3>>& rData) override;

    void WriteDataSetCollective(std::string Path, const std::vector<int>& rData) override;

    void WriteDataSetCollective(std::string Path, const std::vector<double>& rData) override;

    void WriteDataSetCollective(std::string Path,
                                const std::vector<array_1d<double, 3>>& rData) override;

    void ReadDataSet(std::string Path, std::vector<int>& rData, unsigned StartIndex, unsigned BlockSize) override;

    void ReadDataSet(std::string Path,
                     std::vector<double>& rData,
                     unsigned StartIndex,
                     unsigned BlockSize) override;

    void ReadDataSet(std::string Path,
                     std::vector<array_1d<double, 3>>& rData,
                     unsigned StartIndex,
                     unsigned BlockSize) override;

    void ReadDataSetCollective(std::string Path,
                               std::vector<int>& rData,
                               unsigned StartIndex,
                               unsigned BlockSize) override;

    void ReadDataSetCollective(std::string Path,
                               std::vector<double>& rData,
                               unsigned StartIndex,
                               unsigned BlockSize) override;

    void ReadDataSetCollective(std::string Path,
                               std::vector<array_1d<double, 3>>& rData,
                               unsigned StartIndex,
                               unsigned BlockSize) override;
    ///@}

protected:
    ///@name Protected Operations
    ///@{

    ///@}

private:
    ///@name Private Operations
    ///@{
    template <class T>
    void WriteDataSetImpl(std::string Path, const std::vector<T>& rData, DataTransferMode Mode)
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
        hsize_t local_dims[ndims], global_dims[ndims];
        // Set first data space dimension.
        local_dims[0] = rData.size();
        if (Mode == DataTransferMode::collective)
        {
            unsigned send_buf, recv_buf;
            send_buf = local_dims[0];
            MPI_Allreduce(&send_buf, &recv_buf, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
            global_dims[0] = recv_buf;
        }
        else
            global_dims[0] = local_dims[0];
        if (is_array_1d_type)
            local_dims[1] = global_dims[1] = 3; // Set second data space dimension.
        
        hsize_t local_start[ndims];
        if (Mode == DataTransferMode::collective)
        { 
            unsigned send_buf, recv_buf;
            send_buf = local_dims[0];
            // Get global position where local data set ends.
            MPI_Scan(&send_buf, &recv_buf, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
            // Set global position where local data set starts.
            local_start[0] = recv_buf - local_dims[0];
        }
        else
            local_start[0] = 0;
 
        if (is_array_1d_type)
            local_start[1] = 0;
        
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
        hid_t dxpl_id = H5Pcreate(H5P_DATASET_XFER);
        if (Mode == DataTransferMode::collective)
            H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE);
        else
            H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_INDEPENDENT);
        hid_t mspace_id = H5Screate_simple(ndims, local_dims, nullptr);
        hid_t fspace_id = H5Screate_simple(ndims, global_dims, nullptr);
        H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, local_start, nullptr,
                            local_dims, nullptr);
        hid_t dset_id = H5Dcreate(m_file_id, Path.c_str(), dtype_id,
                                  fspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        KRATOS_ERROR_IF(dset_id < 0) << "H5Dcreate failed." << std::endl;
        KRATOS_ERROR_IF(H5Dwrite(dset_id, dtype_id, mspace_id, fspace_id,
                                 dxpl_id, rData.data()) < 0)
            << "H5Dwrite failed." << std::endl;
        KRATOS_ERROR_IF(H5Pclose(dxpl_id) < 0) << "H5Pclose failed." << std::endl;
        KRATOS_ERROR_IF(H5Sclose(mspace_id) < 0) << "H5Sclose failed." << std::endl;
        KRATOS_ERROR_IF(H5Sclose(fspace_id) < 0) << "H5Sclose failed." << std::endl;
        KRATOS_ERROR_IF(H5Dclose(dset_id) < 0) << "H5Dclose failed." << std::endl;
    }

    template <class T>
    void ReadDataSetImpl(std::string Path,
                         std::vector<T>& rData,
                         unsigned StartIndex,
                         unsigned BlockSize,
                         DataTransferMode Mode)
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
        KRATOS_ERROR_IF(file_space_dims.size() != ndims) << "Invalid data set dimension." << std::endl;
        KRATOS_ERROR_IF(StartIndex + BlockSize > file_space_dims[0])
        << "StartIndex (" << StartIndex << ") + BlockSize (" << BlockSize
        << ") > size of data set (" << file_space_dims[0] << ")." << std::endl;
        hsize_t local_mem_dims[ndims];
        // Set first memory space dimension.
        local_mem_dims[0] = BlockSize;
        if (is_array_1d_type)
            local_mem_dims[1] = 3; // Set second dimension.
        if (is_array_1d_type)
            KRATOS_ERROR_IF(file_space_dims[1] != 3)
                << "Invalid data set dimension." << std::endl;

        if (rData.size() != BlockSize)
            rData.resize(BlockSize);

        // Set global position where local data set starts.
        hsize_t local_start[ndims];
        local_start[0] = StartIndex;
        if (is_array_1d_type)
            local_start[1] = 0;

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

        hid_t dxpl_id = H5Pcreate(H5P_DATASET_XFER);
        if (Mode == DataTransferMode::collective)
            H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE);
        else
            H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_INDEPENDENT);
        hid_t dset_id = H5Dopen(m_file_id, Path.c_str(), H5P_DEFAULT);
        KRATOS_ERROR_IF(dset_id < 0) << "H5Dopen failed." << std::endl;
        hid_t file_space_id = H5Dget_space(dset_id);
        hid_t mem_space_id = H5Screate_simple(ndims, local_mem_dims, nullptr);
        KRATOS_ERROR_IF(H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, local_start,
                                            nullptr, local_mem_dims, nullptr) < 0)
            << "H5Sselect_hyperslab failed." << std::endl;
        KRATOS_ERROR_IF(H5Dread(dset_id, dtype_id, mem_space_id, file_space_id,
                                dxpl_id, rData.data()) < 0)
            << "H5Dread failed." << std::endl;
        KRATOS_ERROR_IF(H5Pclose(dxpl_id) < 0) << "H5Pclose failed." << std::endl;
        KRATOS_ERROR_IF(H5Dclose(dset_id) < 0) << "H5Dclose failed." << std::endl;
        KRATOS_ERROR_IF(H5Sclose(file_space_id) < 0) << "H5Sclose failed." << std::endl;
        KRATOS_ERROR_IF(H5Sclose(mem_space_id) < 0) << "H5Sclose failed." << std::endl;
    }
    ///@}
};

///@} addtogroup
} // namespace Kratos.

#endif // KRATOS_HDF5_FILE_MPI_H_INCLUDED defined