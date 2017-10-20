#include "hdf5_file_mpi.h"

#include <regex>
#include <cassert>

namespace Kratos
{
HDF5FileMPI::HDF5FileMPI(Parameters& rParams)
{
    KRATOS_TRY;
    
        Parameters default_params(R"(
                {
                    "file_name" : "PLEASE_SPECIFY_HDF5_FILENAME",
                    "file_access_mode": "exclusive",
                    "file_driver": "mpio",
                    "echo_level" : 0
                })");
    
        rParams.RecursivelyValidateAndAssignDefaults(default_params);
    
        m_file_name = rParams["file_name"].GetString();

        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
        std::string file_access_mode = rParams["file_access_mode"].GetString();
        hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    
        std::string file_driver = rParams["file_driver"].GetString();
        if (file_driver == "mpio")
        {
            KRATOS_ERROR_IF(H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL) < 0)
                << "H5Pset_fapl_mpio failed." << std::endl;
        }
        else
            KRATOS_ERROR << "Unsupported \"file_driver\": " << file_driver << std::endl;
    
        if (file_access_mode == "exclusive")
            m_file_id = H5Fcreate(m_file_name.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, fapl_id);
        else if (file_access_mode == "truncate")
            m_file_id = H5Fcreate(m_file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
        else
        {
            if (file_driver != "core")
                if (rank == 0)
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

HDF5FileMPI::~HDF5FileMPI()
{
}

void HDF5FileMPI::WriteDataSet(std::string Path, const std::vector<int>& rData)
{
    KRATOS_TRY;
    WriteDataSetImpl(Path, rData);
    KRATOS_CATCH("");
}

void HDF5FileMPI::WriteDataSet(std::string Path, const std::vector<double>& rData)
{
    KRATOS_TRY;
    WriteDataSetImpl(Path, rData);
    KRATOS_CATCH("");
}

void HDF5FileMPI::WriteDataSet(std::string Path, const std::vector<array_1d<double, 3>>& rData)
{
    KRATOS_TRY;
    WriteDataSetImpl(Path, rData);
    KRATOS_CATCH("");
}

void HDF5FileMPI::ReadDataSet(std::string Path, std::vector<int>& rData, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetImpl(Path, rData, BlockSize);
    KRATOS_CATCH("");
}

void HDF5FileMPI::ReadDataSet(std::string Path, std::vector<double>& rData, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetImpl(Path, rData, BlockSize);
    KRATOS_CATCH("");
}

void HDF5FileMPI::ReadDataSet(std::string Path, std::vector<array_1d<double, 3>>& rData, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetImpl(Path, rData, BlockSize);
    KRATOS_CATCH("");
}

} // // namespace Kratos.