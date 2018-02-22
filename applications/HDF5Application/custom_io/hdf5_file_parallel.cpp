#include "hdf5_file_parallel.h"

#include "includes/kratos_parameters.h"

namespace Kratos
{
namespace HDF5
{
FileParallel::FileParallel(Parameters& rSettings) : File(rSettings)
{
}

void FileParallel::WriteDataSet(const std::string& rPath, const Vector<int>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(rPath, rData, DataTransferMode::collective, rInfo);
    KRATOS_CATCH("");
}

void FileParallel::WriteDataSet(const std::string& rPath, const Vector<double>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(rPath, rData, DataTransferMode::collective, rInfo);
    KRATOS_CATCH("");
}

void FileParallel::WriteDataSet(const std::string& rPath, const Vector<array_1d<double, 3>>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(rPath, rData, DataTransferMode::collective, rInfo);
    KRATOS_CATCH("");
}

void FileParallel::WriteDataSet(const std::string& rPath, const Matrix<int>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetMatrixImpl(rPath, rData, DataTransferMode::collective, rInfo);
    KRATOS_CATCH("");
}

void FileParallel::WriteDataSet(const std::string& rPath, const Matrix<double>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetMatrixImpl(rPath, rData, DataTransferMode::collective, rInfo);
    KRATOS_CATCH("");
}

void FileParallel::WriteDataSetIndependent(const std::string& rPath, const Vector<int>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(rPath, rData, DataTransferMode::independent, rInfo);
    KRATOS_CATCH("");
}

void FileParallel::WriteDataSetIndependent(const std::string& rPath, const Vector<double>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(rPath, rData, DataTransferMode::independent, rInfo);
    KRATOS_CATCH("");
}

void FileParallel::WriteDataSetIndependent(const std::string& rPath, const Vector<array_1d<double, 3>>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(rPath, rData, DataTransferMode::independent, rInfo);
    KRATOS_CATCH("");
}

void FileParallel::WriteDataSetIndependent(const std::string& rPath, const Matrix<int>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetMatrixImpl(rPath, rData, DataTransferMode::independent, rInfo);
    KRATOS_CATCH("");
}

void FileParallel::WriteDataSetIndependent(const std::string& rPath, const Matrix<double>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetMatrixImpl(rPath, rData, DataTransferMode::independent, rInfo);
    KRATOS_CATCH("");
}

unsigned FileParallel::GetPID() const
{
    int rank, ierr;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    KRATOS_ERROR_IF(ierr != MPI_SUCCESS) << "MPI_Comm_rank failed." << std::endl;
    return static_cast<unsigned>(rank);
}

unsigned FileParallel::GetTotalProcesses() const
{
    int num_proc, ierr;
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    KRATOS_ERROR_IF(ierr != MPI_SUCCESS) << "MPI_Comm_size failed." << std::endl;
    return static_cast<unsigned>(num_proc);
}

void FileParallel::ReadDataSet(const std::string& rPath, Vector<int>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(rPath, rData, StartIndex, BlockSize, DataTransferMode::collective);
    KRATOS_CATCH("");
}

void FileParallel::ReadDataSet(const std::string& rPath, Vector<double>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(rPath, rData, StartIndex, BlockSize, DataTransferMode::collective);
    KRATOS_CATCH("");
}

void FileParallel::ReadDataSet(const std::string& rPath,
                              Vector<array_1d<double, 3>>& rData,
                              unsigned StartIndex,
                              unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(rPath, rData, StartIndex, BlockSize, DataTransferMode::collective);
    KRATOS_CATCH("");
}

void FileParallel::ReadDataSet(const std::string& rPath, Matrix<int>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetMatrixImpl(rPath, rData, StartIndex, BlockSize, DataTransferMode::collective);
    KRATOS_CATCH("");
}

void FileParallel::ReadDataSet(const std::string& rPath, Matrix<double>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetMatrixImpl(rPath, rData, StartIndex, BlockSize, DataTransferMode::collective);
    KRATOS_CATCH("");
}

void FileParallel::ReadDataSetIndependent(const std::string& rPath,
                                        Vector<int>& rData,
                                        unsigned StartIndex,
                                        unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(rPath, rData, StartIndex, BlockSize, DataTransferMode::independent);
    KRATOS_CATCH("");
}

void FileParallel::ReadDataSetIndependent(const std::string& rPath,
                                        Vector<double>& rData,
                                        unsigned StartIndex,
                                        unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(rPath, rData, StartIndex, BlockSize, DataTransferMode::independent);
    KRATOS_CATCH("");
}

void FileParallel::ReadDataSetIndependent(const std::string& rPath,
                                        Vector<array_1d<double, 3>>& rData,
                                        unsigned StartIndex,
                                        unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(rPath, rData, StartIndex, BlockSize, DataTransferMode::independent);
    KRATOS_CATCH("");
}

void FileParallel::ReadDataSetIndependent(const std::string& rPath,
                                              Matrix<int>& rData,
                                              unsigned StartIndex,
                                              unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetMatrixImpl(rPath, rData, StartIndex, BlockSize, DataTransferMode::independent);
    KRATOS_CATCH("");
}

void FileParallel::ReadDataSetIndependent(const std::string& rPath,
                                              Matrix<double>& rData,
                                              unsigned StartIndex,
                                              unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetMatrixImpl(rPath, rData, StartIndex, BlockSize, DataTransferMode::independent);
    KRATOS_CATCH("");
}
} // namespace HDF5.
} // namespace Kratos.