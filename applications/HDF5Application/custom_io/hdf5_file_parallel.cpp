#include "hdf5_file_parallel.h"

namespace Kratos
{
namespace HDF5
{
FileParallel::FileParallel(Parameters& rParams) : File(rParams)
{
}

void FileParallel::WriteDataSet(std::string Path, const Vector<int>& rData)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(Path, rData, DataTransferMode::collective);
    KRATOS_CATCH("");
}

void FileParallel::WriteDataSet(std::string Path, const Vector<double>& rData)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(Path, rData, DataTransferMode::collective);
    KRATOS_CATCH("");
}

void FileParallel::WriteDataSet(std::string Path, const Vector<array_1d<double, 3>>& rData)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(Path, rData, DataTransferMode::collective);
    KRATOS_CATCH("");
}

void FileParallel::WriteDataSet(std::string Path, const Matrix<int>& rData)
{
    KRATOS_TRY;
    WriteDataSetMatrixImpl(Path, rData, DataTransferMode::collective);
    KRATOS_CATCH("");
}

void FileParallel::WriteDataSet(std::string Path, const Matrix<double>& rData)
{
    KRATOS_TRY;
    WriteDataSetMatrixImpl(Path, rData, DataTransferMode::collective);
    KRATOS_CATCH("");
}

void FileParallel::WriteDataPartition(std::string Path, const Vector<int>& rData)
{
    KRATOS_TRY;
    WriteDataPartitionVectorImpl(Path, rData);
    KRATOS_CATCH("");
}

void FileParallel::WriteDataPartition(std::string Path, const Vector<double>& rData)
{
    KRATOS_TRY;
    WriteDataPartitionVectorImpl(Path, rData);
    KRATOS_CATCH("");
}

void FileParallel::WriteDataPartition(std::string Path, const Vector<array_1d<double,3>>& rData)
{
    KRATOS_TRY;
    WriteDataPartitionVectorImpl(Path, rData);
    KRATOS_CATCH("");
}

void FileParallel::WriteDataPartition(std::string Path, const Matrix<int>& rData)
{
    KRATOS_TRY;
    WriteDataPartitionMatrixImpl(Path, rData);
    KRATOS_CATCH("");
}

void FileParallel::WriteDataPartition(std::string Path, const Matrix<double>& rData)
{
    KRATOS_TRY;
    WriteDataPartitionMatrixImpl(Path, rData);
    KRATOS_CATCH("");
}

void FileParallel::WriteDataSetIndependent(std::string Path, const Vector<int>& rData)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(Path, rData, DataTransferMode::independent);
    KRATOS_CATCH("");
}

void FileParallel::WriteDataSetIndependent(std::string Path, const Vector<double>& rData)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(Path, rData, DataTransferMode::independent);
    KRATOS_CATCH("");
}

void FileParallel::WriteDataSetIndependent(std::string Path, const Vector<array_1d<double, 3>>& rData)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(Path, rData, DataTransferMode::independent);
    KRATOS_CATCH("");
}

void FileParallel::WriteDataSetIndependent(std::string Path, const Matrix<int>& rData)
{
    KRATOS_TRY;
    WriteDataSetMatrixImpl(Path, rData, DataTransferMode::independent);
    KRATOS_CATCH("");
}

void FileParallel::WriteDataSetIndependent(std::string Path, const Matrix<double>& rData)
{
    KRATOS_TRY;
    WriteDataSetMatrixImpl(Path, rData, DataTransferMode::independent);
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

void FileParallel::ReadDataSet(std::string Path, Vector<int>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(Path, rData, StartIndex, BlockSize, DataTransferMode::collective);
    KRATOS_CATCH("");
}

void FileParallel::ReadDataSet(std::string Path, Vector<double>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(Path, rData, StartIndex, BlockSize, DataTransferMode::collective);
    KRATOS_CATCH("");
}

void FileParallel::ReadDataSet(std::string Path,
                              Vector<array_1d<double, 3>>& rData,
                              unsigned StartIndex,
                              unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(Path, rData, StartIndex, BlockSize, DataTransferMode::collective);
    KRATOS_CATCH("");
}

void FileParallel::ReadDataSet(std::string Path, Matrix<int>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetMatrixImpl(Path, rData, StartIndex, BlockSize, DataTransferMode::collective);
    KRATOS_CATCH("");
}

void FileParallel::ReadDataSet(std::string Path, Matrix<double>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetMatrixImpl(Path, rData, StartIndex, BlockSize, DataTransferMode::collective);
    KRATOS_CATCH("");
}

void FileParallel::ReadDataSetIndependent(std::string Path,
                                        Vector<int>& rData,
                                        unsigned StartIndex,
                                        unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(Path, rData, StartIndex, BlockSize, DataTransferMode::independent);
    KRATOS_CATCH("");
}

void FileParallel::ReadDataSetIndependent(std::string Path,
                                        Vector<double>& rData,
                                        unsigned StartIndex,
                                        unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(Path, rData, StartIndex, BlockSize, DataTransferMode::independent);
    KRATOS_CATCH("");
}

void FileParallel::ReadDataSetIndependent(std::string Path,
                                        Vector<array_1d<double, 3>>& rData,
                                        unsigned StartIndex,
                                        unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(Path, rData, StartIndex, BlockSize, DataTransferMode::independent);
    KRATOS_CATCH("");
}

void FileParallel::ReadDataSetIndependent(std::string Path,
                                              Matrix<int>& rData,
                                              unsigned StartIndex,
                                              unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetMatrixImpl(Path, rData, StartIndex, BlockSize, DataTransferMode::independent);
    KRATOS_CATCH("");
}

void FileParallel::ReadDataSetIndependent(std::string Path,
                                              Matrix<double>& rData,
                                              unsigned StartIndex,
                                              unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetMatrixImpl(Path, rData, StartIndex, BlockSize, DataTransferMode::independent);
    KRATOS_CATCH("");
}
} // namespace HDF5.
} // namespace Kratos.