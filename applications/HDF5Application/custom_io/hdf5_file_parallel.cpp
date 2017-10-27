#include "hdf5_file_parallel.h"

namespace Kratos
{
HDF5FileParallel::HDF5FileParallel(Parameters& rParams) : HDF5File(rParams)
{
}

void HDF5FileParallel::WriteDataSet(std::string Path, const Vector<int>& rData)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(Path, rData, DataTransferMode::collective);
    KRATOS_CATCH("");
}

void HDF5FileParallel::WriteDataSet(std::string Path, const Vector<double>& rData)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(Path, rData, DataTransferMode::collective);
    KRATOS_CATCH("");
}

void HDF5FileParallel::WriteDataSet(std::string Path, const Vector<array_1d<double, 3>>& rData)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(Path, rData, DataTransferMode::collective);
    KRATOS_CATCH("");
}

void HDF5FileParallel::WriteDataPartition(std::string Path, const Vector<int>& rData)
{
    KRATOS_TRY;
    WriteDataPartitionVectorImpl(Path, rData);
    KRATOS_CATCH("");
}

void HDF5FileParallel::WriteDataPartition(std::string Path, const Vector<double>& rData)
{
    KRATOS_TRY;
    WriteDataPartitionVectorImpl(Path, rData);
    KRATOS_CATCH("");
}

void HDF5FileParallel::WriteDataPartition(std::string Path, const Vector<array_1d<double,3>>& rData)
{
    KRATOS_TRY;
    WriteDataPartitionVectorImpl(Path, rData);
    KRATOS_CATCH("");
}

void HDF5FileParallel::WriteDataSetIndependent(std::string Path, const Vector<int>& rData)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(Path, rData, DataTransferMode::independent);
    KRATOS_CATCH("");
}

void HDF5FileParallel::WriteDataSetIndependent(std::string Path, const Vector<double>& rData)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(Path, rData, DataTransferMode::independent);
    KRATOS_CATCH("");
}

void HDF5FileParallel::WriteDataSetIndependent(std::string Path, const Vector<array_1d<double, 3>>& rData)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(Path, rData, DataTransferMode::independent);
    KRATOS_CATCH("");
}

unsigned HDF5FileParallel::GetPID() const
{
    int rank, ierr;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    KRATOS_ERROR_IF(ierr != MPI_SUCCESS) << "MPI_Comm_rank failed." << std::endl;
    return static_cast<unsigned>(rank);
}

unsigned HDF5FileParallel::GetTotalProcesses() const
{
    int num_proc, ierr;
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    KRATOS_ERROR_IF(ierr != MPI_SUCCESS) << "MPI_Comm_size failed." << std::endl;
    return static_cast<unsigned>(num_proc);
}

void HDF5FileParallel::ReadDataSet(std::string Path, Vector<int>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(Path, rData, StartIndex, BlockSize, DataTransferMode::collective);
    KRATOS_CATCH("");
}

void HDF5FileParallel::ReadDataSet(std::string Path, Vector<double>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(Path, rData, StartIndex, BlockSize, DataTransferMode::collective);
    KRATOS_CATCH("");
}

void HDF5FileParallel::ReadDataSet(std::string Path,
                              Vector<array_1d<double, 3>>& rData,
                              unsigned StartIndex,
                              unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(Path, rData, StartIndex, BlockSize, DataTransferMode::collective);
    KRATOS_CATCH("");
}

void HDF5FileParallel::ReadDataSetIndependent(std::string Path,
                                        Vector<int>& rData,
                                        unsigned StartIndex,
                                        unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(Path, rData, StartIndex, BlockSize, DataTransferMode::independent);
    KRATOS_CATCH("");
}

void HDF5FileParallel::ReadDataSetIndependent(std::string Path,
                                        Vector<double>& rData,
                                        unsigned StartIndex,
                                        unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(Path, rData, StartIndex, BlockSize, DataTransferMode::independent);
    KRATOS_CATCH("");
}

void HDF5FileParallel::ReadDataSetIndependent(std::string Path,
                                        Vector<array_1d<double, 3>>& rData,
                                        unsigned StartIndex,
                                        unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(Path, rData, StartIndex, BlockSize, DataTransferMode::independent);
    KRATOS_CATCH("");
}

} // // namespace Kratos.