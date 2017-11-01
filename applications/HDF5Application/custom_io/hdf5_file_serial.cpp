#include "hdf5_file_serial.h"

namespace Kratos
{
namespace HDF5
{
FileSerial::FileSerial(Parameters& rParams) : File(rParams)
{
}

void FileSerial::WriteDataSet(std::string Path, const Vector<int>& rData)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(Path, rData);
    KRATOS_CATCH("");
}

void FileSerial::WriteDataSet(std::string Path, const Vector<double>& rData)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(Path, rData);
    KRATOS_CATCH("");
}

void FileSerial::WriteDataSet(std::string Path, const Vector<array_1d<double, 3>>& rData)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(Path, rData);
    KRATOS_CATCH("");
}

void FileSerial::WriteDataSet(std::string Path, const Matrix<int>& rData)
{
    KRATOS_TRY;
    WriteDataSetMatrixImpl(Path, rData);
    KRATOS_CATCH("");
}

void FileSerial::WriteDataSet(std::string Path, const Matrix<double>& rData)
{
    KRATOS_TRY;
    WriteDataSetMatrixImpl(Path, rData);
    KRATOS_CATCH("");
}


void FileSerial::WriteDataPartition(std::string Path, const Vector<int>& rData)
{
    KRATOS_TRY;
    WriteDataPartitionVectorImpl(Path, rData);
    KRATOS_CATCH("");
}

void FileSerial::WriteDataPartition(std::string Path, const Vector<double>& rData)
{
    KRATOS_TRY;
    WriteDataPartitionVectorImpl(Path, rData);
    KRATOS_CATCH("");
}

void FileSerial::WriteDataPartition(std::string Path, const Vector<array_1d<double,3>>& rData)
{
    KRATOS_TRY;
    WriteDataPartitionVectorImpl(Path, rData);
    KRATOS_CATCH("");
}

void FileSerial::WriteDataPartition(std::string Path, const Matrix<int>& rData)
{
    KRATOS_TRY;
    WriteDataPartitionMatrixImpl(Path, rData);
    KRATOS_CATCH("");
}

void FileSerial::WriteDataPartition(std::string Path, const Matrix<double>& rData)
{
    KRATOS_TRY;
    WriteDataPartitionMatrixImpl(Path, rData);
    KRATOS_CATCH("");
}

void FileSerial::WriteDataSetIndependent(std::string Path, const Vector<int>& rData)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(Path, rData);
    KRATOS_CATCH("");
}

void FileSerial::WriteDataSetIndependent(std::string Path, const Vector<double>& rData)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(Path, rData);
    KRATOS_CATCH("");
}

void FileSerial::WriteDataSetIndependent(std::string Path, const Vector<array_1d<double, 3>>& rData)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(Path, rData);
    KRATOS_CATCH("");
}

void FileSerial::WriteDataSetIndependent(std::string Path, const Matrix<int>& rData)
{
    KRATOS_TRY;
    WriteDataSetMatrixImpl(Path, rData);
    KRATOS_CATCH("");
}

void FileSerial::WriteDataSetIndependent(std::string Path, const Matrix<double>& rData)
{
    KRATOS_TRY;
    WriteDataSetMatrixImpl(Path, rData);
    KRATOS_CATCH("");
}

unsigned FileSerial::GetPID() const
{
    return 0;
}

unsigned FileSerial::GetTotalProcesses() const
{
    return 1;
}

void FileSerial::ReadDataSet(std::string Path, Vector<int>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(Path, rData, StartIndex, BlockSize);
    KRATOS_CATCH("");
}

void FileSerial::ReadDataSet(std::string Path, Vector<double>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(Path, rData, StartIndex, BlockSize);
    KRATOS_CATCH("");
}

void FileSerial::ReadDataSet(std::string Path, Vector<array_1d<double, 3>>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(Path, rData, StartIndex, BlockSize);
    KRATOS_CATCH("");
}

void FileSerial::ReadDataSet(std::string Path, Matrix<int>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetMatrixImpl(Path, rData, StartIndex, BlockSize);
    KRATOS_CATCH("");
}

void FileSerial::ReadDataSet(std::string Path, Matrix<double>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetMatrixImpl(Path, rData, StartIndex, BlockSize);
    KRATOS_CATCH("");
}

void FileSerial::ReadDataSetIndependent(std::string Path, Vector<int>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(Path, rData, StartIndex, BlockSize);
    KRATOS_CATCH("");
}

void FileSerial::ReadDataSetIndependent(std::string Path, Vector<double>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(Path, rData, StartIndex, BlockSize);
    KRATOS_CATCH("");
}

void FileSerial::ReadDataSetIndependent(std::string Path, Vector<array_1d<double, 3>>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(Path, rData, StartIndex, BlockSize);
    KRATOS_CATCH("");
}

void FileSerial::ReadDataSetIndependent(std::string Path, Matrix<int>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetMatrixImpl(Path, rData, StartIndex, BlockSize);
    KRATOS_CATCH("");
}

void FileSerial::ReadDataSetIndependent(std::string Path, Matrix<double>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetMatrixImpl(Path, rData, StartIndex, BlockSize);
    KRATOS_CATCH("");
}

} // namespace HDF5.
} // namespace Kratos.