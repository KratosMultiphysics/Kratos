#include "hdf5_file_serial.h"

namespace Kratos
{
namespace HDF5
{
FileSerial::FileSerial(Parameters& rSettings) : File(rSettings)
{
}

void FileSerial::WriteDataSet(const std::string& rPath, const Vector<int>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(rPath, rData, rInfo);
    KRATOS_CATCH("");
}

void FileSerial::WriteDataSet(const std::string& rPath, const Vector<double>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(rPath, rData, rInfo);
    KRATOS_CATCH("");
}

void FileSerial::WriteDataSet(const std::string& rPath, const Vector<array_1d<double, 3>>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(rPath, rData, rInfo);
    KRATOS_CATCH("");
}

void FileSerial::WriteDataSet(const std::string& rPath, const Matrix<int>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetMatrixImpl(rPath, rData, rInfo);
    KRATOS_CATCH("");
}

void FileSerial::WriteDataSet(const std::string& rPath, const Matrix<double>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetMatrixImpl(rPath, rData, rInfo);
    KRATOS_CATCH("");
}

void FileSerial::WriteDataSetIndependent(const std::string& rPath, const Vector<int>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(rPath, rData, rInfo);
    KRATOS_CATCH("");
}

void FileSerial::WriteDataSetIndependent(const std::string& rPath, const Vector<double>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(rPath, rData, rInfo);
    KRATOS_CATCH("");
}

void FileSerial::WriteDataSetIndependent(const std::string& rPath, const Vector<array_1d<double, 3>>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(rPath, rData, rInfo);
    KRATOS_CATCH("");
}

void FileSerial::WriteDataSetIndependent(const std::string& rPath, const Matrix<int>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetMatrixImpl(rPath, rData, rInfo);
    KRATOS_CATCH("");
}

void FileSerial::WriteDataSetIndependent(const std::string& rPath, const Matrix<double>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetMatrixImpl(rPath, rData, rInfo);
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

void FileSerial::ReadDataSet(const std::string& rPath, Vector<int>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(rPath, rData, StartIndex, BlockSize);
    KRATOS_CATCH("");
}

void FileSerial::ReadDataSet(const std::string& rPath, Vector<double>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(rPath, rData, StartIndex, BlockSize);
    KRATOS_CATCH("");
}

void FileSerial::ReadDataSet(const std::string& rPath, Vector<array_1d<double, 3>>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(rPath, rData, StartIndex, BlockSize);
    KRATOS_CATCH("");
}

void FileSerial::ReadDataSet(const std::string& rPath, Matrix<int>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetMatrixImpl(rPath, rData, StartIndex, BlockSize);
    KRATOS_CATCH("");
}

void FileSerial::ReadDataSet(const std::string& rPath, Matrix<double>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetMatrixImpl(rPath, rData, StartIndex, BlockSize);
    KRATOS_CATCH("");
}

void FileSerial::ReadDataSetIndependent(const std::string& rPath, Vector<int>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(rPath, rData, StartIndex, BlockSize);
    KRATOS_CATCH("");
}

void FileSerial::ReadDataSetIndependent(const std::string& rPath, Vector<double>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(rPath, rData, StartIndex, BlockSize);
    KRATOS_CATCH("");
}

void FileSerial::ReadDataSetIndependent(const std::string& rPath, Vector<array_1d<double, 3>>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(rPath, rData, StartIndex, BlockSize);
    KRATOS_CATCH("");
}

void FileSerial::ReadDataSetIndependent(const std::string& rPath, Matrix<int>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetMatrixImpl(rPath, rData, StartIndex, BlockSize);
    KRATOS_CATCH("");
}

void FileSerial::ReadDataSetIndependent(const std::string& rPath, Matrix<double>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetMatrixImpl(rPath, rData, StartIndex, BlockSize);
    KRATOS_CATCH("");
}

} // namespace HDF5.
} // namespace Kratos.