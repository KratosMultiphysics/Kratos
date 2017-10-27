#include "hdf5_file_serial.h"

namespace Kratos
{
HDF5FileSerial::HDF5FileSerial(Parameters& rParams) : HDF5File(rParams)
{
}

void HDF5FileSerial::WriteDataSet(std::string Path, const Vector<int>& rData)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(Path, rData);
    KRATOS_CATCH("");
}

void HDF5FileSerial::WriteDataSet(std::string Path, const Vector<double>& rData)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(Path, rData);
    KRATOS_CATCH("");
}

void HDF5FileSerial::WriteDataSet(std::string Path, const Vector<array_1d<double, 3>>& rData)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(Path, rData);
    KRATOS_CATCH("");
}

void HDF5FileSerial::WriteDataPartition(std::string Path, const Vector<int>& rData)
{
    KRATOS_TRY;
    WriteDataPartitionVectorImpl(Path, rData);
    KRATOS_CATCH("");
}

void HDF5FileSerial::WriteDataPartition(std::string Path, const Vector<double>& rData)
{
    KRATOS_TRY;
    WriteDataPartitionVectorImpl(Path, rData);
    KRATOS_CATCH("");
}

void HDF5FileSerial::WriteDataPartition(std::string Path, const Vector<array_1d<double,3>>& rData)
{
    KRATOS_TRY;
    WriteDataPartitionVectorImpl(Path, rData);
    KRATOS_CATCH("");
}

void HDF5FileSerial::WriteDataSetIndependent(std::string Path, const Vector<int>& rData)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(Path, rData);
    KRATOS_CATCH("");
}

void HDF5FileSerial::WriteDataSetIndependent(std::string Path, const Vector<double>& rData)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(Path, rData);
    KRATOS_CATCH("");
}

void HDF5FileSerial::WriteDataSetIndependent(std::string Path, const Vector<array_1d<double, 3>>& rData)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(Path, rData);
    KRATOS_CATCH("");
}

unsigned HDF5FileSerial::GetPID() const
{
    return 0;
}

unsigned HDF5FileSerial::GetTotalProcesses() const
{
    return 1;
}

void HDF5FileSerial::ReadDataSet(std::string Path, Vector<int>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(Path, rData, StartIndex, BlockSize);
    KRATOS_CATCH("");
}

void HDF5FileSerial::ReadDataSet(std::string Path, Vector<double>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(Path, rData, StartIndex, BlockSize);
    KRATOS_CATCH("");
}

void HDF5FileSerial::ReadDataSet(std::string Path, Vector<array_1d<double, 3>>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(Path, rData, StartIndex, BlockSize);
    KRATOS_CATCH("");
}

void HDF5FileSerial::ReadDataSetIndependent(std::string Path, Vector<int>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(Path, rData, StartIndex, BlockSize);
    KRATOS_CATCH("");
}

void HDF5FileSerial::ReadDataSetIndependent(std::string Path, Vector<double>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(Path, rData, StartIndex, BlockSize);
    KRATOS_CATCH("");
}

void HDF5FileSerial::ReadDataSetIndependent(std::string Path, Vector<array_1d<double, 3>>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(Path, rData, StartIndex, BlockSize);
    KRATOS_CATCH("");
}

} // // namespace Kratos.