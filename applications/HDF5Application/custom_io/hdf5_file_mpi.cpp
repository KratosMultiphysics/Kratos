#include "hdf5_file_mpi.h"

#include <regex>
#include <cassert>

namespace Kratos
{
HDF5FileMPI::HDF5FileMPI(Parameters& rParams) : HDF5File(rParams)
{
}

void HDF5FileMPI::WriteDataSet(std::string Path, const std::vector<int>& rData)
{
    KRATOS_TRY;
    WriteDataSetImpl(Path, rData, DataTransferMode::collective);
    KRATOS_CATCH("");
}

void HDF5FileMPI::WriteDataSet(std::string Path, const std::vector<double>& rData)
{
    KRATOS_TRY;
    WriteDataSetImpl(Path, rData, DataTransferMode::collective);
    KRATOS_CATCH("");
}

void HDF5FileMPI::WriteDataSet(std::string Path, const std::vector<array_1d<double, 3>>& rData)
{
    KRATOS_TRY;
    WriteDataSetImpl(Path, rData, DataTransferMode::collective);
    KRATOS_CATCH("");
}

void HDF5FileMPI::WriteDataSetIndependent(std::string Path, const std::vector<int>& rData)
{
    KRATOS_TRY;
    WriteDataSetImpl(Path, rData, DataTransferMode::independent);
    KRATOS_CATCH("");
}

void HDF5FileMPI::WriteDataSetIndependent(std::string Path, const std::vector<double>& rData)
{
    KRATOS_TRY;
    WriteDataSetImpl(Path, rData, DataTransferMode::independent);
    KRATOS_CATCH("");
}

void HDF5FileMPI::WriteDataSetIndependent(std::string Path, const std::vector<array_1d<double, 3>>& rData)
{
    KRATOS_TRY;
    WriteDataSetImpl(Path, rData, DataTransferMode::independent);
    KRATOS_CATCH("");
}

void HDF5FileMPI::ReadDataSet(std::string Path, std::vector<int>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetImpl(Path, rData, StartIndex, BlockSize, DataTransferMode::collective);
    KRATOS_CATCH("");
}

void HDF5FileMPI::ReadDataSet(std::string Path, std::vector<double>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetImpl(Path, rData, StartIndex, BlockSize, DataTransferMode::collective);
    KRATOS_CATCH("");
}

void HDF5FileMPI::ReadDataSet(std::string Path,
                              std::vector<array_1d<double, 3>>& rData,
                              unsigned StartIndex,
                              unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetImpl(Path, rData, StartIndex, BlockSize, DataTransferMode::collective);
    KRATOS_CATCH("");
}

void HDF5FileMPI::ReadDataSetIndependent(std::string Path,
                                        std::vector<int>& rData,
                                        unsigned StartIndex,
                                        unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetImpl(Path, rData, StartIndex, BlockSize, DataTransferMode::independent);
    KRATOS_CATCH("");
}

void HDF5FileMPI::ReadDataSetIndependent(std::string Path,
                                        std::vector<double>& rData,
                                        unsigned StartIndex,
                                        unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetImpl(Path, rData, StartIndex, BlockSize, DataTransferMode::independent);
    KRATOS_CATCH("");
}

void HDF5FileMPI::ReadDataSetIndependent(std::string Path,
                                        std::vector<array_1d<double, 3>>& rData,
                                        unsigned StartIndex,
                                        unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetImpl(Path, rData, StartIndex, BlockSize, DataTransferMode::independent);
    KRATOS_CATCH("");
}

} // // namespace Kratos.