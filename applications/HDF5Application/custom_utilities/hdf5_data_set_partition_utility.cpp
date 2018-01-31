#include "custom_utilities/hdf5_data_set_partition_utility.h"

namespace Kratos
{
namespace HDF5
{

void DataSetPartitionUtility::WritePartition(File& rFile, std::string const& rPath, WriteInfo const& rInfo)
{
    KRATOS_TRY;

    const std::string partition_path = rPath + DataSetPartitionUtility::PartitionSuffix;
    Vector<int> data;
    if (rFile.GetPID() == rFile.GetTotalProcesses() - 1)
    {
        data.resize(2, false);
        data[0] = rInfo.StartIndex;
        data[1] = rInfo.TotalSize;
    }
    else
    {
        data.resize(1, false);
        data[0] = rInfo.StartIndex;
    }

    WriteInfo dummy;
    rFile.WriteDataSet(partition_path, data, dummy);

    KRATOS_CATCH("");
}

void DataSetPartitionUtility::WritePartitionIndependent(File& rFile, std::string const& rPath, Vector<int> const& rPartition)
{
    KRATOS_TRY;

    const std::string partition_path = rPath + DataSetPartitionUtility::PartitionSuffix;
    WriteInfo info;
    rFile.WriteDataSetIndependent(partition_path, rPartition, info);

    KRATOS_CATCH("");
}

bool DataSetPartitionUtility::HasPartition(File& rFile, std::string const& rPath)
{
    KRATOS_TRY;

    const std::string partition_path = rPath + DataSetPartitionUtility::PartitionSuffix;
    return rFile.IsDataSet(partition_path);

    KRATOS_CATCH("");
}

std::tuple<unsigned, unsigned> DataSetPartitionUtility::StartIndexAndBlockSize(File& rFile, std::string const& rPath)
{
    KRATOS_TRY;

    const std::string partition_path = rPath + DataSetPartitionUtility::PartitionSuffix;
    const unsigned my_pid = rFile.GetPID();
    Vector<int> my_partition;
    rFile.ReadDataSet(partition_path, my_partition, my_pid, 2);
    const unsigned start_index = my_partition[0];
    const unsigned block_size = my_partition[1] - my_partition[0];
    return std::make_tuple(start_index, block_size);

    KRATOS_CATCH("");
}

const std::string DataSetPartitionUtility::PartitionSuffix = "_partition";

} // namespace HDF5.
} // namespace Kratos.
