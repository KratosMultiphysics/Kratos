#include "custom_utilities/hdf5_data_set_partition_utility.h"

namespace Kratos
{
namespace HDF5
{

namespace
{
const std::string PartitionSuffix = "_partition";
}

void WritePartitionTable(File& rFile, std::string const& rPath, WriteInfo const& rInfo)
{
    KRATOS_TRY;

    const std::string partition_path = rPath + PartitionSuffix;
    Vector<int> data;
    if (rFile.GetPID() == rFile.GetTotalProcesses() - 1) {
        data.resize(2, false);
        data[0] = rInfo.StartIndex;
        data[1] = rInfo.TotalSize;
    } else  {
        data.resize(1, false);
        data[0] = rInfo.StartIndex;
    }

    WriteInfo dummy;
    rFile.WriteDataSet(partition_path, data, dummy);
    rFile.WriteAttribute(partition_path, "Size", static_cast<int>(rInfo.TotalSize));

    KRATOS_CATCH("");
}

bool HasPartitionTable(File& rFile, std::string const& rPath)
{
    KRATOS_TRY;

    const std::string partition_path = rPath + PartitionSuffix;
    return rFile.IsDataSet(partition_path);

    KRATOS_CATCH("");
}

std::tuple<unsigned, unsigned> StartIndexAndBlockSize(File& rFile, std::string const& rPath)
{
    KRATOS_TRY;

    const std::string partition_path = rPath + PartitionSuffix;

    if (rFile.GetDataCommunicator().IsDistributed()) {
        const unsigned my_pid = rFile.GetPID();
        Vector<int> my_partition;
        rFile.ReadDataSet(partition_path, my_partition, my_pid, 2);
        const unsigned start_index = my_partition[0];
        const unsigned block_size = my_partition[1] - my_partition[0];
        return std::make_tuple(start_index, block_size);
    } else {
        int size;
        rFile.ReadAttribute(partition_path, "Size", size);
        return std::make_tuple(0, size);
    }

    KRATOS_CATCH("");
}


} // namespace HDF5.
} // namespace Kratos.
