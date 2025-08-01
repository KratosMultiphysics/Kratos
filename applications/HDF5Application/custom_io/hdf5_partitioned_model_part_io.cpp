#include "hdf5_partitioned_model_part_io.h"

#include "utilities/openmp_utils.h"
#include "custom_io/hdf5_points_data.h"
#include "custom_io/hdf5_connectivities_data.h"
#include "custom_utilities/factor_elements_and_conditions_utility.h"
#include "custom_utilities/hdf5_data_set_partition_utility.h"
#include "custom_utilities/container_io_utils.h"

namespace Kratos
{
namespace HDF5
{
PartitionedModelPartIO::PartitionedModelPartIO(File::Pointer pFile, std::string const& rPrefix)
: ModelPartIO(pFile, rPrefix)
{
    KRATOS_TRY;

    Check();

    KRATOS_CATCH("");
}

bool PartitionedModelPartIO::ReadNodes(NodesContainerType& rNodes)
{
    KRATOS_TRY;

    rNodes.clear();

    // read local nodes
    BaseType::ReadNodes(rNodes);

    // Read ghost nodes.
    Internals::PointsData<Internals::NodesIO> ghost_points(mPrefix + "/Nodes/Ghost", mpFile);
    ghost_points.Read(rNodes, Internals::NodesIO{});

    return true;

    KRATOS_CATCH("");
}

void PartitionedModelPartIO::WriteNodes(NodesContainerType const& rNodes)
{
    KRATOS_TRY;

    NodesContainerType local_nodes;
    NodesContainerType ghost_nodes;
    const unsigned num_nodes = rNodes.size();
    local_nodes.reserve(num_nodes);
    ghost_nodes.reserve(0.1 * num_nodes);
    File& r_file = *mpFile;

    // Divide nodes into local and global containers.
    int my_pid = r_file.GetPID();
    for (auto it_node = rNodes.begin(); it_node != rNodes.end(); ++it_node)
    {
        if (it_node->FastGetSolutionStepValue(PARTITION_INDEX) == my_pid)
            local_nodes.push_back(*it_node.base());
        else
            ghost_nodes.push_back(*it_node.base());
    }

    // Write local nodes.
    BaseType::WriteNodes(local_nodes);

    // Write ghost nodes.
    Internals::PointsData<Internals::NodesIO> ghost_points(mPrefix + "/Nodes/Ghost", mpFile);
    ghost_points.Write(ghost_nodes, Internals::NodesIO{}, Parameters(R"({})"));
    WritePartitionIndex(mPrefix + "/Nodes/Ghost", ghost_nodes);

    KRATOS_CATCH("");
}

void PartitionedModelPartIO::ReadModelPart(ModelPart& rModelPart)
{
    KRATOS_TRY;

    BaseType::ReadModelPart(rModelPart);
    ReadAndAssignPartitionIndex(mPrefix + "/Nodes/Ghost", rModelPart);

    KRATOS_CATCH("");
}

void PartitionedModelPartIO::Check()
{
    if (mpFile->GetTotalProcesses() == 1)
        KRATOS_ERROR << "Using PartitionedModelPartIO with single process file access." << std::endl;
}

void PartitionedModelPartIO::WritePartitionIndex(const std::string& rPath, NodesContainerType const& rGhostNodes)
{
    KRATOS_TRY;

    Vector<int> partition_ids(rGhostNodes.size());
    const int num_threads = OpenMPUtils::GetNumThreads();
    OpenMPUtils::PartitionVector partition;
    OpenMPUtils::DivideInPartitions(rGhostNodes.size(), num_threads, partition);
#pragma omp parallel
    {
        const int thread_id = OpenMPUtils::ThisThread();
        NodesContainerType::const_iterator it = rGhostNodes.begin() + partition[thread_id];
        for (auto i = partition[thread_id]; i < partition[thread_id + 1]; ++i)
        {
            partition_ids[i] = it->FastGetSolutionStepValue(PARTITION_INDEX);
            ++it;
        }
    }
    WriteInfo info;
    mpFile->WriteDataSet(rPath + "/PARTITION_INDEX", partition_ids, info);

    KRATOS_CATCH("");
}

void PartitionedModelPartIO::ReadAndAssignPartitionIndex(const std::string& rPath, ModelPart& rModelPart) const
{
    KRATOS_TRY;

    const int num_threads = OpenMPUtils::GetNumThreads();
    OpenMPUtils::PartitionVector partition;
    NodesContainerType& r_nodes = rModelPart.Nodes();

    // First assign my partition id to ALL nodes.
    const int my_pid = rModelPart.GetCommunicator().MyPID();
    OpenMPUtils::DivideInPartitions(r_nodes.size(), num_threads, partition);
#pragma omp parallel
    {
        const int thread_id = OpenMPUtils::ThisThread();
        NodesContainerType::const_iterator it = r_nodes.begin() + partition[thread_id];
        for (auto i = partition[thread_id]; i < partition[thread_id + 1]; ++i)
        {
            it->FastGetSolutionStepValue(PARTITION_INDEX) = my_pid;
            ++it;
        }
    }

    // Read and assign partition ids for ghosts. If this is slow, try storing the
    // partition index in the node's non-historical container when the ghost nodes and
    // local nodes are read. Then copy it to the solution step data after the buffer
    // is initialized.
    unsigned start_index, block_size;
    std::tie(start_index, block_size) = StartIndexAndBlockSize(rPath);
    Vector<int> partition_ids, node_ids;
    mpFile->ReadDataSet(rPath + "/PARTITION_INDEX", partition_ids, start_index, block_size);
    mpFile->ReadDataSet(rPath + "/Ids", node_ids, start_index, block_size);
    for (unsigned i = 0; i < node_ids.size(); ++i)
        r_nodes[node_ids[i]].FastGetSolutionStepValue(PARTITION_INDEX) = partition_ids[i];

    KRATOS_CATCH("");
}

std::tuple<unsigned, unsigned> PartitionedModelPartIO::StartIndexAndBlockSize(std::string const& rPath) const
{
    KRATOS_TRY;
    unsigned start_index, block_size;
    std::tie(start_index, block_size) = HDF5::StartIndexAndBlockSize(*mpFile, rPath);
    return std::make_tuple(start_index, block_size);
    KRATOS_CATCH("");
}

void PartitionedModelPartIO::StoreWriteInfo(std::string const& rPath, WriteInfo const& rInfo)
{
    KRATOS_TRY;
    WritePartitionTable(*mpFile, rPath, rInfo);
    const int size = rInfo.TotalSize;
    mpFile->WriteAttribute(rPath, "Size", size);
    KRATOS_CATCH("");
}

} // namespace HDF5.
} // namespace Kratos.
