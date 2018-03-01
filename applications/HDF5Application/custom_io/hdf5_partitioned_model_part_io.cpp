#include "hdf5_partitioned_model_part_io.h"

#include "utilities/openmp_utils.h"
#include "custom_io/hdf5_points_data.h"
#include "custom_io/hdf5_connectivities_data.h"
#include "custom_utilities/factor_elements_and_conditions_utility.h"
#include "custom_utilities/hdf5_data_set_partition_utility.h"

namespace Kratos
{
namespace HDF5
{
PartitionedModelPartIO::PartitionedModelPartIO(File::Pointer pFile, std::string const& rPrefix)
: ModelPartIOBase(pFile, rPrefix)
{
    KRATOS_TRY;

    Check();

    KRATOS_CATCH("");
}

bool PartitionedModelPartIO::ReadNodes(NodesContainerType& rNodes)
{
    KRATOS_TRY;

    rNodes.clear();

    unsigned local_start_index, local_block_size;
    unsigned ghost_start_index, ghost_block_size;

    File& r_file = *mpFile;

    std::tie(local_start_index, local_block_size) = StartIndexAndBlockSize(r_file, mPrefix + "/Nodes/Local");
    std::tie(ghost_start_index, ghost_block_size) = StartIndexAndBlockSize(r_file, mPrefix + "/Nodes/Ghost");
    rNodes.reserve(local_block_size + ghost_block_size);

    // Read local nodes.
    Internals::PointsData local_points;
    local_points.ReadData(r_file, mPrefix + "/Nodes/Local", local_start_index, local_block_size);
    local_points.CreateNodes(rNodes);
    local_points.Clear();
    
    // Read ghost nodes.
    Internals::PointsData ghost_points;
    ghost_points.ReadData(r_file, mPrefix + "/Nodes/Ghost", ghost_start_index, ghost_block_size);
    ghost_points.CreateNodes(rNodes);
    ghost_points.Clear();

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
    Internals::PointsData local_points;
    local_points.SetData(local_nodes);
    WriteInfo info;
    local_points.WriteData(r_file, mPrefix + "/Nodes/Local", info);
    WritePartitionTable(r_file, mPrefix + "/Nodes/Local", info);
    local_points.Clear();

    // Write ghost nodes.
    Internals::PointsData ghost_points;
    ghost_points.SetData(ghost_nodes);
    ghost_points.WriteData(r_file, mPrefix + "/Nodes/Ghost", info);
    WritePartitionTable(r_file, mPrefix + "/Nodes/Ghost", info);
    WritePartitionIndex(mPrefix + "/Nodes/Ghost", ghost_nodes);
    ghost_points.Clear();

    KRATOS_CATCH("");
}

void PartitionedModelPartIO::ReadElements(NodesContainerType& rNodes,
                                          PropertiesContainerType& rProperties,
                                          ElementsContainerType& rElements)
{
    KRATOS_TRY;

    unsigned start_index, block_size;
    rElements.clear();
    std::vector<std::string> group_names;
    mpFile->GetGroupNames(mPrefix + "/Elements", group_names);

    for (const auto& r_name : group_names)
    {
        Internals::ConnectivitiesData connectivities;
        std::tie(start_index, block_size) = StartIndexAndBlockSize(
            *mpFile, mPrefix + "/Elements/" + r_name);
        connectivities.ReadData(*mpFile, mPrefix + "/Elements/" + r_name, start_index, block_size);
        connectivities.CreateEntities(rNodes, rProperties, rElements);
    }

    KRATOS_CATCH("");
}

void PartitionedModelPartIO::WriteElements(ElementsContainerType const& rElements)
{
    KRATOS_TRY;

    std::vector<std::string> names;
    std::vector<ElementsContainerType> factored_elements;
    FactorElements(rElements, names, factored_elements);

    WriteInfo info;
    for (unsigned int i = 0; i < names.size(); ++i)
    {
        Internals::ConnectivitiesData connectivities;
        // For partitioned elements, the local container may be empty. Therefore,
        // we explicitly provide the element name here.
        connectivities.SetData(names[i], factored_elements[i]);
        connectivities.WriteData(*mpFile, mPrefix + "/Elements/" + connectivities.Name(), info);
        const int size = info.TotalSize;
        mpFile->WriteAttribute(mPrefix + "/Elements/" + connectivities.Name(), "Size", size);
        WritePartitionTable(*mpFile, mPrefix + "/Elements/" + connectivities.Name(), info);
    }

    KRATOS_CATCH("");
}

void PartitionedModelPartIO::ReadConditions(NodesContainerType& rNodes,
                                            PropertiesContainerType& rProperties,
                                            ConditionsContainerType& rConditions)
{
    KRATOS_TRY;

    unsigned start_index, block_size;
    rConditions.clear();
    std::vector<std::string> group_names;
    mpFile->GetGroupNames(mPrefix + "/Conditions", group_names);

    for (const auto& r_name : group_names)
    {
        Internals::ConnectivitiesData connectivities;
        std::tie(start_index, block_size) = StartIndexAndBlockSize(
            *mpFile, mPrefix + "/Conditions/" + r_name);
        connectivities.ReadData(*mpFile, mPrefix + "/Conditions/" + r_name, start_index, block_size);
        connectivities.CreateEntities(rNodes, rProperties, rConditions);
    }

    KRATOS_CATCH("");
}

void PartitionedModelPartIO::WriteConditions(ConditionsContainerType const& rConditions)
{
    KRATOS_TRY;

    std::vector<std::string> names;
    std::vector<ConditionsContainerType> factored_conditions;
    FactorConditions(rConditions, names, factored_conditions);

    WriteInfo info;
    for (unsigned i = 0; i < names.size();  ++i)
    {
        Internals::ConnectivitiesData connectivities;
        // For partitioned conditions, the local container may be empty. Therefore,
        // we explicitly provide the condition name here.
        connectivities.SetData(names[i], factored_conditions[i]);
        connectivities.WriteData(*mpFile, mPrefix + "/Conditions/" + connectivities.Name(), info);
        const int size = info.TotalSize;
        mpFile->WriteAttribute(mPrefix + "/Conditions/" + connectivities.Name(), "Size", size);
        WritePartitionTable(*mpFile, mPrefix + "/Conditions/" + connectivities.Name(), info);
    }

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
    std::tie(start_index, block_size) = StartIndexAndBlockSize(*mpFile, rPath);
    Vector<int> partition_ids, node_ids;
    mpFile->ReadDataSet(rPath + "/PARTITION_INDEX", partition_ids, start_index, block_size);
    mpFile->ReadDataSet(rPath + "/Ids", node_ids, start_index, block_size);
    for (unsigned i = 0; i < node_ids.size(); ++i)
        r_nodes[node_ids[i]].FastGetSolutionStepValue(PARTITION_INDEX) = partition_ids[i];

    KRATOS_CATCH("");
}

} // namespace HDF5.
} // namespace Kratos.
