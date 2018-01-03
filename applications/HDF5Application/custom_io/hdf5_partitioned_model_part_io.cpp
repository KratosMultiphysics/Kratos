#include "hdf5_partitioned_model_part_io.h"

#include "utilities/openmp_utils.h"
#include "custom_utilities/hdf5_points_data.h"
#include "custom_io/hdf5_nodal_solution_step_variables_io.h"
#include "custom_io/hdf5_data_value_container_io.h"

namespace Kratos
{
namespace HDF5
{
PartitionedModelPartIO::PartitionedModelPartIO(Parameters Settings, File::Pointer pFile)
: ModelPartIOBase(Settings, pFile)
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

    std::tie(local_start_index, local_block_size) =
        GetPartitionStartIndexAndBlockSize(mPrefix + "/Nodes/Local/Partition");
    std::tie(ghost_start_index, ghost_block_size) =
        GetPartitionStartIndexAndBlockSize(mPrefix + "/Nodes/Ghost/Partition");
    rNodes.reserve(local_block_size + ghost_block_size);

    File& r_file = GetFile();

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

    File& r_file = GetFile();

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
    local_points.WriteData(r_file, mPrefix + "/Nodes/Local");
    r_file.WriteDataPartition(mPrefix + "/Nodes/Local/Partition", local_points.GetIds());
    local_points.Clear();

    // Write ghost nodes.
    Internals::PointsData ghost_points;
    ghost_points.SetData(ghost_nodes);
    ghost_points.WriteData(r_file, mPrefix + "/Nodes/Ghost");
    r_file.WriteDataPartition(mPrefix + "/Nodes/Ghost/Partition", ghost_points.GetIds());
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
    File& r_file = GetFile();

    Internals::ConnectivitiesInput<ElementType> elem_inputs(mElementIO);
    for (auto& r_item : elem_inputs)
    {
        std::tie(start_index, block_size) =
                GetPartitionStartIndexAndBlockSize(r_item.Path + "/Partition");
        r_item.ReadConnectivities(r_file, rNodes, rProperties, start_index, block_size, rElements);
    }

    KRATOS_CATCH("");
}

void PartitionedModelPartIO::WriteElements(ElementsContainerType const& rElements)
{
    KRATOS_TRY;

    Internals::ConnectivitiesOutput<ElementType> elem_outputs(mElementIO, rElements);
    File& r_file = GetFile();
    for (auto& r_item : elem_outputs)
    {
        r_item.WriteConnectivities(r_file);
        r_file.WriteDataPartition(r_item.Path + "/Partition", r_item.Connectivities.GetIds());
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

    File& r_file = GetFile();

    Internals::ConnectivitiesInput<ConditionType> cond_inputs(mConditionIO);
    for (auto& r_item : cond_inputs)
    {
        std::tie(start_index, block_size) =
                GetPartitionStartIndexAndBlockSize(r_item.Path + "/Partition");
        r_item.ReadConnectivities(r_file, rNodes, rProperties, start_index, block_size, rConditions);
    }

    KRATOS_CATCH("");
}

void PartitionedModelPartIO::WriteConditions(ConditionsContainerType const& rConditions)
{
    KRATOS_TRY;

    Internals::ConnectivitiesOutput<ConditionType> cond_outputs(mConditionIO, rConditions);
    File& r_file = GetFile();
    for (auto& r_item : cond_outputs)
    {
        r_item.WriteConnectivities(r_file);
        r_file.WriteDataPartition(r_item.Path + "/Partition", r_item.Connectivities.GetIds());
    }

    KRATOS_CATCH("");
}

void PartitionedModelPartIO::ReadModelPart(ModelPart& rModelPart)
{
    KRATOS_TRY;

    ReadProperties(rModelPart.rProperties());
    Internals::DataValueContainerIO process_info_io(mPrefix + "/ProcessInfo", mpFile);
    process_info_io.ReadDataValueContainer(rModelPart.GetProcessInfo());
    ReadNodes(rModelPart.Nodes());
    ReadElements(rModelPart.Nodes(), rModelPart.rProperties(), rModelPart.Elements());
    ReadConditions(rModelPart.Nodes(), rModelPart.rProperties(), rModelPart.Conditions());
    Internals::NodalSolutionStepVariablesIO nodal_variables_io(mPrefix, mpFile);
    nodal_variables_io.ReadAndAssignVariablesList(rModelPart);
    nodal_variables_io.ReadAndAssignBufferSize(rModelPart);
    ReadAndAssignPartitionIndex(mPrefix + "/Nodes/Ghost", rModelPart);

    KRATOS_CATCH("");
}

void PartitionedModelPartIO::Check()
{
    if (GetFile().GetTotalProcesses() == 1)
        KRATOS_ERROR << "Using PartitionedModelPartIO with single process file access." << std::endl;
}

std::tuple<unsigned, unsigned> PartitionedModelPartIO::GetPartitionStartIndexAndBlockSize(std::string Path) const
{
    unsigned my_pid = GetFile().GetPID();
    Vector<int> my_partition;
    GetFile().ReadDataSet(Path, my_partition, my_pid, 2);
    unsigned start_index = my_partition[0];
    unsigned block_size = my_partition[1] - my_partition[0];
    return std::make_tuple(start_index, block_size);
}

void PartitionedModelPartIO::WritePartitionIndex(std::string Path, NodesContainerType const& rGhostNodes)
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
    mpFile->WriteDataSet(Path + "/PARTITION_INDEX", partition_ids);

    KRATOS_CATCH("");
}

void PartitionedModelPartIO::ReadAndAssignPartitionIndex(std::string Path, ModelPart& rModelPart) const
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
    std::tie(start_index, block_size) = GetPartitionStartIndexAndBlockSize(Path + "/Partition");
    Vector<int> partition_ids, node_ids;
    GetFile().ReadDataSet(Path + "/PARTITION_INDEX", partition_ids, start_index, block_size);
    GetFile().ReadDataSet(Path + "/Ids", node_ids, start_index, block_size);
    for (unsigned i = 0; i < node_ids.size(); ++i)
        r_nodes[node_ids[i]].FastGetSolutionStepValue(PARTITION_INDEX) = partition_ids[i];

    KRATOS_CATCH("");
}

} // namespace HDF5.
} // namespace Kratos.
