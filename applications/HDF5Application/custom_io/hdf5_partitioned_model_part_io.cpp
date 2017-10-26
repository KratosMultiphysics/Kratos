#include "hdf5_partitioned_model_part_io.h"
#include <vector>

namespace Kratos
{

HDF5PartitionedModelPartIO::HDF5PartitionedModelPartIO(Parameters& rParams, HDF5File::Pointer pFile)
: mpFile(pFile)
{
    KRATOS_TRY;
    Check();
    KRATOS_CATCH("");
}

bool HDF5PartitionedModelPartIO::ReadNodes(NodesContainerType& rNodes)
{
    KRATOS_TRY;

    rNodes.clear();
    unsigned local_start_index, local_block_size, ghost_start_index, ghost_block_size;
    std::tie(local_start_index, local_block_size) =
        GetPartitionStartIndexAndBlockSize("/Nodes/Local/Partition");
    std::tie(ghost_start_index, ghost_block_size) =
        GetPartitionStartIndexAndBlockSize("/Nodes/Ghost/Partition");
    rNodes.reserve(local_block_size + ghost_block_size);

    // Read local nodes.
    std::vector<int> local_node_ids;
    std::vector<array_1d<double, 3>> local_node_coords;
    GetFile().ReadDataSet("/Nodes/Local/Id", local_node_ids, local_start_index, local_block_size);
    GetFile().ReadDataSet("/Nodes/Local/Coordinate", local_node_coords, local_start_index, local_block_size);
    for (unsigned i = 0; i < local_block_size; ++i)
    {
        const array_1d<double, 3>& r_coord = local_node_coords[i];
        NodeType::Pointer p_node = boost::make_shared<NodeType>(
            local_node_ids[i], r_coord[0], r_coord[1], r_coord[2]);
        rNodes.push_back(p_node);
    }
    local_node_ids.clear();
    local_node_coords.clear();

    // Read ghost nodes.
    std::vector<int> ghost_node_ids;
    std::vector<array_1d<double, 3>> ghost_node_coords;
    GetFile().ReadDataSet("/Nodes/Ghost/Id", ghost_node_ids, ghost_start_index, ghost_block_size);
    GetFile().ReadDataSet("/Nodes/Ghost/Coordinate", ghost_node_coords, ghost_start_index, ghost_block_size);
    for (unsigned i = 0; i < ghost_block_size; ++i)
    {
        const array_1d<double, 3>& r_coord = ghost_node_coords[i];
        NodeType::Pointer p_node = boost::make_shared<NodeType>(
            ghost_node_ids[i], r_coord[0], r_coord[1], r_coord[2]);
        rNodes.push_back(p_node);
    }
    ghost_node_ids.clear();
    ghost_node_coords.clear();

    return true;
    KRATOS_CATCH("");
}

std::size_t HDF5PartitionedModelPartIO::ReadNodesNumber()
{
    const std::vector<unsigned> dims = GetFile().GetDataDimensions("/Nodes/Local/Id");
    return dims[0];
}

void HDF5PartitionedModelPartIO::WriteNodes(NodesContainerType const& rNodes)
{
    KRATOS_TRY;

    NodesContainerType local_nodes;
    NodesContainerType ghost_nodes;

    const unsigned num_nodes = rNodes.size();
    local_nodes.reserve(num_nodes);
    ghost_nodes.reserve(0.1 * num_nodes);

    // Divide nodes into local and global containers.
    int my_pid = GetFile().GetPID();
    for (auto it_node = rNodes.begin(); it_node != rNodes.end(); ++it_node)
    {
        if (it_node->FastGetSolutionStepValue(PARTITION_INDEX) == my_pid)
            local_nodes.push_back(*it_node.base());
        else
            ghost_nodes.push_back(*it_node.base());
    }

    // Write local nodes.
    std::vector<int> local_node_ids(local_nodes.size());
    std::vector<array_1d<double, 3>> local_node_coords(local_nodes.size());
    unsigned local_pos = 0;
    for (const auto& r_node : local_nodes)
    {
        local_node_ids[local_pos] = r_node.Id();
        local_node_coords[local_pos] = r_node.Coordinates();
        ++local_pos;
    }
    GetFile().WriteDataSet("/Nodes/Local/Id", local_node_ids);
    GetFile().WriteDataSet("/Nodes/Local/Coordinate", local_node_coords);
    GetFile().WriteDataPartition("/Nodes/Local/Partition", local_node_ids);
    local_node_ids.clear();
    local_node_coords.clear();

    // Write ghost nodes.
    std::vector<int> ghost_node_ids(ghost_nodes.size());
    std::vector<int> ghost_node_pids(ghost_nodes.size());
    std::vector<array_1d<double, 3>> ghost_node_coords(ghost_nodes.size());
    unsigned ghost_pos = 0;
    for (const auto& r_node : ghost_nodes)
    {
        ghost_node_ids[ghost_pos] = r_node.Id();
        ghost_node_pids[ghost_pos] = r_node.FastGetSolutionStepValue(PARTITION_INDEX);
        ghost_node_coords[ghost_pos] = r_node.Coordinates();
        ++ghost_pos;
    }
    GetFile().WriteDataSet("/Nodes/Ghost/Id", ghost_node_ids);
    GetFile().WriteDataSet("/Nodes/Ghost/PARTITION_INDEX", ghost_node_pids);
    GetFile().WriteDataSet("/Nodes/Ghost/Coordinate", ghost_node_coords);
    GetFile().WriteDataPartition("/Nodes/Ghost/Partition", ghost_node_ids);
    ghost_node_ids.clear();
    ghost_node_pids.clear();
    ghost_node_coords.clear();
    
    KRATOS_CATCH("");
}

void HDF5PartitionedModelPartIO::ReadElements(NodesContainerType& rNodes,
                                   PropertiesContainerType& rProperties,
                                   ElementsContainerType& rElements)
{
}

std::size_t HDF5PartitionedModelPartIO::ReadElementsConnectivities(ConnectivitiesContainerType& rElementsConnectivities)
{
    return 0;
}

void HDF5PartitionedModelPartIO::WriteElements(ElementsContainerType const& rElements)
{
    
}

void HDF5PartitionedModelPartIO::ReadConditions(NodesContainerType& rNodes,
                                     PropertiesContainerType& rProperties,
                                     ConditionsContainerType& rConditions)
{
}

std::size_t HDF5PartitionedModelPartIO::ReadConditionsConnectivities(ConnectivitiesContainerType& rConditionsConnectivities)
{
    return 0;
}

void HDF5PartitionedModelPartIO::ReadInitialValues(ModelPart& rModelPart)
{
}

void HDF5PartitionedModelPartIO::ReadInitialValues(NodesContainerType& rNodes,
                                        ElementsContainerType& rElements,
                                        ConditionsContainerType& rConditions)
{
}

void HDF5PartitionedModelPartIO::ReadModelPart(ModelPart& rModelPart)
{
}

void HDF5PartitionedModelPartIO::WriteModelPart(ModelPart& rModelPart)
{
}

HDF5File& HDF5PartitionedModelPartIO::GetFile() const
{
    return *mpFile;
}

void HDF5PartitionedModelPartIO::Check()
{
    if (mpFile->GetTotalProcesses() == 1)
        std::cout << "Warning: using HDF5PartitionedModelPartIO with single process file access." << std::endl;
}

std::tuple<unsigned, unsigned> HDF5PartitionedModelPartIO::GetPartitionStartIndexAndBlockSize(std::string Path) const
{
    unsigned my_pid = GetFile().GetPID();
    std::vector<int> my_partition;
    GetFile().ReadDataSet(Path, my_partition, my_pid, 2);
    unsigned start_index = my_partition[0];
    unsigned block_size = my_partition[1] - my_partition[0];
    return std::make_tuple(start_index, block_size);
}

} // namespace Kratos.
