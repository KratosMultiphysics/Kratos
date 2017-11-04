#include "hdf5_partitioned_model_part_io.h"

#include "includes/kratos_components.h"
#include "custom_utilities/hdf5_points_data.h"
#include "custom_utilities/hdf5_connectivities_data.h"
#include "custom_utilities/hdf5_pointer_bins_utility.h"
#include "custom_io/hdf5_properties_io.h"

namespace Kratos
{
namespace HDF5
{
PartitionedModelPartIO::PartitionedModelPartIO(Parameters& rParams, File::Pointer pFile)
: mpFile(pFile)
{
    KRATOS_TRY;

    Parameters default_params(R"(
        {
            "prefix": "",
            "list_of_elements": [],
            "list_of_conditions": []
        })");

    rParams.ValidateAndAssignDefaults(default_params);

    mPrefix = rParams["prefix"].GetString();

    mElementNames.resize(rParams["list_of_elements"].size());
    for (unsigned i = 0; i < mElementNames.size(); ++i)
        mElementNames[i] = rParams["list_of_elements"].GetArrayItem(i).GetString();

    mConditionNames.resize(rParams["list_of_conditions"].size());
    for (unsigned i = 0; i < mConditionNames.size(); ++i)
        mConditionNames[i] = rParams["list_of_conditions"].GetArrayItem(i).GetString();

    Check();

    mElementPointers.resize(mElementNames.size());
    for (unsigned i = 0; i < mElementNames.size(); ++i)
    {
        const ElementType& r_elem = KratosComponents<Element>::Get(mElementNames[i]);
        mElementPointers[i] = &r_elem;
    }

    mConditionPointers.resize(mConditionNames.size());
    for (unsigned i = 0; i < mConditionNames.size(); ++i)
    {
        const ConditionType& r_cond = KratosComponents<Condition>::Get(mConditionNames[i]);
        mConditionPointers[i] = &r_cond;
    }

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
    Detail::PointsData local_points;
    local_points.ReadData(r_file, mPrefix + "/Nodes/Local", local_start_index, local_block_size);
    local_points.CreateNodes(rNodes);
    local_points.Clear();
    
    // Read ghost nodes.
    Detail::PointsData ghost_points;
    ghost_points.ReadData(r_file, mPrefix + "/Nodes/Ghost", ghost_start_index, ghost_block_size);
    ghost_points.CreateNodes(rNodes);
    ghost_points.Clear();

    return true;
    KRATOS_CATCH("");
}

std::size_t PartitionedModelPartIO::ReadNodesNumber()
{
    const std::vector<unsigned> dims = GetFile().GetDataDimensions("/Nodes/Local/Id");
    return dims[0];
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
    rNodes.Sort(); // Avoid inadvertently reordering partway through the writing process.

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
    Detail::PointsData local_points;
    local_points.SetData(local_nodes);
    local_points.WriteData(r_file, mPrefix + "/Nodes/Local");
    r_file.WriteDataPartition(mPrefix + "/Nodes/Local/Partition", local_points.GetIds());
    local_points.Clear();

    // Write ghost nodes.
    Detail::PointsData ghost_points;
    ghost_points.SetData(ghost_nodes);
    ghost_points.WriteData(r_file, mPrefix + "/Nodes/Ghost");
    r_file.WriteDataPartition(mPrefix + "/Nodes/Ghost/Partition", ghost_points.GetIds());
    ghost_points.Clear();

    KRATOS_CATCH("");
}

void PartitionedModelPartIO::ReadProperties(PropertiesContainerType& rProperties)
{
    PropertiesIO prop_io(mPrefix, mpFile);
    prop_io.ReadProperties(rProperties);
}

void PartitionedModelPartIO::WriteProperties(Properties const& rProperties)
{
    PropertiesIO prop_io(mPrefix, mpFile);
    prop_io.WriteProperties(rProperties);
}

void PartitionedModelPartIO::WriteProperties(PropertiesContainerType const& rProperties)
{
    PropertiesIO prop_io(mPrefix, mpFile);
    prop_io.WriteProperties(rProperties);
}

void PartitionedModelPartIO::ReadElements(NodesContainerType& rNodes,
                                          PropertiesContainerType& rProperties,
                                          ElementsContainerType& rElements)
{
    KRATOS_TRY;

    unsigned start_index, block_size;
    rElements.clear();

    File& r_file = GetFile();

    for (unsigned i = 0; i < mElementNames.size(); ++i)
    {
        const std::string elem_path = mPrefix + "/Elements/" + mElementNames[i];
        Detail::ConnectivitiesData connectivities;
        std::tie(start_index, block_size) =
            GetPartitionStartIndexAndBlockSize(elem_path + "/Partition");
        connectivities.ReadData(r_file, elem_path, start_index, block_size);
        const ElementType& r_elem = *mElementPointers[i];
        connectivities.CreateElements(r_elem, rNodes, rProperties, rElements);
    }

    KRATOS_CATCH("");
}

std::size_t PartitionedModelPartIO::ReadElementsConnectivities(ConnectivitiesContainerType& rElementsConnectivities)
{
    return 0;
}

void PartitionedModelPartIO::WriteElements(ElementsContainerType const& rElements)
{
    KRATOS_TRY;

    const unsigned num_elem_types = mElementNames.size();
    Detail::PointerBinsUtility<ElementType> elem_bins(mElementPointers);
    elem_bins.CreateBins(rElements);
    File& r_file = GetFile();
    for (unsigned i_type = 0; i_type < num_elem_types; ++i_type)
    {
        std::string elem_path = mPrefix + "/Elements/" + mElementNames[i_type];
        const ElementType* elem_key = mElementPointers[i_type];
        ConstElementsContainerType& r_elems = elem_bins.GetBin(elem_key);
        Detail::ConnectivitiesData connectivities;
        connectivities.SetData(r_elems);
        connectivities.WriteData(r_file, elem_path);
        r_file.WriteDataPartition(elem_path + "/Partition", connectivities.GetIds());
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

    for (unsigned i = 0; i < mConditionNames.size(); ++i)
    {
        const std::string cond_path = mPrefix + "/Conditions/" + mConditionNames[i];
        Detail::ConnectivitiesData connectivities;
        std::tie(start_index, block_size) =
            GetPartitionStartIndexAndBlockSize(cond_path + "/Partition");
        connectivities.ReadData(r_file, cond_path, start_index, block_size);
        const ConditionType& r_cond = *mConditionPointers[i];
        connectivities.CreateConditions(r_cond, rNodes, rProperties, rConditions);
    }

    KRATOS_CATCH("");
}

void PartitionedModelPartIO::WriteConditions(ConditionsContainerType const& rConditions)
{
    KRATOS_TRY;

    const unsigned num_cond_types = mConditionNames.size();
    Detail::PointerBinsUtility<ConditionType> cond_bins(mConditionPointers);
    cond_bins.CreateBins(rConditions);
    File& r_file = GetFile();
    for (unsigned i_type = 0; i_type < num_cond_types; ++i_type)
    {
        std::string cond_path = mPrefix + "/Conditions/" + mConditionNames[i_type];
        const ConditionType* cond_key = mConditionPointers[i_type];
        ConstConditionsContainerType& r_conds = cond_bins.GetBin(cond_key);
        Detail::ConnectivitiesData connectivities;
        connectivities.SetData(r_conds);
        connectivities.WriteData(r_file, cond_path);
        r_file.WriteDataPartition(cond_path + "/Partition", connectivities.GetIds());
    }

    KRATOS_CATCH("");
}

std::size_t PartitionedModelPartIO::ReadConditionsConnectivities(ConnectivitiesContainerType& rConditionsConnectivities)
{
    return 0;
}

void PartitionedModelPartIO::ReadInitialValues(ModelPart& rModelPart)
{
}

void PartitionedModelPartIO::ReadInitialValues(NodesContainerType& rNodes,
                                               ElementsContainerType& rElements,
                                               ConditionsContainerType& rConditions)
{
}

void PartitionedModelPartIO::ReadModelPart(ModelPart& rModelPart)
{
    ReadProperties(rModelPart.rProperties());
    ReadNodes(rModelPart.Nodes());
    ReadElements(rModelPart.Nodes(), rModelPart.rProperties(), rModelPart.Elements());
    ReadConditions(rModelPart.Nodes(), rModelPart.rProperties(), rModelPart.Conditions());
}

void PartitionedModelPartIO::WriteModelPart(ModelPart& rModelPart)
{
    WriteProperties(rModelPart.rProperties());
    WriteNodes(rModelPart.Nodes());
    WriteElements(rModelPart.Elements());
    WriteConditions(rModelPart.Conditions());
}

File& PartitionedModelPartIO::GetFile() const
{
    return *mpFile;
}

void PartitionedModelPartIO::Check()
{
    if (GetFile().GetTotalProcesses() == 1)
        std::cout << "Warning: using PartitionedModelPartIO with single process file access." << std::endl;
}

std::tuple<unsigned, unsigned> PartitionedModelPartIO::GetPartitionStartIndexAndBlockSize(std::string Path) const
{
    unsigned my_pid = GetFile().GetPID();
    HDF5::File::Vector<int> my_partition;
    GetFile().ReadDataSet(Path, my_partition, my_pid, 2);
    unsigned start_index = my_partition[0];
    unsigned block_size = my_partition[1] - my_partition[0];
    return std::make_tuple(start_index, block_size);
}

} // namespace HDF5.
} // namespace Kratos.
