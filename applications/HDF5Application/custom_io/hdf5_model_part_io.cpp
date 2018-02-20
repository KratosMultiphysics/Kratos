#include "custom_io/hdf5_model_part_io.h"

#include "custom_io/hdf5_points_data.h"
#include "custom_io/hdf5_connectivities_data.h"
#include "custom_utilities/factor_elements_and_conditions_utility.h"

namespace Kratos
{
namespace HDF5
{
ModelPartIO::ModelPartIO(File::Pointer pFile, std::string const& rPrefix)
: ModelPartIOBase(pFile, rPrefix)
{
    KRATOS_TRY;

    Check();

    KRATOS_CATCH("");
}

bool ModelPartIO::ReadNodes(NodesContainerType& rNodes)
{
    KRATOS_TRY;

    rNodes.clear();
    int start_index, block_size;
    std::tie(start_index, block_size) =
        StartIndexAndBlockSize(mPrefix + "/Nodes/Local");
    Internals::PointsData points;
    points.ReadData(*mpFile, mPrefix + "/Nodes/Local", start_index, block_size);
    points.CreateNodes(rNodes);
    return true;

    KRATOS_CATCH("");
}

void ModelPartIO::WriteNodes(NodesContainerType const& rNodes)
{
    KRATOS_TRY;

    Internals::PointsData points;
    points.SetData(rNodes);
    WriteInfo info;
    points.WriteData(*mpFile, mPrefix + "/Nodes/Local", info);
    StoreWriteInfo(mPrefix + "/Nodes/Local", info);

    KRATOS_CATCH("");
}

void ModelPartIO::ReadElements(NodesContainerType& rNodes,
                               PropertiesContainerType& rProperties,
                               ElementsContainerType& rElements)
{
    KRATOS_TRY;

    rElements.clear();
    std::vector<std::string> group_names;
    mpFile->GetGroupNames(mPrefix + "/Elements", group_names);

    for (const auto& r_name : group_names)
    {
        int start_index, block_size;
        std::tie(start_index, block_size) = StartIndexAndBlockSize(mPrefix + "/Elements/" + r_name);
        Internals::ConnectivitiesData connectivities;
        connectivities.ReadData(*mpFile, mPrefix + "/Elements/" + r_name, start_index, block_size);
        connectivities.CreateEntities(rNodes, rProperties, rElements);
    }

    KRATOS_CATCH("");
}

void ModelPartIO::WriteElements(ElementsContainerType const& rElements)
{
    KRATOS_TRY;

    WriteInfo info;
    for (const auto& r_elems : FactorElements(rElements))
    {
        Internals::ConnectivitiesData connectivities;
        connectivities.SetData(r_elems);
        connectivities.WriteData(*mpFile, mPrefix + "/Elements/" + connectivities.Name(), info);
        StoreWriteInfo(mPrefix + "/Elements/" + connectivities.Name(), info);
    }

    KRATOS_CATCH("");
}

void ModelPartIO::ReadConditions(NodesContainerType& rNodes,
                                 PropertiesContainerType& rProperties,
                                 ConditionsContainerType& rConditions)
{
    KRATOS_TRY;

    rConditions.clear();
    std::vector<std::string> group_names;
    mpFile->GetGroupNames(mPrefix + "/Conditions", group_names);

    for (const auto& r_name : group_names)
    {
        int start_index, block_size;
        std::tie(start_index, block_size) = StartIndexAndBlockSize(mPrefix + "/Conditions/" + r_name);
        Internals::ConnectivitiesData connectivities;
        connectivities.ReadData(*mpFile, mPrefix + "/Conditions/" + r_name, start_index, block_size);
        connectivities.CreateEntities(rNodes, rProperties, rConditions);
    }

    KRATOS_CATCH("");
}

void ModelPartIO::WriteConditions(ConditionsContainerType const& rConditions)
{
    KRATOS_TRY;

    WriteInfo info;
    for (const auto& r_conds : FactorConditions(rConditions))
    {
        Internals::ConnectivitiesData connectivities;
        connectivities.SetData(r_conds);
        connectivities.WriteData(*mpFile, mPrefix + "/Conditions/" + connectivities.Name(), info);
        StoreWriteInfo(mPrefix + "/Conditions/" + connectivities.Name(), info);
    }

    KRATOS_CATCH("");
}

void ModelPartIO::Check()
{
    KRATOS_TRY;

    KRATOS_ERROR_IF(mpFile->GetTotalProcesses() != 1)
        << "Serial IO expects file access by a single process only." << std::endl;

    KRATOS_CATCH("");
}

} // namespace HDF5.
} // namespace Kratos.
