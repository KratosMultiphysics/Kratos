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

    const std::size_t num_nodes = ReadNodesNumber();
    Internals::PointsData points;

    points.ReadData(*mpFile, mPrefix + "/Nodes/Local", 0, num_nodes);
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
        int num_elems;
        mpFile->ReadAttribute(mPrefix + "/Elements/" + r_name, "Size", num_elems);
        Internals::ConnectivitiesData connectivities;
        connectivities.ReadData(*mpFile, mPrefix + "/Elements/" + r_name, 0, num_elems);
        connectivities.CreateEntities(rNodes, rProperties, rElements);
    }

    KRATOS_CATCH("");
}

void ModelPartIO::WriteElements(ElementsContainerType const& rElements)
{
    KRATOS_TRY;

    std::vector<ElementsContainerType> factored_elements = FactorElements(rElements);

    WriteInfo info;
    for (const auto& r_elems : factored_elements)
    {
        Internals::ConnectivitiesData connectivities;
        connectivities.SetData(r_elems);
        connectivities.WriteData(*mpFile, mPrefix + "/Elements/" + connectivities.Name(), info);
        const int size = info.TotalSize;
        mpFile->WriteAttribute(mPrefix + "/Elements/" + connectivities.Name(), "Size", size);
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
        int num_conds;
        mpFile->ReadAttribute(mPrefix + "/Conditions/" + r_name, "Size", num_conds);
        Internals::ConnectivitiesData connectivities;
        connectivities.ReadData(*mpFile, mPrefix + "/Conditions/" + r_name, 0, num_conds);
        connectivities.CreateEntities(rNodes, rProperties, rConditions);
    }

    KRATOS_CATCH("");
}

void ModelPartIO::WriteConditions(ConditionsContainerType const& rConditions)
{
    KRATOS_TRY;

    std::vector<ConditionsContainerType> factored_conditions = FactorConditions(rConditions);

    WriteInfo info;
    for (const auto& r_conds : factored_conditions)
    {
        Internals::ConnectivitiesData connectivities;
        connectivities.SetData(r_conds);
        connectivities.WriteData(*mpFile, mPrefix + "/Conditions/" + connectivities.Name(), info);
        const int size = info.TotalSize;
        mpFile->WriteAttribute(mPrefix + "/Conditions/" + connectivities.Name(), "Size", size);
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
