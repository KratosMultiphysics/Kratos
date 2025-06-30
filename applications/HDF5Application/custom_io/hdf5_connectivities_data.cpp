//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  license: HDF5Application/license.txt
//
//  Main author:    Michael Andre, https://github.com/msandre
//                  Suneth Warnakulasuriya
//

// System includes

// Project includes
#include "includes/kratos_components.h"
#include "utilities/compare_elements_and_conditions_utility.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "custom_io/hdf5_file.h"
#include "custom_utilities/hdf5_data_set_partition_utility.h"

// Include base h
#include "custom_io/hdf5_connectivities_data.h"

namespace Kratos
{
namespace HDF5
{
namespace Internals
{

template<class TContainerType>
ConnectivitiesData<TContainerType>::ConnectivitiesData(
    const std::string& rPrefix,
    File::Pointer pFile)
    : mpFile(pFile),
      mPrefix(rPrefix)
{
}


template<class TContainerType>
void ConnectivitiesData<TContainerType>::Read(
    const std::string& rEntityName,
    NodesContainerType& rNodes,
    PropertiesContainerType& rProperties,
    TContainerType& rEntities)
{
    KRATOS_TRY;

    KRATOS_ERROR_IF_NOT(KratosComponents<EntityType>::Has(rEntityName))
        << "Entity \"" << rEntityName << "\" not found in registered entities.\n";

    const auto& entity_group_path = mPrefix + "/" + rEntityName;

    KRATOS_ERROR_IF_NOT(mpFile->HasPath(entity_group_path))
        << "Entity data for \"" << rEntityName << "\" not found in file [ expected path in file = "
        << entity_group_path << " ].\n";

    IndexType start_index, block_size;
    std::tie(start_index, block_size) = StartIndexAndBlockSize(*mpFile, entity_group_path);

    Vector<int> entity_ids, entity_property_ids;
    Matrix<int> connectivities;

    std::string file_entity_name;
    mpFile->ReadAttribute(entity_group_path, "Name", file_entity_name);

    KRATOS_ERROR_IF_NOT(file_entity_name == rEntityName)
        << "Entity name mismatch [ expected entity name = "
        << rEntityName << ", entity name in file = " << file_entity_name << " ].\n";

    mpFile->ReadDataSet(entity_group_path + "/Ids", entity_ids, start_index, block_size);
    mpFile->ReadDataSet(entity_group_path + "/PropertiesIds", entity_property_ids, start_index, block_size);
    mpFile->ReadDataSet(entity_group_path + "/Connectivities", connectivities, start_index, block_size);

    const auto& r_entity = KratosComponents<EntityType>::Get(rEntityName);
    const unsigned geometry_size = r_entity.GetGeometry().size();

    KRATOS_ERROR_IF(geometry_size != connectivities.size2())
        << "Entity \"" << rEntityName << "\""
        << " has a different connectivities size than the file data [ expected connectivity size = "
        << geometry_size << ", connectivity size in file = " << connectivities.size2() << " ].\n";

    const unsigned num_new_elems = entity_ids.size();

    rEntities.reserve(rEntities.size() + num_new_elems);
    typename EntityType::NodesArrayType nodes(geometry_size);

    for (unsigned i = 0; i < num_new_elems; ++i) {
        for (unsigned j = 0; j < geometry_size; ++j) {
            const int node_id = connectivities(i, j);
            nodes(j) = rNodes(node_id);
        }
        auto p_elem = r_entity.Create(entity_ids[i], nodes, rProperties(entity_property_ids[i]));
        rEntities.push_back(p_elem);
    }

    KRATOS_CATCH("");
}

template<class TContainerType>
void ConnectivitiesData<TContainerType>::Write(
    const TContainerType& rEntities,
    const bool WriteProperties)
{
    KRATOS_TRY;

    if (mpFile->GetDataCommunicator().SumAll(rEntities.size()) == 0) {
        // do nothing if the all the ranks have no entities.
        return;
    }

    std::string entity_name = "";
    unsigned int geometry_size = 0;
    if (!rEntities.empty()) {
        CompareElementsAndConditionsUtility::GetRegisteredName(rEntities.front(), entity_name);
        // no need of doing communication to get the correct geometry size
        // for all ranks, because the connectivity matrix will be correctly
        // sized in the HDF5File write call by MPI Communication.
        // Otherwise, this will be done twice which is redundant.
        geometry_size = rEntities.front().GetGeometry().size();
    }

    const auto& r_data_communicator = mpFile->GetDataCommunicator();

    // now do the communication between ranks to find the name, because there may be ranks which
    // are empty, but every rank still need the entity name to proceed further.
    std::vector<char> char_entity_name(entity_name.begin(), entity_name.end());
    const auto& global_char_entity_names = r_data_communicator.AllGatherv(char_entity_name);
    if (entity_name.empty()) {
        for (const auto& rank_char_entity_name : global_char_entity_names) {
            if (!rank_char_entity_name.empty()) {
                entity_name = std::string(rank_char_entity_name.begin(), rank_char_entity_name.end());
                break;
            }
        }
    }

    // now cross check all names
    for (const auto& rank_char_entity_name : global_char_entity_names) {
        std::string rank_name(rank_char_entity_name.begin(), rank_char_entity_name.end());
        // here we have to check current entity name is equal to other ranks entity names.
        // but if a rank has not entities, then the AllGatherv above will have an empty string, but
        // in the above if block, we will make the empty rank's entity_name also non emtpy. Therefore
        // an additional check is included to check whether the rank_name is also empty.
        KRATOS_ERROR_IF(entity_name != rank_name && !rank_name.empty()) << "All the ranks should have the same entities.";
    }

    const auto& entity_group_path = mPrefix + "/" + entity_name;

    const unsigned int num_entities = rEntities.size();

    Vector<int> entity_ids, entity_property_ids;
    Matrix<int> connectivities;

    entity_ids.resize(num_entities, false);
    entity_property_ids.resize(num_entities, false);
    connectivities.resize(num_entities, geometry_size, false);

    // Fill arrays and perform checks.
    IndexPartition<IndexType>(num_entities).for_each([&rEntities, &entity_ids, &entity_property_ids, &connectivities, geometry_size](const auto Index) {
        const auto& r_entity = *(rEntities.begin() + Index);

        // Check that the element and geometry types are the same.
        KRATOS_ERROR_IF_NOT(GeometricalObject::IsSame(r_entity, rEntities.front()))
            << "Element #" << r_entity.Id() << " is not the same as #" << rEntities.front().Id() << '!' << std::endl;

        // Fill ids.
        entity_ids[Index] = r_entity.Id();
        entity_property_ids[Index] = r_entity.GetProperties().Id();

        // Fill connectivities.
        const auto& r_geom = r_entity.GetGeometry();
        for (unsigned k = 0; k < geometry_size; ++k) {
            connectivities(Index, k) = r_geom[k].Id();
        }
    });

    int ws_dim, num_nodes;
    if (KratosComponents<EntityType>::Has(entity_name)) {
        const auto& r_geom = KratosComponents<EntityType>::Get(entity_name).GetGeometry();
        ws_dim = r_geom.WorkingSpaceDimension();
        num_nodes = r_geom.size();
    } else {
        KRATOS_ERROR << "Invalid name \"" << entity_name << "\"!" << std::endl;
    }

    WriteInfo info;

    mpFile->WriteDataSet(entity_group_path + "/Ids", entity_ids, info);
    mpFile->WriteDataSet(entity_group_path + "/Connectivities", connectivities, info);
    if (WriteProperties) {
        mpFile->WriteDataSet(entity_group_path + "/PropertiesIds", entity_property_ids, info);
    }

    mpFile->WriteAttribute(entity_group_path, "Name", entity_name);
    mpFile->WriteAttribute(entity_group_path, "WorkingSpaceDimension", ws_dim);
    mpFile->WriteAttribute(entity_group_path, "NumberOfNodes", num_nodes);
    mpFile->WriteAttribute(entity_group_path, "Size", static_cast<int>(info.TotalSize));

    WritePartitionTable(*mpFile, entity_group_path, info);

    KRATOS_CATCH("");
}

// template instantiations
template class KRATOS_API(HDF5_APPLICATION) ConnectivitiesData<ModelPart::ConditionsContainerType>;
template class KRATOS_API(HDF5_APPLICATION) ConnectivitiesData<ModelPart::ElementsContainerType>;

} // namespace Internals.
} // namespace HDF5.
} // namespace Kratos.