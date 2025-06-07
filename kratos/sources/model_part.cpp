//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//
//

// System includes
#include <sstream>
#include <type_traits>

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/exception.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

namespace Kratos
{

KRATOS_CREATE_LOCAL_FLAG(ModelPart, ALL_ENTITIES, 0);
KRATOS_CREATE_LOCAL_FLAG(ModelPart, OVERWRITE_ENTITIES, 1);

namespace ModelPartHelperUtilities
{

template<class TContainerGetter>
void AddEntitiesFromIds(
    TContainerGetter&& rContainerGetter,
    ModelPart* pModelPart,
    const std::vector<IndexType>& rEntityIds)
{
    KRATOS_TRY

    using container_type = std::remove_pointer_t<decltype(std::declval<TContainerGetter>()(std::declval<ModelPart*>()))>;

    if(pModelPart->IsSubModelPart()) { //does nothing if we are on the top model part
        //obtain from the root model part the corresponding list of nodes
        ModelPart* root_model_part = &pModelPart->GetRootModelPart();
        const auto& r_container = *rContainerGetter(root_model_part);

        // we first sort and unique the given entity ids.
        std::vector<typename container_type::pointer> entities;
        entities.resize(rEntityIds.size());

        // we are doing this in parallel because the most complex one is the finding of the entity for given id.
        IndexPartition<IndexType>(rEntityIds.size()).for_each([pModelPart, &r_container, &entities,  &rEntityIds](const auto Index) {
            const auto entity_id = rEntityIds[Index];

            auto it = r_container.find(entity_id);

            KRATOS_ERROR_IF(it == r_container.end())
                << "while adding " << ModelPart::Container<container_type>::GetEntityName() << "s to submodelpart "
                << pModelPart->FullName() << ", the "
                << ModelPart::Container<container_type>::GetEntityName()
                << " with Id " << entity_id << " does not exist in the root model part";

            entities[Index] = *(it.base());
        });

        // aux is created to avoid doing sort and unique for all the parent model parts and current model part
        // when insertion is done. since we dont use the entities anymore we can use the move operator here.
        // since aux is an empty container. it will call push_back in the insert.
        container_type aux;
        aux.insert(std::move(entities));

        ModelPart* current_part = pModelPart;
        while(current_part->IsSubModelPart()) {
            // this will directly call the PVS::insert overloaded for the PVS::iterator type.
            rContainerGetter(current_part)->insert(aux.begin(), aux.end());
            current_part = &(current_part->GetParentModelPart());
        }
    }

    KRATOS_CATCH("");
}

template<class TContainerType>
void RemoveEntities(
    ModelPart::MeshType& rMesh,
    const Flags& rIdentifierFlag)
{
    KRATOS_TRY

    auto& r_container = ModelPart::Container<TContainerType>::GetContainer(rMesh);

    //count the nodes to be erase
    const auto erase_count = block_for_each<SumReduction<unsigned int>>(r_container, [&rIdentifierFlag](const auto& rEntity) -> unsigned int {
        return rEntity.Is(rIdentifierFlag);
    });

    TContainerType temp_entities;
    temp_entities.reserve(r_container.size() - erase_count);
    temp_entities.swap(r_container);

    for(auto i_entity = temp_entities.begin() ; i_entity != temp_entities.end() ; ++i_entity) {
        if (i_entity->IsNot(rIdentifierFlag)) {
            // we can safely insert them at the end with the correct hint here since, the original r_mesh
            // is sorted and unique.
            r_container.insert(r_container.end(), std::move(*(i_entity.base())));
        }
    }

    KRATOS_CATCH("");
}

} // namespace ModelPartHelperUtilities

/// Default constructor.
ModelPart::ModelPart(VariablesList::Pointer pVariablesList, Model& rOwnerModel) : ModelPart("Default", pVariablesList, rOwnerModel) { }

/// Constructor with name
ModelPart::ModelPart(std::string const& NewName,VariablesList::Pointer pVariablesList, Model& rOwnerModel) : ModelPart(NewName, 1, pVariablesList, rOwnerModel) { }

/// Constructor with name and bufferSize
ModelPart::ModelPart(std::string const& NewName, IndexType NewBufferSize,VariablesList::Pointer pVariablesList, Model& rOwnerModel)
    : DataValueContainer()
    , Flags()
    , mBufferSize(NewBufferSize)
    , mpProcessInfo(new ProcessInfo())
    , mGeometries()
    , mpVariablesList(pVariablesList)
    , mpCommunicator(new Communicator)
    , mpParentModelPart(NULL)
    , mSubModelParts()
    , mrModel(rOwnerModel)
{
    KRATOS_ERROR_IF(NewName.empty()) << "Please don't use empty names (\"\") when creating a ModelPart" << std::endl;

    KRATOS_ERROR_IF_NOT(NewName.find('.') == std::string::npos) << "Please don't use names containing (\".\") when creating a ModelPart (used in \"" << NewName << "\")" << std::endl;

    mName = NewName;
    MeshType mesh;
    mMeshes.push_back(Kratos::make_shared<MeshType>(mesh.Clone()));
    mpCommunicator->SetLocalMesh(pGetMesh());  // assigning the current mesh to the local mesh of communicator for openmp cases
}

/// Destructor.
ModelPart::~ModelPart()
{
    Clear();
}

void ModelPart::Clear()
{
    KRATOS_TRY

    // Call recursively clear of all submodel parts
    for (auto& r_sub_model_part : mSubModelParts) {
        r_sub_model_part.Clear();
    }

    // Clear sub model parts list
    mSubModelParts.clear();

    // Clear meshes
    for(auto& r_mesh : mMeshes) {
        r_mesh.Clear();
    }

    // Clear meshes list
    mMeshes.clear();
    mMeshes.emplace_back(Kratos::make_shared<MeshType>());

    // Clear geometries
    mGeometries.Clear();

    mTables.clear();

    mpCommunicator->Clear();

    this->AssignFlags(Flags());

    KRATOS_CATCH("");
}

void ModelPart::Reset()
{
    KRATOS_TRY

    // Clears the model part
    Clear();

    // construct a new variable list and process info. Old data ptrs is not destroyed
    // since, same data may be shared with some other model parts as well.
    mpVariablesList = Kratos::make_intrusive<VariablesList>();
    mpProcessInfo = Kratos::make_shared<ProcessInfo>();
    mBufferSize = 0;

    KRATOS_CATCH("");
}

ModelPart::IndexType ModelPart::CreateSolutionStep()
{
    KRATOS_THROW_ERROR(std::logic_error, "This method needs updating and is not working. Pooyan", "")
    return 0;
}

ModelPart::IndexType ModelPart::CloneSolutionStep()
{
    KRATOS_ERROR_IF(IsSubModelPart()) << "Calling the method of the sub model part "
        << Name() << " please call the one of the root model part: "
        << GetRootModelPart().Name() << std::endl;

    auto& r_nodes = Nodes();
    IndexPartition<size_t>(r_nodes.size()).for_each([&](size_t i){
        auto node_iterator = r_nodes.begin() + i;
        node_iterator->CloneSolutionStepData();
    });

    mpProcessInfo->CloneSolutionStepInfo();

    mpProcessInfo->ClearHistory(mBufferSize);

    return 0;
}

ModelPart::IndexType ModelPart::CloneTimeStep()
{
    KRATOS_ERROR_IF(IsSubModelPart()) << "Calling the method of the sub model part "
        << Name() << " please call the one of the root model part: "
        << GetRootModelPart().Name() << std::endl;

    IndexType new_index = CloneSolutionStep();
    mpProcessInfo->SetAsTimeStepInfo();

    return new_index;
}


ModelPart::IndexType ModelPart::CreateTimeStep(double NewTime)
{
    KRATOS_ERROR_IF(IsSubModelPart()) << "Calling the method of the sub model part "
        << Name() << " please call the one of the root model part: "
        << GetRootModelPart().Name() << std::endl;

    IndexType new_index = CreateSolutionStep();
    mpProcessInfo->SetAsTimeStepInfo(NewTime);

    return new_index;
}

ModelPart::IndexType ModelPart::CloneTimeStep(double NewTime)
{
    KRATOS_ERROR_IF(IsSubModelPart()) << "Calling the method of the sub model part "
        << Name() << " please call the one of the root model part: "
        << GetRootModelPart().Name() << std::endl;

    IndexType new_index = CloneSolutionStep();
    mpProcessInfo->SetAsTimeStepInfo(NewTime);

    return new_index;
}

void ModelPart::OverwriteSolutionStepData(IndexType SourceSolutionStepIndex, IndexType DestinationSourceSolutionStepIndex)
{
    KRATOS_ERROR_IF(IsSubModelPart()) << "Calling the method of the sub model part "
        << Name() << " please call the one of the root model part: "
        << GetRootModelPart().Name() << std::endl;

    for (NodeIterator node_iterator = NodesBegin(); node_iterator != NodesEnd(); node_iterator++)
        node_iterator->OverwriteSolutionStepData(SourceSolutionStepIndex, DestinationSourceSolutionStepIndex);

}

void ModelPart::ReduceTimeStep(ModelPart& rModelPart, double NewTime)
{
    KRATOS_TRY

    //ATTENTION: this function does not touch the coordinates of the nodes.
    //It just resets the database values to the values at the beginning of the time step

    KRATOS_ERROR_IF(IsSubModelPart()) << "Calling the method of the sub model part "
        << Name() << " please call the one of the root model part: "
        << GetRootModelPart().Name() << std::endl;

    rModelPart.OverwriteSolutionStepData(1, 0);
    rModelPart.GetProcessInfo().SetCurrentTime(NewTime);

    KRATOS_CATCH("error in reducing the time step")

}

/** Inserts a node in the mesh with ThisIndex.
*/
void ModelPart::AddNode(ModelPart::NodeType::Pointer pNewNode, ModelPart::IndexType ThisIndex)
{
    if (IsSubModelPart())
    {
        mpParentModelPart->AddNode(pNewNode, ThisIndex);
        GetMesh(ThisIndex).AddNode(pNewNode);
    }
    else
    {
        auto existing_node_it = this->GetMesh(ThisIndex).Nodes().find(pNewNode->Id());
        if( existing_node_it == GetMesh(ThisIndex).NodesEnd()) //node did not exist
        {
            GetMesh(ThisIndex).AddNode(pNewNode);
        }
        else //node did exist already
        {
            if(&(*existing_node_it) != (pNewNode.get()))//check if the pointee coincides
                KRATOS_ERROR << "attempting to add pNewNode with Id :" << pNewNode->Id() << ", unfortunately a (different) node with the same Id already exists" << std::endl;
        }
    }
}

/** Inserts a list of nodes in a submodelpart provided their Id. Does nothing if applied to the top model part
*/
void ModelPart::AddNodes(std::vector<IndexType> const& rNodeIds, IndexType ThisIndex)
{
    ModelPartHelperUtilities::AddEntitiesFromIds([](ModelPart* pModelPart) { return &pModelPart->Nodes(); }, this, rNodeIds);
}

/** Inserts a node in the mesh with ThisIndex.
*/
ModelPart::NodeType::Pointer ModelPart::CreateNewNode(int Id, double x, double y, double z, VariablesList::Pointer pNewVariablesList, ModelPart::IndexType ThisIndex)
{
    KRATOS_TRY
    if (IsSubModelPart())
    {
        NodeType::Pointer p_new_node = mpParentModelPart->CreateNewNode(Id, x, y, z, pNewVariablesList, ThisIndex);
        GetMesh(ThisIndex).AddNode(p_new_node);

        return p_new_node;
    }

    //verify if the node exists and eventually give back the existing node
    auto& root_nodes = this->Nodes(); //note that if we are here than we are working with the root model_part
    auto existing_node_it = root_nodes.find(Id);
    if( existing_node_it != root_nodes.end())
    {
        //the node already exists - now check if the position we ask for coincides with the one of the existing one
        double distance = std::sqrt( std::pow( existing_node_it->X() - x,2) + std::pow(existing_node_it->Y() - y,2) + std::pow(existing_node_it->Z() - z,2) );

        KRATOS_ERROR_IF(distance > std::numeric_limits<double>::epsilon()*1000)
            << "trying to create a node with Id " << Id << " however a node with the same Id already exists in the root model part. Existing node coordinates are " << existing_node_it->Coordinates() << " coordinates of the nodes we are attempting to create are :" << x << " " << y << " " << z;

        //if the node we attempt to create is in the same position as the one that is already there, we return the old one
        return *(existing_node_it.base());
    }

    //create a new node
    NodeType::Pointer p_new_node = Kratos::make_intrusive< NodeType >( Id, x, y, z );

    // Giving model part's variables list to the node
    p_new_node->SetSolutionStepVariablesList(pNewVariablesList);

    //set buffer size
    p_new_node->SetBufferSize(mBufferSize);

    //add the new node to the list of nodes
    GetMesh(ThisIndex).AddNode(p_new_node);

    return p_new_node;
    KRATOS_CATCH("")
}

ModelPart::NodeType::Pointer ModelPart::CreateNewNode(ModelPart::IndexType Id, double x, double y, double z, ModelPart::IndexType ThisIndex)
{
    return CreateNewNode(Id, x, y, z, mpVariablesList, ThisIndex);
}

ModelPart::NodeType::Pointer ModelPart::CreateNewNode(ModelPart::IndexType Id, double x, double y, double z, double* pThisData, ModelPart::IndexType ThisIndex)
{
    KRATOS_TRY
    if (IsSubModelPart())

    {
        NodeType::Pointer p_new_node = mpParentModelPart->CreateNewNode(Id, x, y, z, pThisData, ThisIndex);
        GetMesh(ThisIndex).AddNode(p_new_node);

        return p_new_node;
    }
    //verify if the node exists and eventually give back the existing node
    NodesContainerType::iterator existing_node_it = this->GetMesh(ThisIndex).Nodes().find(Id);
    if( existing_node_it != GetMesh(ThisIndex).NodesEnd())
    {
        //the node already exists - now check if the position we ask for coincides with the one of the existing one
        double distance = std::sqrt( std::pow( existing_node_it->X() - x,2) + std::pow(existing_node_it->Y() - y,2) + std::pow(existing_node_it->Z() - z,2) );

        KRATOS_ERROR_IF(distance > std::numeric_limits<double>::epsilon()*1000)
            << "trying to create a node with Id " << Id << " however a node with the same Id already exists in the root model part. Existing node coordinates are " << existing_node_it->Coordinates() << " coordinates of the nodes we are attempting to create are :" << x << " " << y << " " << z;

        //if the node we attempt to create is in the same position as the one that is already there, we return the old one
        return *(existing_node_it.base());
    }

    //create a new node
    NodeType::Pointer p_new_node = Kratos::make_intrusive< NodeType >( Id, x, y, z, mpVariablesList, pThisData, mBufferSize);
    //add the new node to the list of nodes
    GetMesh(ThisIndex).AddNode(p_new_node);

    return p_new_node;
    KRATOS_CATCH("")
}

ModelPart::NodeType::Pointer ModelPart::CreateNewNode(ModelPart::IndexType NodeId, ModelPart::NodeType const& rSourceNode, ModelPart::IndexType ThisIndex)
{
    return CreateNewNode(NodeId, rSourceNode.X(), rSourceNode.Y(), rSourceNode.Z(), mpVariablesList, ThisIndex);
}

void ModelPart::AssignNode(ModelPart::NodeType::Pointer pThisNode, ModelPart::IndexType ThisIndex)
{
    if (IsSubModelPart())
    {
        mpParentModelPart->AssignNode(pThisNode, ThisIndex);

        //add the new node to the list of nodes
        GetMesh(ThisIndex).AddNode(pThisNode);

        return;
    }

    // Giving model part's variables list to the node
    pThisNode->SetSolutionStepVariablesList(mpVariablesList);

    //set buffer size
    pThisNode->SetBufferSize(mBufferSize);

    //add the new node to the list of nodes
    GetMesh(ThisIndex).AddNode(pThisNode);

}

/** Remove the node with given Id from mesh with ThisIndex in this modelpart and all its subs.
*/
void ModelPart::RemoveNode(ModelPart::IndexType NodeId, ModelPart::IndexType ThisIndex)
{
    GetMesh(ThisIndex).RemoveNode(NodeId);

    for (SubModelPartIterator i_sub_model_part = SubModelPartsBegin(); i_sub_model_part != SubModelPartsEnd(); i_sub_model_part++)
        i_sub_model_part->RemoveNode(NodeId, ThisIndex);
}

/** Remove given node from mesh with ThisIndex in this modelpart and all its subs.
*/
void ModelPart::RemoveNode(ModelPart::NodeType& ThisNode, ModelPart::IndexType ThisIndex)
{
    GetMesh(ThisIndex).RemoveNode(ThisNode);

    for (SubModelPartIterator i_sub_model_part = SubModelPartsBegin(); i_sub_model_part != SubModelPartsEnd(); i_sub_model_part++)
        i_sub_model_part->RemoveNode(ThisNode, ThisIndex);
}

/** Remove given node from mesh with ThisIndex in this modelpart and all its subs.
*/
void ModelPart::RemoveNode(ModelPart::NodeType::Pointer pThisNode, ModelPart::IndexType ThisIndex)
{
    GetMesh(ThisIndex).RemoveNode(pThisNode);

    for (SubModelPartIterator i_sub_model_part = SubModelPartsBegin(); i_sub_model_part != SubModelPartsEnd(); i_sub_model_part++)
        i_sub_model_part->RemoveNode(pThisNode, ThisIndex);
}

/** Remove the node with given Id from mesh with ThisIndex in parents, itself and children.
*/
void ModelPart::RemoveNodeFromAllLevels(ModelPart::IndexType NodeId, ModelPart::IndexType ThisIndex)
{
    if (IsSubModelPart())
    {
        mpParentModelPart->RemoveNodeFromAllLevels(NodeId, ThisIndex);
        return;
    }
    RemoveNode(NodeId, ThisIndex);
}

/** Remove given node from mesh with ThisIndex in parents, itself and children.
*/
void ModelPart::RemoveNodeFromAllLevels(ModelPart::NodeType& ThisNode, ModelPart::IndexType ThisIndex)
{
    if (IsSubModelPart())
    {
        mpParentModelPart->RemoveNode(ThisNode, ThisIndex);
        return;
    }
    RemoveNode(ThisNode, ThisIndex);
}

/** Remove given node from mesh with ThisIndex in parents, itself and children.
*/
void ModelPart::RemoveNodeFromAllLevels(ModelPart::NodeType::Pointer pThisNode, ModelPart::IndexType ThisIndex)
{
    if (IsSubModelPart())
    {
        mpParentModelPart->RemoveNode(pThisNode, ThisIndex);
        return;
    }
    RemoveNode(pThisNode, ThisIndex);
}

void ModelPart::RemoveNodes(Flags IdentifierFlag)
{
    // This method is optimized to free the memory
    // Loop over all the local meshes (Is this still necessary with Submodelparts?)
    for(auto& r_mesh: this->GetMeshes()) {
        ModelPartHelperUtilities::RemoveEntities<ModelPart::NodesContainerType>(r_mesh, IdentifierFlag);
    }

    if (IsDistributed()) {
        // Mark the IdentifierFlag across partitions coherently
        this->GetCommunicator().SynchronizeOrNodalFlags(IdentifierFlag);

        // Remove the nodes from the mpi-interfaces in case there is any
        ModelPartHelperUtilities::RemoveEntities<ModelPart::NodesContainerType>(this->GetCommunicator().LocalMesh(), IdentifierFlag);
        for(auto& r_mesh: this->GetCommunicator().LocalMeshes()) {
            ModelPartHelperUtilities::RemoveEntities<ModelPart::NodesContainerType>(r_mesh, IdentifierFlag);
        }

        ModelPartHelperUtilities::RemoveEntities<ModelPart::NodesContainerType>(this->GetCommunicator().GhostMesh(), IdentifierFlag);
        for(auto& r_mesh: this->GetCommunicator().GhostMeshes()) {
            ModelPartHelperUtilities::RemoveEntities<ModelPart::NodesContainerType>(r_mesh, IdentifierFlag);
        }

        ModelPartHelperUtilities::RemoveEntities<ModelPart::NodesContainerType>(this->GetCommunicator().InterfaceMesh(), IdentifierFlag);
        for(auto& r_mesh: this->GetCommunicator().InterfaceMeshes()) {
            ModelPartHelperUtilities::RemoveEntities<ModelPart::NodesContainerType>(r_mesh, IdentifierFlag);
        }
    }

    // Now recursively remove the nodes in the submodelparts
    for (auto& r_sub_model_part : SubModelParts()) {
        r_sub_model_part.RemoveNodes(IdentifierFlag);
    }
}

void ModelPart::RemoveNodesFromAllLevels(Flags IdentifierFlag)
{
    ModelPart& root_model_part = GetRootModelPart();
    root_model_part.RemoveNodes(IdentifierFlag);
}

ModelPart& ModelPart::GetRootModelPart()
{
    if (IsSubModelPart())
        return mpParentModelPart->GetRootModelPart();
    else
        return *this;
}

const ModelPart& ModelPart::GetRootModelPart() const
{
    if (IsSubModelPart())
        return mpParentModelPart->GetRootModelPart();
    else
        return *this;
}

void ModelPart::SetNodalSolutionStepVariablesList()
{
    KRATOS_ERROR_IF(IsSubModelPart()) << "Calling the method of the sub model part "
        << Name() << " please call the one of the root model part: "
        << GetRootModelPart().Name() << std::endl;

    // Iterate over nodes
    auto& r_nodes_array = this->Nodes();
    block_for_each(r_nodes_array,[&](NodeType& rNode) {
        rNode.SetSolutionStepVariablesList(mpVariablesList);
    });
}

/** Inserts a Table
*/
void ModelPart::AddTable(ModelPart::IndexType TableId, ModelPart::TableType::Pointer pNewTable)
{
    if (IsSubModelPart())
        mpParentModelPart->AddTable(TableId, pNewTable);

    mTables.insert(TableId, pNewTable);
}

/** Remove the Table with given Id from current mesh.
*/
void ModelPart::RemoveTable(ModelPart::IndexType TableId)
{
    mTables.erase(TableId);

    for (SubModelPartIterator i_sub_model_part = SubModelPartsBegin(); i_sub_model_part != SubModelPartsEnd(); i_sub_model_part++)
        i_sub_model_part->RemoveTable(TableId);
}

/** Remove the Table with given Id from current mesh in parents, itself and all children.
*/
void ModelPart::RemoveTableFromAllLevels(ModelPart::IndexType TableId)
{
    if (IsSubModelPart())
    {
        mpParentModelPart->RemoveTableFromAllLevels(TableId);
        return;
    }

    RemoveTable(TableId);
}

/***********************************************************************************/
/***********************************************************************************/

ModelPart::SizeType ModelPart::NumberOfProperties(IndexType ThisIndex) const
{
    return GetMesh(ThisIndex).NumberOfProperties();
}

/** Inserts a properties in the mesh with ThisIndex.
*/
void ModelPart::AddProperties(ModelPart::PropertiesType::Pointer pNewProperties, ModelPart::IndexType ThisIndex)
{
    if (IsSubModelPart())
    {
        mpParentModelPart->AddProperties(pNewProperties, ThisIndex);
    }

    auto existing_prop_it = GetMesh(ThisIndex).Properties().find(pNewProperties->Id());
    if( existing_prop_it != GetMesh(ThisIndex).Properties().end() )
    {
        KRATOS_ERROR_IF( &(*existing_prop_it) != pNewProperties.get() )
            << "trying to add a property with existing Id within the model part : " << Name() << ", property Id is :" << pNewProperties->Id();
    }
    else
    {
        GetMesh(ThisIndex).AddProperties(pNewProperties);
    }
}

/***********************************************************************************/
/***********************************************************************************/

bool ModelPart::HasProperties(
    IndexType PropertiesId,
    IndexType MeshIndex
    ) const
{
    auto pprop_it = GetMesh(MeshIndex).Properties().find(PropertiesId);
    if(pprop_it != GetMesh(MeshIndex).Properties().end()) { // Property does exist
        return true;
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

bool ModelPart::RecursivelyHasProperties(
    IndexType PropertiesId,
    IndexType MeshIndex
    ) const
{
    auto pprop_it = GetMesh(MeshIndex).Properties().find(PropertiesId);
    if(pprop_it != GetMesh(MeshIndex).Properties().end()) { // Property does exist
        return true;
    } else {
        if(IsSubModelPart()) {
            return mpParentModelPart->RecursivelyHasProperties(PropertiesId, MeshIndex);
        } else {
            return false;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

ModelPart::PropertiesType::Pointer ModelPart::CreateNewProperties(
    IndexType PropertiesId,
    IndexType MeshIndex
    )
{
    auto pprop_it = GetMesh(MeshIndex).Properties().find(PropertiesId);
    if(pprop_it != GetMesh(MeshIndex).Properties().end()) { // Property does exist
        KRATOS_ERROR << "Property #" << PropertiesId << " already existing. Please use pGetProperties() instead" << std::endl;
    } else {
        if(IsSubModelPart()) {
            PropertiesType::Pointer pprop =  mpParentModelPart->CreateNewProperties(PropertiesId, MeshIndex);
            GetMesh(MeshIndex).AddProperties(pprop);
            return pprop;
        } else {
            PropertiesType::Pointer pnew_property = Kratos::make_shared<PropertiesType>(PropertiesId);
            GetMesh(MeshIndex).AddProperties(pnew_property);
            return pnew_property;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

ModelPart::PropertiesType::Pointer ModelPart::pGetProperties(
    IndexType PropertiesId,
    IndexType MeshIndex
    )
{
    auto pprop_it = GetMesh(MeshIndex).Properties().find(PropertiesId);
    if(pprop_it != GetMesh(MeshIndex).Properties().end()) { // Property does exist
        return *(pprop_it.base());
    } else {
        if(IsSubModelPart()) {
            PropertiesType::Pointer pprop =  mpParentModelPart->pGetProperties(PropertiesId, MeshIndex);
            GetMesh(MeshIndex).AddProperties(pprop);
            return pprop;
        } else {
            KRATOS_WARNING("ModelPart") << "Property " << PropertiesId << " does not exist!. Creating and adding new property. Please use CreateNewProperties() instead" << std::endl;
            PropertiesType::Pointer pnew_property = Kratos::make_shared<PropertiesType>(PropertiesId);
            GetMesh(MeshIndex).AddProperties(pnew_property);
            return pnew_property;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

const ModelPart::PropertiesType::Pointer ModelPart::pGetProperties(
    IndexType PropertiesId,
    IndexType MeshIndex
    ) const
{
    auto pprop_it = GetMesh(MeshIndex).Properties().find(PropertiesId);
    if(pprop_it != GetMesh(MeshIndex).Properties().end()) { // Property does exist
        return *(pprop_it.base());
    } else {
        if(IsSubModelPart()) {
            PropertiesType::Pointer pprop =  mpParentModelPart->pGetProperties(PropertiesId, MeshIndex);
            return pprop;
        } else {
            KRATOS_ERROR << "Property " << PropertiesId << " does not exist!. This is constant model part and cannot be created a new one" << std::endl;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

ModelPart::PropertiesType& ModelPart::GetProperties(
    IndexType PropertiesId,
    IndexType MeshIndex
    )
{
    auto pprop_it = GetMesh(MeshIndex).Properties().find(PropertiesId);
    if(pprop_it != GetMesh(MeshIndex).Properties().end()) { // Property does exist
        return *pprop_it;
    } else {
        if(IsSubModelPart()) {
            PropertiesType::Pointer pprop =  mpParentModelPart->pGetProperties(PropertiesId, MeshIndex);
            GetMesh(MeshIndex).AddProperties(pprop);
            return *pprop;
        } else {
            KRATOS_WARNING("ModelPart") << "Property " << PropertiesId << " does not exist!. Creating and adding new property. Please use CreateNewProperties() instead" << std::endl;
            PropertiesType::Pointer pnew_property = Kratos::make_shared<PropertiesType>(PropertiesId);
            GetMesh(MeshIndex).AddProperties(pnew_property);
            return *pnew_property;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

const ModelPart::PropertiesType& ModelPart::GetProperties(
    IndexType PropertiesId,
    IndexType MeshIndex
    ) const
{
    auto pprop_it = GetMesh(MeshIndex).Properties().find(PropertiesId);
    if(pprop_it != GetMesh(MeshIndex).Properties().end()) { // Property does exist
        return *pprop_it;
    } else {
        if(IsSubModelPart()) {
            PropertiesType::Pointer pprop =  mpParentModelPart->pGetProperties(PropertiesId, MeshIndex);
            return *pprop;
        } else {
            KRATOS_ERROR << "Property " << PropertiesId << " does not exist!. This is constant model part and cannot be created a new one" << std::endl;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

bool ModelPart::HasProperties(
    const std::string& rAddress,
    IndexType MeshIndex
    ) const
{
    const std::vector<IndexType> component_name = TrimComponentName(rAddress);
    if (HasProperties(component_name[0], MeshIndex)) {
        bool has_properties = true;
        Properties::Pointer p_prop = pGetProperties(component_name[0], MeshIndex);
        for (IndexType i = 1; i < component_name.size(); ++i) {
            if (p_prop->HasSubProperties(component_name[i])) {
                p_prop = p_prop->pGetSubProperties(component_name[i]);
            } else {
                return false;
            }
        }
        return has_properties;
    } else {
        return false;
    }
}

/***********************************************************************************/
/***********************************************************************************/

Properties::Pointer ModelPart::pGetProperties(
    const std::string& rAddress,
    IndexType MeshIndex
    )
{
    const std::vector<IndexType> component_name = TrimComponentName(rAddress);
    if (HasProperties(component_name[0], MeshIndex)) {
        Properties::Pointer p_prop = pGetProperties(component_name[0], MeshIndex);
        for (IndexType i = 1; i < component_name.size(); ++i) {
            if (p_prop->HasSubProperties(component_name[i])) {
                p_prop = p_prop->pGetSubProperties(component_name[i]);
            } else {
                KRATOS_ERROR << "Index is wrong, does not correspond with any sub Properties Id: " << rAddress << std::endl;
            }
        }
        return p_prop;
    } else {
        KRATOS_ERROR << "First index is wrong, does not correspond with any sub Properties Id: " << component_name[0] << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

const Properties::Pointer ModelPart::pGetProperties(
    const std::string& rAddress,
    IndexType MeshIndex
    ) const
{
    const std::vector<IndexType> component_name = TrimComponentName(rAddress);
    if (HasProperties(component_name[0], MeshIndex)) {
        Properties::Pointer p_prop = pGetProperties(component_name[0], MeshIndex);
        for (IndexType i = 1; i < component_name.size(); ++i) {
            if (p_prop->HasSubProperties(component_name[i])) {
                p_prop = p_prop->pGetSubProperties(component_name[i]);
            } else {
                KRATOS_ERROR << "Index is wrong, does not correspond with any sub Properties Id: " << rAddress << std::endl;
            }
        }
        return p_prop;
    } else {
        KRATOS_ERROR << "First index is wrong, does not correspond with any sub Properties Id: " << component_name[0] << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

Properties& ModelPart::GetProperties(
    const std::string& rAddress,
    IndexType MeshIndex
    )
{
    return *pGetProperties(rAddress, MeshIndex);
}

/***********************************************************************************/
/***********************************************************************************/

const Properties& ModelPart::GetProperties(
    const std::string& rAddress,
    IndexType MeshIndex
    ) const
{
    return *pGetProperties(rAddress, MeshIndex);
}

/** Remove the Properties with given Id from mesh with ThisIndex in this modelpart and all its subs.
*/
void ModelPart::RemoveProperties(ModelPart::IndexType PropertiesId, IndexType ThisIndex)
{
    GetMesh(ThisIndex).RemoveProperties(PropertiesId);

    for (SubModelPartIterator i_sub_model_part = SubModelPartsBegin(); i_sub_model_part != SubModelPartsEnd(); i_sub_model_part++)
        i_sub_model_part->RemoveProperties(PropertiesId, ThisIndex);
}

/** Remove given Properties from mesh with ThisIndex in this modelpart and all its subs.
*/
void ModelPart::RemoveProperties(ModelPart::PropertiesType& ThisProperties, ModelPart::IndexType ThisIndex)
{
    GetMesh(ThisIndex).RemoveProperties(ThisProperties);

    for (SubModelPartIterator i_sub_model_part = SubModelPartsBegin(); i_sub_model_part != SubModelPartsEnd(); i_sub_model_part++)
        i_sub_model_part->RemoveProperties(ThisProperties, ThisIndex);
}

/** Remove given Properties from mesh with ThisIndex in this modelpart and all its subs.
*/
void ModelPart::RemoveProperties(ModelPart::PropertiesType::Pointer pThisProperties, ModelPart::IndexType ThisIndex)
{
    GetMesh(ThisIndex).RemoveProperties(pThisProperties);

    for (SubModelPartIterator i_sub_model_part = SubModelPartsBegin(); i_sub_model_part != SubModelPartsEnd(); i_sub_model_part++)
        i_sub_model_part->RemoveProperties(pThisProperties, ThisIndex);
}

/** Remove the Properties with given Id from mesh with ThisIndex in parents, itself and children.
*/
void ModelPart::RemovePropertiesFromAllLevels(ModelPart::IndexType PropertiesId, IndexType ThisIndex)
{
    if (IsSubModelPart())
    {
        mpParentModelPart->RemovePropertiesFromAllLevels(PropertiesId, ThisIndex);
        return;
    }

    RemoveProperties(PropertiesId, ThisIndex);
}

/** Remove given Properties from mesh with ThisIndex in parents, itself and children.
*/
void ModelPart::RemovePropertiesFromAllLevels(ModelPart::PropertiesType& ThisProperties, ModelPart::IndexType ThisIndex)
{
    if (IsSubModelPart())
    {
        mpParentModelPart->RemoveProperties(ThisProperties, ThisIndex);
    }

    RemoveProperties(ThisProperties, ThisIndex);
}

/** Remove given Properties from mesh with ThisIndex in parents, itself and children.
*/
void ModelPart::RemovePropertiesFromAllLevels(ModelPart::PropertiesType::Pointer pThisProperties, ModelPart::IndexType ThisIndex)
{
    if (IsSubModelPart())
    {
        mpParentModelPart->RemoveProperties(pThisProperties, ThisIndex);
    }

    RemoveProperties(pThisProperties, ThisIndex);
}

/** Inserts a element in the mesh with ThisIndex.
*/
void ModelPart::AddElement(ModelPart::ElementType::Pointer pNewElement, ModelPart::IndexType ThisIndex)
{
    if (IsSubModelPart())
    {
        mpParentModelPart->AddElement(pNewElement, ThisIndex);
        GetMesh(ThisIndex).AddElement(pNewElement);
    }
    else
    {
        auto existing_element_it = this->GetMesh(ThisIndex).Elements().find(pNewElement->Id());
        if( existing_element_it == GetMesh(ThisIndex).ElementsEnd()) //node did not exist
        {
            GetMesh(ThisIndex).AddElement(pNewElement);
        }
        else //node did exist already
        {
            KRATOS_ERROR_IF(&(*existing_element_it) != (pNewElement.get()))//check if the pointee coincides
                << "attempting to add pNewElement with Id :" << pNewElement->Id() << ", unfortunately a (different) element with the same Id already exists" << std::endl;
        }
    }
}

/** Inserts a list of conditions to a submodelpart provided their Id. Does nothing if applied to the top model part
*/
void ModelPart::AddElements(std::vector<IndexType> const& ElementIds, IndexType ThisIndex)
{
    ModelPartHelperUtilities::AddEntitiesFromIds([](ModelPart* pModelPart) { return &pModelPart->Elements(); }, this, ElementIds);
}

/** Inserts an element in the mesh with ThisIndex.
*/
ModelPart::ElementType::Pointer ModelPart::CreateNewElement(std::string ElementName,
        ModelPart::IndexType Id, std::vector<ModelPart::IndexType> ElementNodeIds,
        ModelPart::PropertiesType::Pointer pProperties, ModelPart::IndexType ThisIndex)
{
    KRATOS_TRY
    if (IsSubModelPart())
    {
        ElementType::Pointer p_new_element = mpParentModelPart->CreateNewElement(ElementName, Id, ElementNodeIds, pProperties, ThisIndex);
        GetMesh(ThisIndex).AddElement(p_new_element);
        return p_new_element;
    }

    Geometry< Node >::PointsArrayType pElementNodes;

    for (unsigned int i = 0; i < ElementNodeIds.size(); i++)
    {
        pElementNodes.push_back(pGetNode(ElementNodeIds[i]));
    }

    return CreateNewElement(ElementName, Id, pElementNodes, pProperties, ThisIndex);
    KRATOS_CATCH("");
}

/** Inserts an element in the mesh with ThisIndex.
*/
ModelPart::ElementType::Pointer ModelPart::CreateNewElement(std::string ElementName,
        ModelPart::IndexType Id, Geometry< Node >::PointsArrayType pElementNodes,
        ModelPart::PropertiesType::Pointer pProperties, ModelPart::IndexType ThisIndex)
{
    KRATOS_TRY
    if (IsSubModelPart())
    {
        ElementType::Pointer p_new_element = mpParentModelPart->CreateNewElement(ElementName, Id, pElementNodes, pProperties, ThisIndex);
        GetMesh(ThisIndex).AddElement(p_new_element);
        return p_new_element;
    }

    auto existing_element_iterator = GetMesh(ThisIndex).Elements().find(Id);
    KRATOS_ERROR_IF(existing_element_iterator != GetMesh(ThisIndex).ElementsEnd() )
        << "trying to construct an element with ID " << Id << " however an element with the same Id already exists";


    //create the new element
    ElementType const& r_clone_element = KratosComponents<ElementType>::Get(ElementName);
    Element::Pointer p_element = r_clone_element.Create(Id, pElementNodes, pProperties);

    //add the new element
    GetMesh(ThisIndex).AddElement(p_element);

    return p_element;
    KRATOS_CATCH("")
}

/** Inserts an element in the mesh with ThisIndex.
*/
ModelPart::ElementType::Pointer ModelPart::CreateNewElement(std::string ElementName,
        ModelPart::IndexType Id, typename GeometryType::Pointer pGeometry,
        ModelPart::PropertiesType::Pointer pProperties, ModelPart::IndexType ThisIndex)
{
    KRATOS_TRY
    if (IsSubModelPart())
    {
        ElementType::Pointer p_new_element = mpParentModelPart->CreateNewElement(ElementName, Id, pGeometry, pProperties, ThisIndex);
        GetMesh(ThisIndex).AddElement(p_new_element);
        return p_new_element;
    }

    auto existing_element_iterator = GetMesh(ThisIndex).Elements().find(Id);
    KRATOS_ERROR_IF(existing_element_iterator != GetMesh(ThisIndex).ElementsEnd() )
        << "trying to construct an element with ID " << Id << " however an element with the same Id already exists";


    //create the new element
    ElementType const& r_clone_element = KratosComponents<ElementType>::Get(ElementName);
    Element::Pointer p_element = r_clone_element.Create(Id, pGeometry, pProperties);

    //add the new element
    GetMesh(ThisIndex).AddElement(p_element);

    return p_element;
    KRATOS_CATCH("")
}

/** Remove the element with given Id from mesh with ThisIndex in this modelpart and all its subs.
*/
void ModelPart::RemoveElement(ModelPart::IndexType ElementId, ModelPart::IndexType ThisIndex)
{
    GetMesh(ThisIndex).RemoveElement(ElementId);

    for (SubModelPartIterator i_sub_model_part = SubModelPartsBegin(); i_sub_model_part != SubModelPartsEnd(); i_sub_model_part++)
        i_sub_model_part->RemoveElement(ElementId, ThisIndex);
}

/** Remove given element from mesh with ThisIndex in this modelpart and all its subs.
*/
void ModelPart::RemoveElement(ModelPart::ElementType& ThisElement, ModelPart::IndexType ThisIndex)
{
    GetMesh(ThisIndex).RemoveElement(ThisElement);

    for (SubModelPartIterator i_sub_model_part = SubModelPartsBegin(); i_sub_model_part != SubModelPartsEnd(); i_sub_model_part++)
        i_sub_model_part->RemoveElement(ThisElement, ThisIndex);
}

/** Remove given element from mesh with ThisIndex in this modelpart and all its subs.
*/
void ModelPart::RemoveElement(ModelPart::ElementType::Pointer pThisElement, ModelPart::IndexType ThisIndex)
{
    GetMesh(ThisIndex).RemoveElement(pThisElement);

    for (SubModelPartIterator i_sub_model_part = SubModelPartsBegin(); i_sub_model_part != SubModelPartsEnd(); i_sub_model_part++)
        i_sub_model_part->RemoveElement(pThisElement, ThisIndex);
}

/** Remove the element with given Id from mesh with ThisIndex in parents, itself and children.
*/
void ModelPart::RemoveElementFromAllLevels(ModelPart::IndexType ElementId, ModelPart::IndexType ThisIndex)
{
    if (IsSubModelPart())
    {
        mpParentModelPart->RemoveElement(ElementId, ThisIndex);
        return;
    }

    RemoveElement(ElementId, ThisIndex);
}

/** Remove given element from mesh with ThisIndex in parents, itself and children.
*/
void ModelPart::RemoveElementFromAllLevels(ModelPart::ElementType& ThisElement, ModelPart::IndexType ThisIndex)
{
    if (IsSubModelPart())
    {
        mpParentModelPart->RemoveElement(ThisElement, ThisIndex);
        return;
    }

    RemoveElement(ThisElement, ThisIndex);
}

/** Remove given element from mesh with ThisIndex in parents, itself and children.
*/
void ModelPart::RemoveElementFromAllLevels(ModelPart::ElementType::Pointer pThisElement, ModelPart::IndexType ThisIndex)
{
    if (IsSubModelPart())
    {
        mpParentModelPart->RemoveElement(pThisElement, ThisIndex);
        return;
    }

    RemoveElement(pThisElement, ThisIndex);
}

void ModelPart::RemoveElements(Flags IdentifierFlag)
{
    // This method is optimized to free the memory
    //loop over all the meshes
    for(auto& r_mesh : this->GetMeshes()) {
        ModelPartHelperUtilities::RemoveEntities<ModelPart::ElementsContainerType>(r_mesh, IdentifierFlag);
    }

    //now recursively remove the elements in the submodelparts
    for (auto& r_sub_model_part : this->SubModelParts()) {
        r_sub_model_part.RemoveElements(IdentifierFlag);
    }
}

void ModelPart::RemoveElementsFromAllLevels(Flags IdentifierFlag)
{
    ModelPart& root_model_part = GetRootModelPart();
    root_model_part.RemoveElements(IdentifierFlag);
}

/*
    Functions for Master-Slave Constraint
*/

/** Inserts a master-slave constraint in the current mesh.
 */
void ModelPart::AddMasterSlaveConstraint(ModelPart::MasterSlaveConstraintType::Pointer pNewMasterSlaveConstraint, IndexType ThisIndex)
{
    if (IsSubModelPart())
    {
        // First add it to the parent modelpart
        mpParentModelPart->AddMasterSlaveConstraint(pNewMasterSlaveConstraint, ThisIndex);
        GetMesh(ThisIndex).AddMasterSlaveConstraint(pNewMasterSlaveConstraint);
    }
    else
    {
        auto existing_constraint_it = GetMesh(ThisIndex).MasterSlaveConstraints().find(pNewMasterSlaveConstraint->Id());
        if( existing_constraint_it == GetMesh(ThisIndex).MasterSlaveConstraintsEnd()) //master-slave constraint did not exist
        {
            GetMesh(ThisIndex).AddMasterSlaveConstraint(pNewMasterSlaveConstraint);
        }
        else //master-slave constraint did exist already
        {
            KRATOS_ERROR_IF(&(*existing_constraint_it) != (pNewMasterSlaveConstraint.get()))//check if the pointee coincides
                << "attempting to add Master-Slave constraint with Id :" << pNewMasterSlaveConstraint->Id() << ", unfortunately a (different) condition with the same Id already exists" << std::endl;
        }
    }
}

/** Inserts a list of master-slave constraints to a submodelpart provided their Id. Does nothing if applied to the top model part
 */
void ModelPart::AddMasterSlaveConstraints(std::vector<IndexType> const& MasterSlaveConstraintIds, IndexType ThisIndex)
{
    ModelPartHelperUtilities::AddEntitiesFromIds([](ModelPart* pModelPart) { return &pModelPart->MasterSlaveConstraints(); }, this, MasterSlaveConstraintIds);
}

/// @brief Construct a new @ref MasterSlaveConstraint and insert it into the specified @ref Mesh.
/// @note The constraint is created by the root @ref ModelPart and inserted into the root mesh as well.
/// @throws if a constraint with the same ID already exists in the target mesh.
ModelPart::MasterSlaveConstraintType::Pointer ModelPart::CreateNewMasterSlaveConstraint(const std::string& ConstraintName,
                                                                                        IndexType Id,
                                                                                        ModelPart::DofsVectorType& rMasterDofsVector,
                                                                                        ModelPart::DofsVectorType& rSlaveDofsVector,
                                                                                        const ModelPart::MatrixType& RelationMatrix,
                                                                                        const ModelPart::VectorType& ConstantVector,
                                                                                        IndexType ThisIndex)
{

    KRATOS_TRY
    MeshType& r_mesh = GetMesh(ThisIndex);
    ModelPart::MasterSlaveConstraintType::Pointer p_new_constraint;

    if (IsSubModelPart()) {
        // Defer constraint construction to the root model part
        p_new_constraint = mpParentModelPart->CreateNewMasterSlaveConstraint(
            ConstraintName,
            Id,
            rMasterDofsVector,
            rSlaveDofsVector,
            RelationMatrix,
            ConstantVector,
            ThisIndex);

        // Add the constraint
        if (&r_mesh != &mpParentModelPart->GetMesh(ThisIndex)) {
            KRATOS_ERROR_IF_NOT(r_mesh.AddMasterSlaveConstraint(p_new_constraint))
                << "trying to insert a master-slave constraint with ID "
                << Id << " but a constraint with the same ID already exists\n";
        }
    } else /*IsSubModelPart*/ {
        // Construct the new constraint
        ModelPart::MasterSlaveConstraintType const& r_registered_constraint = KratosComponents<MasterSlaveConstraintType>::Get(ConstraintName);
        p_new_constraint = r_registered_constraint.Create(
            Id,
            rMasterDofsVector,
            rSlaveDofsVector,
            RelationMatrix,
            ConstantVector);

        KRATOS_ERROR_IF_NOT(r_mesh.AddMasterSlaveConstraint(p_new_constraint))
            << "trying to insert a master-slave constraint with ID "
            << Id << " but a constraint with the same ID already exists\n";
    }

    return p_new_constraint;
    KRATOS_CATCH("")

}

ModelPart::MasterSlaveConstraintType::Pointer ModelPart::CreateNewMasterSlaveConstraint(const std::string& ConstraintName,
                                                                                    ModelPart::IndexType Id,
                                                                                    ModelPart::NodeType& rMasterNode,
                                                                                    const ModelPart::DoubleVariableType& rMasterVariable,
                                                                                    ModelPart::NodeType& rSlaveNode,
                                                                                    const ModelPart::DoubleVariableType& rSlaveVariable,
                                                                                    const double Weight,
                                                                                    const double Constant,
                                                                                    IndexType ThisIndex)
{
    KRATOS_ERROR_IF_NOT(rMasterNode.HasDofFor(rMasterVariable))
        << "master node " << rMasterNode.Id() << " has no variable " << rMasterVariable.Name() << "\n";

    KRATOS_ERROR_IF_NOT(rSlaveNode.HasDofFor(rSlaveVariable))
        << "slave node " << rSlaveNode.Id() << " has no variable " << rSlaveVariable.Name() << "\n";

    KRATOS_TRY

    ModelPart::MasterSlaveConstraintType::Pointer p_new_constraint;
    MeshType& r_mesh = this->GetMesh(ThisIndex);

    if (IsSubModelPart()) {
        // Defer constraint construction to the root model part
        p_new_constraint = mpParentModelPart->CreateNewMasterSlaveConstraint(
            ConstraintName,
            Id,
            rMasterNode,
            rMasterVariable,
            rSlaveNode,
            rSlaveVariable,
            Weight,
            Constant,
            ThisIndex);

        // Insert the constraint
        if (&r_mesh != &mpParentModelPart->GetMesh(ThisIndex)) {
            KRATOS_ERROR_IF_NOT(r_mesh.AddMasterSlaveConstraint(p_new_constraint))
                << "trying to insert a master-slave constraint with ID "
                << Id << " but a constraint with the same ID already exists\n";
        }
    } else { /*IsSubModelPart*/
        // Construct the new constraint
        ModelPart::MasterSlaveConstraintType const& r_clone_constraint = KratosComponents<MasterSlaveConstraintType>::Get(ConstraintName);
        p_new_constraint = r_clone_constraint.Create(
            Id,
            rMasterNode,
            rMasterVariable,
            rSlaveNode,
            rSlaveVariable,
            Weight,
            Constant);

        // Insert the constraint
        KRATOS_ERROR_IF_NOT(r_mesh.AddMasterSlaveConstraint(p_new_constraint))
            << "trying to insert a master-slave constraint with ID "
            << Id << " but a constraint with the same ID already exists\n";
    }

    return p_new_constraint;

    KRATOS_CATCH("")

}

/** Remove the master-slave constraint with given Id from this modelpart and all its subs.
*/
void ModelPart::RemoveMasterSlaveConstraint(ModelPart::IndexType MasterSlaveConstraintId,  IndexType ThisIndex)
{
    GetMesh(ThisIndex).RemoveMasterSlaveConstraint(MasterSlaveConstraintId);

    for (SubModelPartIterator i_sub_model_part = SubModelPartsBegin(); i_sub_model_part != SubModelPartsEnd(); i_sub_model_part++)
        i_sub_model_part->RemoveMasterSlaveConstraint(MasterSlaveConstraintId, ThisIndex);
}

/** Remove given master-slave constraint from this modelpart and all its subs.
*/
void ModelPart::RemoveMasterSlaveConstraint(ModelPart::MasterSlaveConstraintType& ThisMasterSlaveConstraint, IndexType ThisIndex)
{
    GetMesh(ThisIndex).RemoveMasterSlaveConstraint(ThisMasterSlaveConstraint);
    for (SubModelPartIterator i_sub_model_part = SubModelPartsBegin(); i_sub_model_part != SubModelPartsEnd(); i_sub_model_part++)
        i_sub_model_part->RemoveMasterSlaveConstraint(ThisMasterSlaveConstraint, ThisIndex);
}

/** Remove the master-slave constraint with given Id from mesh with ThisIndex in parents, itself and children.
*/
void ModelPart::RemoveMasterSlaveConstraintFromAllLevels(ModelPart::IndexType MasterSlaveConstraintId, IndexType ThisIndex)
{

    if (IsSubModelPart()){
        mpParentModelPart->RemoveMasterSlaveConstraintFromAllLevels(MasterSlaveConstraintId, ThisIndex);
    }
    RemoveMasterSlaveConstraint(MasterSlaveConstraintId, ThisIndex);

}

/** Remove given master-slave constraint from mesh with ThisIndex in parents, itself and children.
*/
void ModelPart::RemoveMasterSlaveConstraintFromAllLevels(ModelPart::MasterSlaveConstraintType& ThisMasterSlaveConstraint, IndexType ThisIndex)
{

    if (IsSubModelPart()){
        mpParentModelPart->RemoveMasterSlaveConstraintFromAllLevels(ThisMasterSlaveConstraint, ThisIndex);
    }
    RemoveMasterSlaveConstraint(ThisMasterSlaveConstraint, ThisIndex);
}

/***********************************************************************************/
/***********************************************************************************/

void ModelPart::RemoveMasterSlaveConstraints(Flags IdentifierFlag)
{
    // This method is optimized to free the memory loop over all the meshes
    for(auto& r_mesh : this->GetMeshes()) {
        ModelPartHelperUtilities::RemoveEntities<ModelPart::MasterSlaveConstraintContainerType>(r_mesh, IdentifierFlag);
    }

    // Now recursively remove the constraints in the submodelparts
    for (auto& r_sub_model_part : this->SubModelParts()) {
        r_sub_model_part.RemoveMasterSlaveConstraints(IdentifierFlag);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void ModelPart::RemoveMasterSlaveConstraintsFromAllLevels(Flags IdentifierFlag)
{
    ModelPart& root_model_part = GetRootModelPart();
    root_model_part.RemoveMasterSlaveConstraints(IdentifierFlag);
}

/** Returns the MasterSlaveConstraint::Pointer  corresponding to it's identifier */
ModelPart::MasterSlaveConstraintType::Pointer ModelPart::pGetMasterSlaveConstraint(ModelPart::IndexType MasterSlaveConstraintId, IndexType ThisIndex)
{
        return GetMesh(ThisIndex).pGetMasterSlaveConstraint(MasterSlaveConstraintId);
}

/** Returns a reference MasterSlaveConstraint corresponding to it's identifier */
ModelPart::MasterSlaveConstraintType& ModelPart::GetMasterSlaveConstraint(ModelPart::IndexType MasterSlaveConstraintId,  IndexType ThisIndex)
{
        return GetMesh(ThisIndex).GetMasterSlaveConstraint(MasterSlaveConstraintId);
}

const ModelPart::MasterSlaveConstraintType& ModelPart::GetMasterSlaveConstraint(ModelPart::IndexType MasterSlaveConstraintId, IndexType ThisIndex) const
{
    return GetMesh(ThisIndex).GetMasterSlaveConstraint(MasterSlaveConstraintId);
}

/** Inserts a condition in the mesh with ThisIndex.
*/
void ModelPart::AddCondition(ModelPart::ConditionType::Pointer pNewCondition, ModelPart::IndexType ThisIndex)
{
    if (IsSubModelPart())
    {
        mpParentModelPart->AddCondition(pNewCondition, ThisIndex);
        GetMesh(ThisIndex).AddCondition(pNewCondition);
    }
    else
    {
        auto existing_condition_it = this->GetMesh(ThisIndex).Conditions().find(pNewCondition->Id());
        if( existing_condition_it == GetMesh(ThisIndex).ConditionsEnd()) //node did not exist
        {
            GetMesh(ThisIndex).AddCondition(pNewCondition);
        }
        else //node did exist already
        {
            KRATOS_ERROR_IF(&(*existing_condition_it) != (pNewCondition.get()))//check if the pointee coincides
                << "attempting to add pNewCondition with Id :" << pNewCondition->Id() << ", unfortunately a (different) condition with the same Id already exists" << std::endl;
        }
    }
}

/** Inserts a list of conditions to a submodelpart provided their Id. Does nothing if applied to the top model part
*/
void ModelPart::AddConditions(std::vector<IndexType> const& ConditionIds, IndexType ThisIndex)
{
    ModelPartHelperUtilities::AddEntitiesFromIds([](ModelPart* pModelPart) { return &pModelPart->Conditions(); }, this, ConditionIds);
}

/** Inserts a condition in the mesh with ThisIndex.
*/
ModelPart::ConditionType::Pointer ModelPart::CreateNewCondition(std::string ConditionName,
        ModelPart::IndexType Id, std::vector<IndexType> ConditionNodeIds,
        ModelPart::PropertiesType::Pointer pProperties, ModelPart::IndexType ThisIndex)
{
    KRATOS_TRY

    if (IsSubModelPart()) {
        ConditionType::Pointer p_new_condition = mpParentModelPart->CreateNewCondition(ConditionName, Id, ConditionNodeIds, pProperties, ThisIndex);
        GetMesh(ThisIndex).AddCondition(p_new_condition);
        return p_new_condition;
    }

    Geometry< Node >::PointsArrayType pConditionNodes;

    for (unsigned int i = 0; i < ConditionNodeIds.size(); i++)
    {
        pConditionNodes.push_back(pGetNode(ConditionNodeIds[i]));
    }

    return CreateNewCondition(ConditionName, Id, pConditionNodes, pProperties, ThisIndex);
    KRATOS_CATCH("")
}

/** Inserts a condition in the mesh with ThisIndex.
*/
ModelPart::ConditionType::Pointer ModelPart::CreateNewCondition(std::string ConditionName,
        ModelPart::IndexType Id, Geometry< Node >::PointsArrayType pConditionNodes,
        ModelPart::PropertiesType::Pointer pProperties, ModelPart::IndexType ThisIndex)
{
    KRATOS_TRY
    if (IsSubModelPart())
    {
        ConditionType::Pointer p_new_condition = mpParentModelPart->CreateNewCondition(ConditionName, Id, pConditionNodes, pProperties, ThisIndex);
        GetMesh(ThisIndex).AddCondition(p_new_condition);
        return p_new_condition;
    }

    auto existing_condition_iterator = GetMesh(ThisIndex).Conditions().find(Id);
    KRATOS_ERROR_IF(existing_condition_iterator != GetMesh(ThisIndex).ConditionsEnd() )
        << "trying to construct a condition with ID " << Id << " however a condition with the same Id already exists";

    //get the condition
    ConditionType const& r_clone_condition = KratosComponents<ConditionType>::Get(ConditionName);
    ConditionType::Pointer p_condition = r_clone_condition.Create(Id, pConditionNodes, pProperties);

    //add the new element
    GetMesh(ThisIndex).AddCondition(p_condition);

    return p_condition;
    KRATOS_CATCH("")
}

/** Inserts a condition in the mesh with ThisIndex.
*/
ModelPart::ConditionType::Pointer ModelPart::CreateNewCondition(std::string ConditionName,
        ModelPart::IndexType Id, typename GeometryType::Pointer pGeometry,
        ModelPart::PropertiesType::Pointer pProperties, ModelPart::IndexType ThisIndex)
{
    KRATOS_TRY
    if (IsSubModelPart())
    {
        ConditionType::Pointer p_new_condition = mpParentModelPart->CreateNewCondition(ConditionName, Id, pGeometry, pProperties, ThisIndex);
        GetMesh(ThisIndex).AddCondition(p_new_condition);
        return p_new_condition;
    }

    auto existing_condition_iterator = GetMesh(ThisIndex).Conditions().find(Id);
    KRATOS_ERROR_IF(existing_condition_iterator != GetMesh(ThisIndex).ConditionsEnd() )
        << "trying to construct a condition with ID " << Id << " however a condition with the same Id already exists";

    //get the condition
    ConditionType const& r_clone_condition = KratosComponents<ConditionType>::Get(ConditionName);
    ConditionType::Pointer p_condition = r_clone_condition.Create(Id, pGeometry, pProperties);

    //add the new element
    GetMesh(ThisIndex).AddCondition(p_condition);

    return p_condition;
    KRATOS_CATCH("")
}

/** Remove the condition with given Id from mesh with ThisIndex in this modelpart and all its subs.
*/
void ModelPart::RemoveCondition(ModelPart::IndexType ConditionId, ModelPart::IndexType ThisIndex)
{
    GetMesh(ThisIndex).RemoveCondition(ConditionId);

    for (SubModelPartIterator i_sub_model_part = SubModelPartsBegin(); i_sub_model_part != SubModelPartsEnd(); i_sub_model_part++)
        i_sub_model_part->RemoveCondition(ConditionId, ThisIndex);
}

/** Remove given condition from mesh with ThisIndex in this modelpart and all its subs.
*/
void ModelPart::RemoveCondition(ModelPart::ConditionType& ThisCondition, ModelPart::IndexType ThisIndex)
{
    GetMesh(ThisIndex).RemoveCondition(ThisCondition);

    for (SubModelPartIterator i_sub_model_part = SubModelPartsBegin(); i_sub_model_part != SubModelPartsEnd(); i_sub_model_part++)
        i_sub_model_part->RemoveCondition(ThisCondition, ThisIndex);
}

/** Remove given condition from mesh with ThisIndex in this modelpart and all its subs.
*/
void ModelPart::RemoveCondition(ModelPart::ConditionType::Pointer pThisCondition, ModelPart::IndexType ThisIndex)
{
    GetMesh(ThisIndex).RemoveCondition(pThisCondition);

    for (SubModelPartIterator i_sub_model_part = SubModelPartsBegin(); i_sub_model_part != SubModelPartsEnd(); i_sub_model_part++)
        i_sub_model_part->RemoveCondition(pThisCondition, ThisIndex);
}

/** Remove the condition with given Id from mesh with ThisIndex in parents, itself and children.
*/
void ModelPart::RemoveConditionFromAllLevels(ModelPart::IndexType ConditionId, ModelPart::IndexType ThisIndex)
{
    if (IsSubModelPart())
    {
        mpParentModelPart->RemoveCondition(ConditionId, ThisIndex);
        return;
    }

    RemoveCondition(ConditionId, ThisIndex);
}

/** Remove given condition from mesh with ThisIndex in parents, itself and children.
*/
void ModelPart::RemoveConditionFromAllLevels(ModelPart::ConditionType& ThisCondition, ModelPart::IndexType ThisIndex)
{
    if (IsSubModelPart())
    {
        mpParentModelPart->RemoveCondition(ThisCondition, ThisIndex);
        return;
    }

    RemoveCondition(ThisCondition, ThisIndex);
}

/** Remove given condition from mesh with ThisIndex in parents, itself and children.
*/
void ModelPart::RemoveConditionFromAllLevels(ModelPart::ConditionType::Pointer pThisCondition, ModelPart::IndexType ThisIndex)
{
    if (IsSubModelPart())
    {
        mpParentModelPart->RemoveCondition(pThisCondition, ThisIndex);
        return;
    }

    RemoveCondition(pThisCondition, ThisIndex);
}

void ModelPart::RemoveConditions(Flags IdentifierFlag)
{
    // This method is optimized to free the memory
    //loop over all the meshes
    for(auto& r_mesh : this->GetMeshes())  {
        ModelPartHelperUtilities::RemoveEntities<ModelPart::ConditionsContainerType>(r_mesh, IdentifierFlag);
    }

    //now recursively remove the conditions in the submodelparts
    for (auto& r_sub_model_part : this->SubModelParts()) {
        r_sub_model_part.RemoveConditions(IdentifierFlag);
    }
}

void ModelPart::RemoveConditionsFromAllLevels(Flags IdentifierFlag)
{
    ModelPart& root_model_part = GetRootModelPart();
    root_model_part.RemoveConditions(IdentifierFlag);
}

///@}
///@name Geometry Container
///@{

ModelPart::GeometryType::Pointer ModelPart::CreateNewGeometry(
    const std::string& rGeometryTypeName,
    const std::vector<IndexType>& rGeometryNodeIds
    )
{
    if (IsSubModelPart()) {
        GeometryType::Pointer p_new_geometry = mpParentModelPart->CreateNewGeometry(rGeometryTypeName, rGeometryNodeIds);
        this->AddGeometry(p_new_geometry);
        return p_new_geometry;
    }

    GeometryType::PointsArrayType p_geometry_nodes;
    for (IndexType i = 0; i < rGeometryNodeIds.size(); ++i) {
        p_geometry_nodes.push_back(pGetNode(rGeometryNodeIds[i]));
    }

    return CreateNewGeometry(rGeometryTypeName, p_geometry_nodes);
}

ModelPart::GeometryType::Pointer ModelPart::CreateNewGeometry(
    const std::string& rGeometryTypeName,
    GeometryType::PointsArrayType pGeometryNodes
    )
{
    KRATOS_TRY

        if (IsSubModelPart()) {
            GeometryType::Pointer p_new_geometry = mpParentModelPart->CreateNewGeometry(rGeometryTypeName, pGeometryNodes);
            this->AddGeometry(p_new_geometry);
            return p_new_geometry;
        }

    // Create the new geometry
    GeometryType const& r_clone_geometry = KratosComponents<GeometryType>::Get(rGeometryTypeName);
    GeometryType::Pointer p_geometry = r_clone_geometry.Create(pGeometryNodes);

    //add the new geometry
    this->AddGeometry(p_geometry);

    return p_geometry;

    KRATOS_CATCH("")
}

ModelPart::GeometryType::Pointer ModelPart::CreateNewGeometry(
    const std::string& rGeometryTypeName,
    GeometryType::Pointer pGeometry
    )
{
    KRATOS_TRY

        if (IsSubModelPart()) {
            GeometryType::Pointer p_new_geometry = mpParentModelPart->CreateNewGeometry(rGeometryTypeName, pGeometry);
            this->AddGeometry(p_new_geometry);
            return p_new_geometry;
        }

    // Create the new geometry
    GeometryType const& r_clone_geometry = KratosComponents<GeometryType>::Get(rGeometryTypeName);
    GeometryType::Pointer p_geometry = r_clone_geometry.Create(*pGeometry);

    //add the new geometry
    this->AddGeometry(p_geometry);

    return p_geometry;

    KRATOS_CATCH("")
}

ModelPart::GeometryType::Pointer ModelPart::CreateNewGeometry(
    const std::string& rGeometryTypeName,
    const IndexType GeometryId,
    const std::vector<IndexType>& rGeometryNodeIds
    )
{
    if (IsSubModelPart()) {
        GeometryType::Pointer p_new_geometry = mpParentModelPart->CreateNewGeometry(rGeometryTypeName, GeometryId, rGeometryNodeIds);
        this->AddGeometry(p_new_geometry);
        return p_new_geometry;
    }

    GeometryType::PointsArrayType p_geometry_nodes;
    for (IndexType i = 0; i < rGeometryNodeIds.size(); ++i) {
        p_geometry_nodes.push_back(pGetNode(rGeometryNodeIds[i]));
    }

    return CreateNewGeometry(rGeometryTypeName, GeometryId, p_geometry_nodes);
}

ModelPart::GeometryType::Pointer ModelPart::CreateNewGeometry(
    const std::string& rGeometryTypeName,
    const IndexType GeometryId,
    GeometryType::PointsArrayType pGeometryNodes
    )
{
    KRATOS_TRY

    if (IsSubModelPart()) {
        GeometryType::Pointer p_new_geometry = mpParentModelPart->CreateNewGeometry(rGeometryTypeName, GeometryId, pGeometryNodes);
        this->AddGeometry(p_new_geometry);
        return p_new_geometry;
    }

    // Check if the geometry already exists
    if (this->HasGeometry(GeometryId)) {
        // Get the existing geometry with the same Id
        const auto p_existing_geom = this->pGetGeometry(GeometryId);

        // Check if the existing geometry has the same type
        KRATOS_ERROR_IF_NOT(GeometryType::HasSameGeometryType(*p_existing_geom, KratosComponents<GeometryType>::Get(rGeometryTypeName)))
            << "Attempting to add geometry with Id: " << GeometryId << ". A different geometry with the same Id already exists." << std::endl;

        // Check if the connectivities (nodes) of the existing geometry match that of the new
        for (IndexType i = 0; i < p_existing_geom->PointsNumber(); ++i) {
            KRATOS_ERROR_IF_NOT((p_existing_geom->operator()(i)).get() == &(pGeometryNodes[i]))
                << "Attempting to add a new geometry with Id: " << GeometryId << ". A same type geometry with same Id but different connectivities already exists." << std::endl;
        }

        // Return the existing geometry
        return p_existing_geom;
    }

    // Create the new geometry
    GeometryType const& r_clone_geometry = KratosComponents<GeometryType>::Get(rGeometryTypeName);
    GeometryType::Pointer p_geometry = r_clone_geometry.Create(GeometryId, pGeometryNodes);

    // Add the new geometry
    this->AddGeometry(p_geometry);

    // Return the new geometry
    return p_geometry;

    KRATOS_CATCH("")
}

ModelPart::GeometryType::Pointer ModelPart::CreateNewGeometry(
    const std::string& rGeometryTypeName,
    const IndexType GeometryId,
    GeometryType::Pointer pGeometry
    )
{
    KRATOS_TRY

    if (IsSubModelPart()) {
        GeometryType::Pointer p_new_geometry = mpParentModelPart->CreateNewGeometry(rGeometryTypeName, GeometryId, pGeometry);
        this->AddGeometry(p_new_geometry);
        return p_new_geometry;
    }

    // Check if the geometry already exists
    if (this->HasGeometry(GeometryId)) {
        // Get the existing geometry with the same Id
        const auto p_existing_geom = this->pGetGeometry(GeometryId);

        // Check if the existing geometry has the same type
        KRATOS_ERROR_IF_NOT(GeometryType::HasSameGeometryType(*p_existing_geom, KratosComponents<GeometryType>::Get(rGeometryTypeName)))
            << "Attempting to add geometry with Id: " << GeometryId << ". A different geometry with the same Id already exists." << std::endl;

        // Check if the connectivities (nodes) of the existing geometry match that of the new
        for (IndexType i = 0; i < p_existing_geom->PointsNumber(); ++i) {
            KRATOS_ERROR_IF_NOT((p_existing_geom->operator()(i)).get() == ((*pGeometry)(i)).get())
                << "Attempting to add a new geometry with Id: " << GeometryId << ". A same type geometry with same Id but different connectivities already exists." << std::endl;
        }

        // Return the existing geometry
        return p_existing_geom;
    }

    // Create the new geometry
    GeometryType const& r_clone_geometry = KratosComponents<GeometryType>::Get(rGeometryTypeName);
    GeometryType::Pointer p_geometry = r_clone_geometry.Create(GeometryId, *pGeometry);

    // Add the new geometry
    this->AddGeometry(p_geometry);

    // Return the new geometry
    return p_geometry;

    KRATOS_CATCH("")
}

ModelPart::GeometryType::Pointer ModelPart::CreateNewGeometry(
    const std::string& rGeometryTypeName,
    const std::string& rGeometryIdentifierName,
    const std::vector<IndexType>& rGeometryNodeIds
    )
{
    if (IsSubModelPart()) {
        GeometryType::Pointer p_new_geometry = mpParentModelPart->CreateNewGeometry(rGeometryTypeName, rGeometryIdentifierName, rGeometryNodeIds);
        this->AddGeometry(p_new_geometry);
        return p_new_geometry;
    }

    GeometryType::PointsArrayType p_geometry_nodes;
    for (IndexType i = 0; i < rGeometryNodeIds.size(); ++i) {
        p_geometry_nodes.push_back(pGetNode(rGeometryNodeIds[i]));
    }

    return CreateNewGeometry(rGeometryTypeName, rGeometryIdentifierName, p_geometry_nodes);
}

ModelPart::GeometryType::Pointer ModelPart::CreateNewGeometry(
    const std::string& rGeometryTypeName,
    const std::string& rGeometryIdentifierName,
    GeometryType::PointsArrayType pGeometryNodes
    )
{
    KRATOS_TRY

    if (IsSubModelPart()) {
        GeometryType::Pointer p_new_geometry = mpParentModelPart->CreateNewGeometry(rGeometryTypeName, rGeometryIdentifierName, pGeometryNodes);
        this->AddGeometry(p_new_geometry);
        return p_new_geometry;
    }

    // Check if the geometry already exists
    if (this->HasGeometry(rGeometryIdentifierName)) {
        // Get the existing geometry with the same Id
        const auto p_existing_geom = this->pGetGeometry(rGeometryIdentifierName);

        // Check if the existing geometry has the same type
        KRATOS_ERROR_IF_NOT(GeometryType::HasSameGeometryType(*p_existing_geom, KratosComponents<GeometryType>::Get(rGeometryTypeName)))
            << "Attempting to add geometry with Id: " << rGeometryIdentifierName << ". A different geometry with the same Id already exists." << std::endl;

        // Check if the connectivities (nodes) of the existing geometry match that of the new
        for (IndexType i = 0; i < p_existing_geom->PointsNumber(); ++i) {
            KRATOS_ERROR_IF_NOT((p_existing_geom->operator()(i)).get() == &(pGeometryNodes[i]))
                << "Attempting to add a new geometry with Id: " << rGeometryIdentifierName << ". A same type geometry with same Id but different connectivities already exists." << std::endl;
        }

        // Return the existing geometry
        return p_existing_geom;
    }

    // Create the new geometry
    GeometryType const& r_clone_geometry = KratosComponents<GeometryType>::Get(rGeometryTypeName);
    GeometryType::Pointer p_geometry = r_clone_geometry.Create(rGeometryIdentifierName, pGeometryNodes);

    //add the new geometry
    this->AddGeometry(p_geometry);

    return p_geometry;

    KRATOS_CATCH("")
}

ModelPart::GeometryType::Pointer ModelPart::CreateNewGeometry(
    const std::string& rGeometryTypeName,
    const std::string& rGeometryIdentifierName,
    GeometryType::Pointer pGeometry
    )
{
    KRATOS_TRY

    if (IsSubModelPart()) {
        GeometryType::Pointer p_new_geometry = mpParentModelPart->CreateNewGeometry(rGeometryTypeName, rGeometryIdentifierName, pGeometry);
        this->AddGeometry(p_new_geometry);
        return p_new_geometry;
    }

    // Check if the geometry already exists
    if (this->HasGeometry(rGeometryIdentifierName)) {
        // Get the existing geometry with the same Id
        const auto p_existing_geom = this->pGetGeometry(rGeometryIdentifierName);

        // Check if the existing geometry has the same type
        KRATOS_ERROR_IF_NOT(GeometryType::HasSameGeometryType(*p_existing_geom, KratosComponents<GeometryType>::Get(rGeometryTypeName)))
            << "Attempting to add geometry with Id: " << rGeometryIdentifierName << ". A different geometry with the same Id already exists." << std::endl;

        // Check if the connectivities (nodes) of the existing geometry match that of the new
        for (IndexType i = 0; i < p_existing_geom->PointsNumber(); ++i) {
            KRATOS_ERROR_IF_NOT((p_existing_geom->operator()(i)).get() == ((*pGeometry)(i)).get())
                << "Attempting to add a new geometry with Id: " << rGeometryIdentifierName << ". A same type geometry with same Id but different connectivities already exists." << std::endl;
        }

        // Return the existing geometry
        return p_existing_geom;
    }

    // Create the new geometry
    GeometryType const& r_clone_geometry = KratosComponents<GeometryType>::Get(rGeometryTypeName);
    GeometryType::Pointer p_geometry = r_clone_geometry.Create(rGeometryIdentifierName, *pGeometry);

    //add the new geometry
    this->AddGeometry(p_geometry);

    return p_geometry;

    KRATOS_CATCH("")
}

/// Adds a geometry to the geometry container.
void ModelPart::AddGeometry(
    typename GeometryType::Pointer pNewGeometry)
{
    if (IsSubModelPart()) {
        if (!mpParentModelPart->HasGeometry(pNewGeometry->Id())) {
            mpParentModelPart->AddGeometry(pNewGeometry);
        }
    }
    /// Check if geometry id already used, is done within the geometry container.
    mGeometries.AddGeometry(pNewGeometry);
}

/** Inserts a list of geometries to a submodelpart provided their Id. Does nothing if applied to the top model part
 */
void ModelPart::AddGeometries(std::vector<IndexType> const& GeometriesIds)
{
    ModelPartHelperUtilities::AddEntitiesFromIds([](ModelPart* pModelPart) { return &pModelPart->Geometries(); }, this, GeometriesIds);
}

/// Removes a geometry by id.
void ModelPart::RemoveGeometry(
    const IndexType GeometryId)
{
    mGeometries.RemoveGeometry(GeometryId);

    for (SubModelPartIterator i_sub_model_part = SubModelPartsBegin();
        i_sub_model_part != SubModelPartsEnd();
        ++i_sub_model_part)
        i_sub_model_part->RemoveGeometry(GeometryId);
}

/// Removes a geometry by name.
void ModelPart::RemoveGeometry(
    std::string GeometryName)
{
    mGeometries.RemoveGeometry(GeometryName);

    for (SubModelPartIterator i_sub_model_part = SubModelPartsBegin();
        i_sub_model_part != SubModelPartsEnd();
        ++i_sub_model_part)
        i_sub_model_part->RemoveGeometry(GeometryName);
}

/// Removes a geometry by id from all root and sub model parts.
void ModelPart::RemoveGeometryFromAllLevels(const IndexType GeometryId)
{
    if (IsSubModelPart())
    {
        mpParentModelPart->RemoveGeometry(GeometryId);
        return;
    }

    RemoveGeometry(GeometryId);
}

/// Removes a geometry by name from all root and sub model parts.
void ModelPart::RemoveGeometryFromAllLevels(std::string GeometryName)
{
    if (IsSubModelPart())
    {
        mpParentModelPart->RemoveGeometry(GeometryName);
        return;
    }

    RemoveGeometry(GeometryName);
}

///@}
///@name Sub Model Parts
///@{

ModelPart& ModelPart::CreateSubModelPart(std::string const& NewSubModelPartName)
{
    const auto delim_pos = NewSubModelPartName.find('.');
    const std::string& sub_model_part_name = NewSubModelPartName.substr(0, delim_pos);

    if (delim_pos == std::string::npos) {
        KRATOS_ERROR_IF(mSubModelParts.find(NewSubModelPartName) != mSubModelParts.end())
            << "There is an already existing sub model part with name \"" << NewSubModelPartName
            << "\" in model part: \"" << FullName() << "\"" << std::endl;

        ModelPart* praw = new ModelPart(NewSubModelPartName, this->mpVariablesList, this->GetModel());
        Kratos::shared_ptr<ModelPart> p_model_part(praw); //we need to construct first a raw pointer
        p_model_part->SetParentModelPart(this);
        p_model_part->mBufferSize = this->mBufferSize;
        p_model_part->mpProcessInfo = this->mpProcessInfo;
        mSubModelParts.insert(p_model_part);
        return *p_model_part;
    } else {
        ModelPart *p;
        SubModelPartIterator i = mSubModelParts.find(sub_model_part_name);
        if (i == mSubModelParts.end()) {
            p = &CreateSubModelPart(sub_model_part_name);
        } else {
            p = &(*i);
        }
        return p->CreateSubModelPart(NewSubModelPartName.substr(delim_pos + 1));
    }
}

ModelPart& ModelPart::GetSubModelPart(std::string const& SubModelPartName)
{
    const auto delim_pos = SubModelPartName.find('.');
    const std::string& sub_model_part_name = SubModelPartName.substr(0, delim_pos);

    SubModelPartIterator i = mSubModelParts.find(sub_model_part_name);
    if (i == mSubModelParts.end()) {
        ErrorNonExistingSubModelPart(sub_model_part_name);
    }

    if (delim_pos == std::string::npos) {
        return *i;
    } else {
        return i->GetSubModelPart(SubModelPartName.substr(delim_pos + 1));
    }
}

const ModelPart& ModelPart::GetSubModelPart(std::string const& SubModelPartName) const
{
    const auto delim_pos = SubModelPartName.find('.');
    const std::string& r_sub_model_part_name = SubModelPartName.substr(0, delim_pos);

    const auto i = mSubModelParts.find(r_sub_model_part_name);
    if (i == mSubModelParts.end()) {
        ErrorNonExistingSubModelPart(r_sub_model_part_name);
    }

    if (delim_pos == std::string::npos) {
        return *i;
    } else {
        return i->GetSubModelPart(SubModelPartName.substr(delim_pos + 1));
    }
}

ModelPart* ModelPart::pGetSubModelPart(std::string const& SubModelPartName)
{
    const auto delim_pos = SubModelPartName.find('.');
    const std::string& sub_model_part_name = SubModelPartName.substr(0, delim_pos);

    SubModelPartIterator i = mSubModelParts.find(sub_model_part_name);
    if (i == mSubModelParts.end()) {
        ErrorNonExistingSubModelPart(sub_model_part_name);
    }

    if (delim_pos == std::string::npos) {
        return  (i.base()->second).get();
    } else {
        return i->pGetSubModelPart(SubModelPartName.substr(delim_pos + 1));
    }
}

/** Remove a sub modelpart with given name.
*/
void ModelPart::RemoveSubModelPart(std::string const& ThisSubModelPartName)
{
    const auto delim_pos = ThisSubModelPartName.find('.');
    const std::string& sub_model_part_name = ThisSubModelPartName.substr(0, delim_pos);

    SubModelPartIterator i = mSubModelParts.find(sub_model_part_name);
    if (delim_pos == std::string::npos) {
        if (i == mSubModelParts.end()) {
            std::stringstream warning_msg;
            warning_msg << "Trying to remove sub model part with name \"" << ThisSubModelPartName
                    << "\" in model part \"" << FullName() << "\" which does not exist.\n"
                    << "The the following sub model parts are available:";
            for (const auto& r_avail_smp_name : GetSubModelPartNames()) {
                warning_msg << "\n\t" "\"" << r_avail_smp_name << "\"";
            }
            KRATOS_WARNING("ModelPart") << warning_msg.str() << std::endl;
        } else {
            mSubModelParts.erase(ThisSubModelPartName);
        }
    } else {
        if (i == mSubModelParts.end()) {
            ErrorNonExistingSubModelPart(sub_model_part_name);
        }

        return i->RemoveSubModelPart(ThisSubModelPartName.substr(delim_pos + 1));
    }
}

/** Remove given sub model part.
*/
void ModelPart::RemoveSubModelPart(ModelPart& ThisSubModelPart)
{
    std::string name = ThisSubModelPart.Name();
    // finding the sub model part
    SubModelPartIterator i_sub_model_part = mSubModelParts.find(name);

    KRATOS_ERROR_IF(i_sub_model_part == mSubModelParts.end()) << "The sub model part  \"" << name << "\" does not exist in the \"" << Name() << "\" model part to be removed" << std::endl;

    mSubModelParts.erase(name);
}


ModelPart& ModelPart::GetParentModelPart()
{
    if (IsSubModelPart()) {
        return *mpParentModelPart;
    } else {
        return *this;
    }
}

const ModelPart& ModelPart::GetParentModelPart() const
{
    if (IsSubModelPart()) {
        return *mpParentModelPart;
    } else {
        return *this;
    }
}

bool ModelPart::HasSubModelPart(std::string const& ThisSubModelPartName) const
{
    const auto delim_pos = ThisSubModelPartName.find('.');
    const std::string& sub_model_part_name = ThisSubModelPartName.substr(0, delim_pos);

    auto i = mSubModelParts.find(sub_model_part_name);
    if (i == mSubModelParts.end()) {
        return false;
    } else {
        if (delim_pos != std::string::npos) {
            return i->HasSubModelPart(ThisSubModelPartName.substr(delim_pos + 1));
        } else {
            return true;
        }
    }
}

std::vector<std::string> ModelPart::GetSubModelPartNames() const
{
    std::vector<std::string> SubModelPartsNames;
    SubModelPartsNames.reserve(NumberOfSubModelParts());

    for(auto& r_sub_model_part : mSubModelParts) {
        SubModelPartsNames.push_back(r_sub_model_part.Name());
    }

    return SubModelPartsNames;
}

void ModelPart::SetBufferSize(ModelPart::IndexType NewBufferSize)
{
    KRATOS_ERROR_IF(IsSubModelPart()) << "Calling the method of the sub model part "
        << Name() << " please call the one of the root model part: "
        << GetRootModelPart().Name() << std::endl;

    for(auto& r_sub_model_part : mSubModelParts) {
        r_sub_model_part.SetBufferSizeSubModelParts(NewBufferSize);
    }

    mBufferSize = NewBufferSize;

    auto& r_nodes = Nodes();
    IndexPartition<size_t>(r_nodes.size()).for_each([&](size_t i){
        auto node_iterator = r_nodes.begin() + i;
        node_iterator->SetBufferSize(mBufferSize);
    });

}

void ModelPart::SetBufferSizeSubModelParts(ModelPart::IndexType NewBufferSize)
{
    for(auto& r_sub_model_part : mSubModelParts) {
        r_sub_model_part.SetBufferSizeSubModelParts(NewBufferSize);
    }

    mBufferSize = NewBufferSize;
}

/// run input validation
int ModelPart::Check() const
{
    KRATOS_TRY

    const ProcessInfo& r_current_process_info = this->GetProcessInfo();

    // Checks for all of the elements
    block_for_each(this->Elements(), [&r_current_process_info](const Element& rElement){
        rElement.Check(r_current_process_info);
    });

    // Checks for all of the conditions
    block_for_each(this->Conditions(), [&r_current_process_info](const Condition& rCondition){
        rCondition.Check(r_current_process_info);
    });

    // Checks for all of the constraints
    block_for_each(this->MasterSlaveConstraints(), [&r_current_process_info](const MasterSlaveConstraint& rConstraint){
        rConstraint.Check(r_current_process_info);
    });

    return 0;


    KRATOS_CATCH("");
}

/// Turn back information as a string.
std::string ModelPart::Info() const
{
    return "-" + mName + "- model part";
}

/// Print information about this object.

void ModelPart::PrintInfo(std::ostream& rOStream) const
{
    rOStream << Info();
}

/// Print object's data.

void ModelPart::PrintData(std::ostream& rOStream) const
{
    DataValueContainer::PrintData(rOStream);

    if (!IsSubModelPart()) {
        rOStream  << "    Buffer Size : " << mBufferSize << std::endl;
    }
    rOStream << "    Number of tables : " << NumberOfTables() << std::endl;
    rOStream << "    Number of sub model parts : " << NumberOfSubModelParts() << std::endl;
    if (!IsSubModelPart()) {
        if (IsDistributed()) {
            rOStream << "    Distributed; Communicator has " << mpCommunicator->TotalProcesses() << " total processes" << std::endl;
        }
        mpProcessInfo->PrintData(rOStream);
    }
    rOStream << std::endl;
    rOStream << "    Number of Geometries  : " << mGeometries.NumberOfGeometries() << std::endl;
    for (IndexType i = 0; i < mMeshes.size(); i++) {
        rOStream << "    Mesh " << i << " :" << std::endl;
        GetMesh(i).PrintData(rOStream, "    ");
    }
    rOStream << std::endl;

    // Printing the submodelparts by their names in alphabetical order
    std::vector< std::string > submodel_part_names;
    submodel_part_names.reserve(NumberOfSubModelParts());
    for (const auto& r_sub_model_part : mSubModelParts) {
        submodel_part_names.push_back(r_sub_model_part.Name());
    }
    std::sort(submodel_part_names.begin(),submodel_part_names.end());

    for (const auto& r_sub_model_part_name : submodel_part_names) {
        const auto& r_sub_model_part = *(mSubModelParts.find(r_sub_model_part_name));
        r_sub_model_part.PrintInfo(rOStream, "    ");
        rOStream << std::endl;
        r_sub_model_part.PrintData(rOStream, "    ");
    }
}

/// Print information about this object.

void ModelPart::PrintInfo(std::ostream& rOStream, std::string const& PrefixString) const
{
    rOStream << PrefixString << Info();
}

/// Print object's data.

void ModelPart::PrintData(std::ostream& rOStream, std::string const& PrefixString) const
{
    if (!IsSubModelPart()) {
        rOStream << PrefixString << "    Buffer Size : " << mBufferSize << std::endl;
    }
    rOStream << PrefixString << "    Number of tables : " << NumberOfTables() << std::endl;
    rOStream << PrefixString << "    Number of sub model parts : " << NumberOfSubModelParts() << std::endl;

    if (!IsSubModelPart()) {
        mpProcessInfo->PrintData(rOStream);
    }
    rOStream << std::endl;
    rOStream << PrefixString << "    Number of Geometries  : " << mGeometries.NumberOfGeometries() << std::endl;

    for (IndexType i = 0; i < mMeshes.size(); i++) {
        rOStream << PrefixString << "    Mesh " << i << " :" << std::endl;
        GetMesh(i).PrintData(rOStream, PrefixString + "    ");
    }

    // Printing the submodelparts by their names in alphabetical order
    std::vector< std::string > submodel_part_names;
    submodel_part_names.reserve(NumberOfSubModelParts());
    for (const auto& r_sub_model_part : mSubModelParts) {
        submodel_part_names.push_back(r_sub_model_part.Name());
    }
    std::sort(submodel_part_names.begin(),submodel_part_names.end());

    for (const auto& r_sub_model_part_name : submodel_part_names) {
        const auto& r_sub_model_part = *(mSubModelParts.find(r_sub_model_part_name));
        r_sub_model_part.PrintInfo(rOStream, PrefixString + "    ");
        rOStream << std::endl;
        r_sub_model_part.PrintData(rOStream, PrefixString + "    ");
    }
}

void ModelPart::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DataValueContainer);
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Flags );
    rSerializer.save("Name", mName);
    rSerializer.save("Buffer Size", mBufferSize);
    rSerializer.save("ProcessInfo", mpProcessInfo);
    rSerializer.save("Tables", mTables);
    rSerializer.save("Variables List", mpVariablesList);
    rSerializer.save("Meshes", mMeshes);
    rSerializer.save("Geometries", mGeometries);

    rSerializer.save("NumberOfSubModelParts", NumberOfSubModelParts());

    for (SubModelPartConstantIterator i_sub_model_part = SubModelPartsBegin(); i_sub_model_part != SubModelPartsEnd(); i_sub_model_part++)
         rSerializer.save("SubModelPartName", i_sub_model_part->Name());

    for (SubModelPartConstantIterator i_sub_model_part = SubModelPartsBegin(); i_sub_model_part != SubModelPartsEnd(); i_sub_model_part++)
        rSerializer.save("SubModelPart", *(i_sub_model_part));
}

void ModelPart::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DataValueContainer);
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Flags );
    std::string ModelPartName;
    rSerializer.load("Name", ModelPartName);

    KRATOS_ERROR_IF(ModelPartName != mName) //checking if the name is correct
        << "trying to load a model part called :   " << ModelPartName << "    into an object named :   " << mName << " the two names should coincide but do not" << std::endl;

    rSerializer.load("Buffer Size", mBufferSize);
    rSerializer.load("ProcessInfo", mpProcessInfo);
    rSerializer.load("Tables", mTables);
    rSerializer.load("Variables List", mpVariablesList);
    rSerializer.load("Meshes", mMeshes);
    rSerializer.load("Geometries", mGeometries);

    SizeType number_of_submodelparts;
    rSerializer.load("NumberOfSubModelParts", number_of_submodelparts);

    std::vector< std::string > submodel_part_names;
    for(SizeType i=0; i<number_of_submodelparts; ++i)
    {
        std::string name;
        rSerializer.load("SubModelPartName",name);
        submodel_part_names.push_back(name);
    }

    for(const auto& name : submodel_part_names)
    {
        auto& subpart = CreateSubModelPart(name);
        rSerializer.load("SubModelPart",subpart);
    }

    for (SubModelPartIterator i_sub_model_part = SubModelPartsBegin(); i_sub_model_part != SubModelPartsEnd(); i_sub_model_part++)
        i_sub_model_part->SetParentModelPart(this);
}


void ModelPart::ErrorNonExistingSubModelPart(const std::string& rSubModelPartName) const
{
    std::stringstream err_msg;
    err_msg << "There is no sub model part with name \"" << rSubModelPartName
            << "\" in model part \"" << FullName() << "\"\n"
            << "The following sub model parts are available:";
    for (const auto& r_avail_smp_name : GetSubModelPartNames()) {
        err_msg << "\n\t" << "\""<<r_avail_smp_name << "\"";
    }
    KRATOS_ERROR << err_msg.str() << std::endl;
}

}  // namespace Kratos.
