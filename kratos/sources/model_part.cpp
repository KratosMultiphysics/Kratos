//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//
//

// System includes


// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/exception.h"

namespace Kratos
{
KRATOS_CREATE_LOCAL_FLAG(ModelPart, ALL_ENTITIES, 0);
KRATOS_CREATE_LOCAL_FLAG(ModelPart, OVERWRITE_ENTITIES, 1);

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

    KRATOS_ERROR_IF_NOT(NewName.find(".") == std::string::npos) << "Please don't use names containing (\".\") when creating a ModelPart (used in \"" << NewName << "\")" << std::endl;

    mName = NewName;
    MeshType mesh;
    mMeshes.push_back(Kratos::make_shared<MeshType>(mesh.Clone()));
    mpCommunicator->SetLocalMesh(pGetMesh());  // assigning the current mesh to the local mesh of communicator for openmp cases
}

/// Destructor.
ModelPart::~ModelPart()
{
    mpCommunicator->Clear();

    for(auto i_mesh = mMeshes.begin() ; i_mesh != mMeshes.end() ; i_mesh++)
      i_mesh->Clear();


//     if (!IsSubModelPart())
//       delete mpVariablesList;
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

    const int nnodes = static_cast<int>(Nodes().size());
    auto nodes_begin = NodesBegin();
    #pragma omp parallel for firstprivate(nodes_begin,nnodes)
    for(int i = 0; i<nnodes; ++i)
    {
        auto node_iterator = nodes_begin + i;
        node_iterator->CloneSolutionStepData();
    }

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
void ModelPart::AddNodes(std::vector<IndexType> const& NodeIds, IndexType ThisIndex)
{
    KRATOS_TRY
    if(IsSubModelPart()) //does nothing if we are on the top model part
    {
        //obtain from the root model part the corresponding list of nodes
        ModelPart* root_model_part = &this->GetRootModelPart();
        ModelPart::NodesContainerType  aux;
        aux.reserve(NodeIds.size());
        for(unsigned int i=0; i<NodeIds.size(); i++)
        {
            ModelPart::NodesContainerType::iterator it = root_model_part->Nodes().find(NodeIds[i]);
            if(it!=root_model_part->NodesEnd())
                aux.push_back(*(it.base()));
            else
                KRATOS_ERROR << "while adding nodes to submodelpart, the node with Id " << NodeIds[i] << " does not exist in the root model part";
        }

        ModelPart* current_part = this;
        while(current_part->IsSubModelPart())
        {
            for(auto it = aux.begin(); it!=aux.end(); it++)
                current_part->Nodes().push_back( *(it.base()) );

            current_part->Nodes().Unique();

            current_part = &(current_part->GetParentModelPart());
        }
    }

    KRATOS_CATCH("");
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
    // Lambda to remove nodes from a mesh
    auto remove_nodes_from_mesh = [&](ModelPart::MeshType& r_mesh) {
        //count the nodes to be erase
        const unsigned int nnodes = r_mesh.Nodes().size();
        unsigned int erase_count = 0;
        #pragma omp parallel for reduction(+:erase_count)
        for(int i=0; i<static_cast<int>(nnodes); ++i) {
            ModelPart::NodesContainerType::iterator i_node = r_mesh.NodesBegin() + i;

            if( i_node->IsNot(IdentifierFlag) )
                erase_count++;
        }

        ModelPart::NodesContainerType temp_nodes_container;
        temp_nodes_container.reserve(r_mesh.Nodes().size() - erase_count);

        temp_nodes_container.swap(r_mesh.Nodes());

        for(ModelPart::NodesContainerType::iterator i_node = temp_nodes_container.begin() ; i_node != temp_nodes_container.end() ; ++i_node) {
            if( i_node->IsNot(IdentifierFlag) )
                (r_mesh.Nodes()).push_back(std::move(*(i_node.base())));
        }
    };

    // This method is optimized to free the memory
    // Loop over all the local meshes (Is this still necessary with Submodelparts?)
    for(auto& r_mesh: this->GetMeshes()) {
        remove_nodes_from_mesh(r_mesh);
    }

    if (IsDistributed()) {
        // Mark the IdentifierFlag across partitions coherently
        this->GetCommunicator().SynchronizeOrNodalFlags(IdentifierFlag);

        // Remove the nodes from the mpi-interfaces in case there is any
        remove_nodes_from_mesh(this->GetCommunicator().LocalMesh());
        for(auto& r_mesh: this->GetCommunicator().LocalMeshes()) {
            remove_nodes_from_mesh(r_mesh);
        }

        remove_nodes_from_mesh(this->GetCommunicator().GhostMesh());
        for(auto& r_mesh: this->GetCommunicator().GhostMeshes()) {
            remove_nodes_from_mesh(r_mesh);
        }

        remove_nodes_from_mesh(this->GetCommunicator().InterfaceMesh());
        for(auto& r_mesh: this->GetCommunicator().InterfaceMeshes()) {
            remove_nodes_from_mesh(r_mesh);
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

    for (NodeIterator i_node = NodesBegin(); i_node != NodesEnd(); ++i_node)
        i_node->SetSolutionStepVariablesList(mpVariablesList);
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

ModelPart::PropertiesType::Pointer ModelPart::pGetProperties(
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

ModelPart::PropertiesType& ModelPart::GetProperties(
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
    KRATOS_TRY
    if(IsSubModelPart()) //does nothing if we are on the top model part
    {
        //obtain from the root model part the corresponding list of nodes
        ModelPart* root_model_part = &this->GetRootModelPart();
        ModelPart::ElementsContainerType  aux;
        aux.reserve(ElementIds.size());
        for(unsigned int i=0; i<ElementIds.size(); i++)
        {
            ModelPart::ElementsContainerType::iterator it = root_model_part->Elements().find(ElementIds[i]);
            if(it!=root_model_part->ElementsEnd())
                aux.push_back(*(it.base()));
            else
                KRATOS_ERROR << "the element with Id " << ElementIds[i] << " does not exist in the root model part";
        }

        ModelPart* current_part = this;
        while(current_part->IsSubModelPart())
        {
            for(auto it = aux.begin(); it!=aux.end(); it++)
                current_part->Elements().push_back( *(it.base()) );

            current_part->Elements().Unique();

            current_part = &(current_part->GetParentModelPart());
        }
    }
    KRATOS_CATCH("");
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

    Geometry< Node < 3 > >::PointsArrayType pElementNodes;

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
        ModelPart::IndexType Id, Geometry< Node < 3 > >::PointsArrayType pElementNodes,
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
    auto& meshes = this->GetMeshes();
    for(ModelPart::MeshesContainerType::iterator i_mesh = meshes.begin() ; i_mesh != meshes.end() ; i_mesh++)
    {
        //count the elements to be erase
        const unsigned int nelements = i_mesh->Elements().size();
        unsigned int erase_count = 0;
        #pragma omp parallel for reduction(+:erase_count)
        for(int i=0; i<static_cast<int>(nelements); ++i)
        {
            auto i_elem = i_mesh->ElementsBegin() + i;

            if( i_elem->IsNot(IdentifierFlag) )
                erase_count++;
        }

        ModelPart::ElementsContainerType temp_elements_container;
        temp_elements_container.reserve(i_mesh->Elements().size() - erase_count);

        temp_elements_container.swap(i_mesh->Elements());

        for(ModelPart::ElementsContainerType::iterator i_elem = temp_elements_container.begin() ; i_elem != temp_elements_container.end() ; i_elem++)
        {
            if( i_elem->IsNot(IdentifierFlag) )
                (i_mesh->Elements()).push_back(std::move(*(i_elem.base())));
        }
    }

    //now recursively remove the elements in the submodelparts
    for (SubModelPartIterator i_sub_model_part = SubModelPartsBegin(); i_sub_model_part != SubModelPartsEnd(); i_sub_model_part++)
        i_sub_model_part->RemoveElements(IdentifierFlag);
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
    KRATOS_TRY
    if(IsSubModelPart()) //does nothing if we are on the top model part
    {
        //obtain from the root model part the corresponding list of constraints
        ModelPart* root_model_part = &this->GetRootModelPart();
        ModelPart::MasterSlaveConstraintContainerType  aux;
        aux.reserve(MasterSlaveConstraintIds.size());
        for(unsigned int i=0; i<MasterSlaveConstraintIds.size(); i++)
        {
            ModelPart::MasterSlaveConstraintContainerType::iterator it = root_model_part->MasterSlaveConstraints().find(MasterSlaveConstraintIds[i]);
            if(it!=root_model_part->MasterSlaveConstraintsEnd())
                aux.push_back(*(it.base()));
            else
                KRATOS_ERROR << "the master-slave constraint with Id " << MasterSlaveConstraintIds[i] << " does not exist in the root model part";
        }

        ModelPart* current_part = this;
        while(current_part->IsSubModelPart())
        {
            for(auto it = aux.begin(); it!=aux.end(); it++)
                current_part->MasterSlaveConstraints().push_back( *(it.base()) );

            current_part->MasterSlaveConstraints().Unique();

            current_part = &(current_part->GetParentModelPart());
        }
    }
    KRATOS_CATCH("");
}

/** Inserts an master-slave constraint in the current mesh.
 */
ModelPart::MasterSlaveConstraintType::Pointer ModelPart::CreateNewMasterSlaveConstraint(const std::string& ConstraintName,
                                                                                    IndexType Id,
                                                                                    ModelPart::DofsVectorType& rMasterDofsVector,
                                                                                    ModelPart::DofsVectorType& rSlaveDofsVector,
                                                                                    const ModelPart::MatrixType& RelationMatrix,
                                                                                    const ModelPart::VectorType& ConstantVector,
                                                                                    IndexType ThisIndex)
{

    KRATOS_TRY
    if (IsSubModelPart())
    {
        ModelPart::MasterSlaveConstraintType::Pointer p_new_constraint = mpParentModelPart->CreateNewMasterSlaveConstraint(ConstraintName, Id, rMasterDofsVector,
                                                                                                                    rSlaveDofsVector,
                                                                                                                    RelationMatrix,
                                                                                                                    ConstantVector,
                                                                                                                    ThisIndex);
        GetMesh(ThisIndex).AddMasterSlaveConstraint(p_new_constraint);
        GetMesh(ThisIndex).MasterSlaveConstraints().Unique();

        return p_new_constraint;
    }

    auto existing_constraint_iterator = GetMesh(ThisIndex).MasterSlaveConstraints().find(Id);
    KRATOS_ERROR_IF(existing_constraint_iterator != GetMesh(ThisIndex).MasterSlaveConstraintsEnd() )
        << "trying to construct an master-slave constraint with ID " << Id << " however a constraint with the same Id already exists";


    //create the new element
    ModelPart::MasterSlaveConstraintType const& r_clone_constraint = KratosComponents<MasterSlaveConstraintType>::Get(ConstraintName);
    ModelPart::MasterSlaveConstraintType::Pointer p_new_constraint = r_clone_constraint.Create(Id, rMasterDofsVector,
                                                                                        rSlaveDofsVector,
                                                                                        RelationMatrix,
                                                                                        ConstantVector);

    GetMesh(ThisIndex).AddMasterSlaveConstraint(p_new_constraint);
    GetMesh(ThisIndex).MasterSlaveConstraints().Unique();

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

    KRATOS_TRY
    if (rMasterNode.HasDofFor(rMasterVariable) && rSlaveNode.HasDofFor(rSlaveVariable) )
    {
        if (IsSubModelPart())
        {
                ModelPart::MasterSlaveConstraintType::Pointer p_new_constraint = mpParentModelPart->CreateNewMasterSlaveConstraint(ConstraintName, Id, rMasterNode,
                                                                                                                            rMasterVariable,
                                                                                                                            rSlaveNode,
                                                                                                                            rSlaveVariable,
                                                                                                                            Weight,
                                                                                                                            Constant,
                                                                                                                            ThisIndex);

                GetMesh(ThisIndex).AddMasterSlaveConstraint(p_new_constraint);
                GetMesh(ThisIndex).MasterSlaveConstraints().Unique();
                return p_new_constraint;
        }

        KRATOS_ERROR_IF(GetMesh(ThisIndex).HasMasterSlaveConstraint(Id))
            << "trying to construct an master-slave constraint with ID " << Id << " however a constraint with the same Id already exists";


            //create the new element
        ModelPart::MasterSlaveConstraintType const& r_clone_constraint = KratosComponents<MasterSlaveConstraintType>::Get(ConstraintName);
        ModelPart::MasterSlaveConstraintType::Pointer p_new_constraint = r_clone_constraint.Create(Id, rMasterNode,
                                                                                                rMasterVariable,
                                                                                                rSlaveNode,
                                                                                                rSlaveVariable,
                                                                                                Weight,
                                                                                                Constant);

        GetMesh(ThisIndex).AddMasterSlaveConstraint(p_new_constraint);
        GetMesh(ThisIndex).MasterSlaveConstraints().Unique();
        return p_new_constraint;
    } else
    {
        KRATOS_ERROR << "Master or Slave node does not have requested DOF " <<std::endl;
    }

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
    auto& meshes = this->GetMeshes();
    for(auto it_mesh = meshes.begin() ; it_mesh != meshes.end() ; it_mesh++) {
        // Count the constraints to be erase
        const SizeType nconstraints = it_mesh->MasterSlaveConstraints().size();
        SizeType erase_count = 0;
        #pragma omp parallel for reduction(+:erase_count)
        for(int i=0; i<static_cast<int>(nconstraints); ++i) {
            auto it_const = it_mesh->MasterSlaveConstraintsBegin() + i;

            if( it_const->IsNot(IdentifierFlag) )
                erase_count++;
        }

        ModelPart::MasterSlaveConstraintContainerType temp_constraints_container;
        temp_constraints_container.reserve(it_mesh->MasterSlaveConstraints().size() - erase_count);

        temp_constraints_container.swap(it_mesh->MasterSlaveConstraints());

        for(auto it_const = temp_constraints_container.begin() ; it_const != temp_constraints_container.end(); it_const++) {
            if( it_const->IsNot(IdentifierFlag) )
                (it_mesh->MasterSlaveConstraints()).push_back(std::move(*(it_const.base())));
        }
    }

    // Now recursively remove the constraints in the submodelparts
    for (SubModelPartIterator i_sub_model_part = SubModelPartsBegin(); i_sub_model_part != SubModelPartsEnd(); i_sub_model_part++)
        i_sub_model_part->RemoveMasterSlaveConstraints(IdentifierFlag);
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
    KRATOS_TRY
    if(IsSubModelPart()) //does nothing if we are on the top model part
    {
        //obtain from the root model part the corresponding list of nodes
        ModelPart* root_model_part = &this->GetRootModelPart();
        ModelPart::ConditionsContainerType  aux;
        aux.reserve(ConditionIds.size());
        for(unsigned int i=0; i<ConditionIds.size(); i++)
        {
            ModelPart::ConditionsContainerType::iterator it = root_model_part->Conditions().find(ConditionIds[i]);
            if(it!=root_model_part->ConditionsEnd())
                aux.push_back(*(it.base()));
            else
                KRATOS_ERROR << "the condition with Id " << ConditionIds[i] << " does not exist in the root model part";
        }

        ModelPart* current_part = this;
        while(current_part->IsSubModelPart())
        {
            for(auto it = aux.begin(); it!=aux.end(); it++)
                current_part->Conditions().push_back( *(it.base()) );

            current_part->Conditions().Unique();

            current_part = &(current_part->GetParentModelPart());
        }
    }
    KRATOS_CATCH("");
}

/** Inserts a condition in the mesh with ThisIndex.
*/
ModelPart::ConditionType::Pointer ModelPart::CreateNewCondition(std::string ConditionName,
        ModelPart::IndexType Id, std::vector<IndexType> ConditionNodeIds,
        ModelPart::PropertiesType::Pointer pProperties, ModelPart::IndexType ThisIndex)
{
    KRATOS_TRY
    Geometry< Node < 3 > >::PointsArrayType pConditionNodes;

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
        ModelPart::IndexType Id, Geometry< Node < 3 > >::PointsArrayType pConditionNodes,
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
    auto& meshes = this->GetMeshes();
    for(ModelPart::MeshesContainerType::iterator i_mesh = meshes.begin() ; i_mesh != meshes.end() ; i_mesh++)
    {
        //count the conditions to be erase
        const unsigned int nconditions = i_mesh->Conditions().size();
        unsigned int erase_count = 0;
        #pragma omp parallel for reduction(+:erase_count)
        for(int i=0; i<static_cast<int>(nconditions); ++i)
        {
            auto i_cond = i_mesh->ConditionsBegin() + i;

            if( i_cond->IsNot(IdentifierFlag) )
                erase_count++;
        }

        ModelPart::ConditionsContainerType temp_conditions_container;
        temp_conditions_container.reserve(i_mesh->Conditions().size() - erase_count);

        temp_conditions_container.swap(i_mesh->Conditions());

        for(ModelPart::ConditionsContainerType::iterator i_cond = temp_conditions_container.begin() ; i_cond != temp_conditions_container.end() ; i_cond++)
        {
            if( i_cond->IsNot(IdentifierFlag) )
                (i_mesh->Conditions()).push_back(std::move(*(i_cond.base())));
        }
    }

    //now recursively remove the conditions in the submodelparts
    for (SubModelPartIterator i_sub_model_part = SubModelPartsBegin(); i_sub_model_part != SubModelPartsEnd(); i_sub_model_part++)
        i_sub_model_part->RemoveConditions(IdentifierFlag);
}

void ModelPart::RemoveConditionsFromAllLevels(Flags IdentifierFlag)
{
    ModelPart& root_model_part = GetRootModelPart();
    root_model_part.RemoveConditions(IdentifierFlag);
}

///@}
///@name Geometry Container
///@{

/// Adds a geometry to the geometry container.
void ModelPart::AddGeometry(
    typename GeometryType::Pointer pNewGeometry)
{
    if (IsSubModelPart())
    {
        mpParentModelPart->AddGeometry(pNewGeometry);
    }
    /// Check if geometry id already used, is done within the geometry container.
    mGeometries.AddGeometry(pNewGeometry);
}


/// Removes a geometry by id.
void ModelPart::RemoveGeometry(
    IndexType GeometryId)
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
void ModelPart::RemoveGeometryFromAllLevels(IndexType GeometryId)
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

ModelPart&  ModelPart::CreateSubModelPart(std::string const& NewSubModelPartName)
{
    // Here a warning would be enough. To be disscussed. Pooyan.
    KRATOS_ERROR_IF(mSubModelParts.find(NewSubModelPartName) != mSubModelParts.end())
        << "There is an already existing sub model part with name \"" << NewSubModelPartName
        << "\" in model part: \"" << Name() << "\"" << std::endl;

    ModelPart* praw = new ModelPart(NewSubModelPartName, this->mpVariablesList, this->GetModel());
    Kratos::shared_ptr<ModelPart>  p_model_part(praw); //we need to construct first a raw pointer
    p_model_part->SetParentModelPart(this);
    p_model_part->mBufferSize = this->mBufferSize;
    p_model_part->mpProcessInfo = this->mpProcessInfo;
    mSubModelParts.insert(p_model_part);
    return *p_model_part;
}

ModelPart& ModelPart::GetSubModelPart(std::string const& SubModelPartName)
{
    SubModelPartIterator i = mSubModelParts.find(SubModelPartName);
    KRATOS_ERROR_IF(i == mSubModelParts.end()) << "There is no sub model part with name: \"" << SubModelPartName << "\" in model part\"" << Name() << "\"" << std::endl;

    return *i;
}

ModelPart* ModelPart::pGetSubModelPart(std::string const& SubModelPartName)
{
    SubModelPartIterator i = mSubModelParts.find(SubModelPartName);
    KRATOS_ERROR_IF(i == mSubModelParts.end()) << "There is no sub model part with name: \"" << SubModelPartName << "\" in model part\"" << Name() << "\"" << std::endl;

    return (i.base()->second).get();
}

/** Remove a sub modelpart with given name.
*/
void ModelPart::RemoveSubModelPart(std::string const& ThisSubModelPartName)
{
    // finding the sub model part
    SubModelPartIterator i_sub_model_part = mSubModelParts.find(ThisSubModelPartName);

    if (i_sub_model_part == mSubModelParts.end())
        return; // TODO: send a warning here. Pooyan.

    // now erase the pointer from the list
    mSubModelParts.erase(ThisSubModelPartName);
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
    return (mSubModelParts.find(ThisSubModelPartName) != mSubModelParts.end());
}

std::vector<std::string> ModelPart::GetSubModelPartNames()
{
    std::vector<std::string> SubModelPartsNames;

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

    auto nodes_begin = NodesBegin();
    const int nnodes = static_cast<int>(Nodes().size());
    #pragma omp parallel for firstprivate(nodes_begin,nnodes)
    for(int i = 0; i<nnodes; ++i)
    {
        auto node_iterator = nodes_begin + i;
        node_iterator->SetBufferSize(mBufferSize);
    }

}

void ModelPart::SetBufferSizeSubModelParts(ModelPart::IndexType NewBufferSize)
{
    for(auto& r_sub_model_part : mSubModelParts) {
        r_sub_model_part.SetBufferSizeSubModelParts(NewBufferSize);
    }

    mBufferSize = NewBufferSize;
}

/// run input validation
int ModelPart::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY
    int err = 0;
    for (ElementConstantIterator elem_iterator = ElementsBegin(); elem_iterator != ElementsEnd(); elem_iterator++)
    {
        const auto& r_elem = *elem_iterator;
        err = r_elem.Check(rCurrentProcessInfo);
    }
    for (ConditionConstantIterator condition_iterator = ConditionsBegin(); condition_iterator != ConditionsEnd(); condition_iterator++)
    {
        const auto& r_cond = *condition_iterator;
        err = r_cond.Check(rCurrentProcessInfo);
    }
    for (MasterSlaveConstraintConstantIteratorType constraint_iterator = MasterSlaveConstraintsBegin();
            constraint_iterator != MasterSlaveConstraintsEnd(); constraint_iterator++)
    {
        const auto& r_constraint = *constraint_iterator;
        err = r_constraint.Check(rCurrentProcessInfo);
    }
    return err;
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

}  // namespace Kratos.
