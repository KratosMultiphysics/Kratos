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
#include "boost/make_shared.hpp"


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/exception.h"

namespace Kratos
{
KRATOS_CREATE_LOCAL_FLAG(ModelPart, ALL_ENTITIES, 0);
KRATOS_CREATE_LOCAL_FLAG(ModelPart, OVERWRITE_ENTITIES, 1);

/// Default constructor.
ModelPart::ModelPart()
    : DataValueContainer()
    , Flags()
    , mBufferSize(1)
    , mpProcessInfo(new ProcessInfo())
    , mIndices(1, 0)
    , mpVariablesList(new VariablesList)
    , mpCommunicator(new Communicator)
    , mpParentModelPart(NULL)
    , mSubModelParts()
{
    mName = "Default";
    MeshType mesh;
    mMeshes.push_back(Kratos::make_shared<MeshType>(mesh.Clone()));
    mpCommunicator->SetLocalMesh(pGetMesh());  // assigning the current mesh to the local mesh of communicator for openmp cases
}

/// Constructor with name
ModelPart::ModelPart(std::string const& NewName)
    : DataValueContainer()
    , Flags()
    , mBufferSize(1)
    , mpProcessInfo(new ProcessInfo())
    , mIndices(1, 0)
    , mpVariablesList(new VariablesList)
    , mpCommunicator(new Communicator)
    , mpParentModelPart(NULL)
    , mSubModelParts()
{
    KRATOS_ERROR_IF( NewName.empty() ) << "Please don't use empty names (\"\") when creating a ModelPart" << std::endl;
    mName = NewName;
    MeshType mesh;
    mMeshes.push_back(Kratos::make_shared<MeshType>(mesh.Clone()));
    mpCommunicator->SetLocalMesh(pGetMesh());  // assigning the current mesh to the local mesh of communicator for openmp cases
}

/// Constructor with name and bufferSize
ModelPart::ModelPart(std::string const& NewName, IndexType NewBufferSize)
    : DataValueContainer()
    , Flags()
    , mBufferSize(NewBufferSize)
    , mpProcessInfo(new ProcessInfo())
    , mIndices(NewBufferSize, 0)
    , mpVariablesList(new VariablesList)
    , mpCommunicator(new Communicator)
    , mpParentModelPart(NULL)
    , mSubModelParts()
{
    KRATOS_ERROR_IF( NewName.empty() ) << "Please don't use empty names (\"\") when creating a ModelPart" << std::endl;
    mName = NewName;
    MeshType mesh;
    mMeshes.push_back(Kratos::make_shared<MeshType>(mesh.Clone()));
    mpCommunicator->SetLocalMesh(pGetMesh());  // assigning the current mesh to the local mesh of communicator for openmp cases
}

/// Destructor.
ModelPart::~ModelPart()
{
    if (!IsSubModelPart()){

      for (NodeIterator i_node = NodesBegin(); i_node != NodesEnd(); i_node++)
	{
	  if (i_node->pGetVariablesList() == mpVariablesList)
	    i_node->ClearSolutionStepsData();
	}
    }

    mpCommunicator->Clear();

    for(auto i_mesh = mMeshes.begin() ; i_mesh != mMeshes.end() ; i_mesh++)
      i_mesh->Clear();


    if (!IsSubModelPart())
      delete mpVariablesList;
}

ModelPart::IndexType ModelPart::CreateSolutionStep()
{
    KRATOS_THROW_ERROR(std::logic_error, "This method needs updating and is not working. Pooyan", "")
    return 0;
}

ModelPart::IndexType ModelPart::CloneSolutionStep()
{
    if (IsSubModelPart())
        //Todo KRATOS_THROW_ERROR(std::logic_error, "Calling the method of the sub model part ", Name())
        KRATOS_ERROR << "Calling the CloneSolutionStep method of the sub model part " << Name()
                     << " please call the one of the parent modelpart : " << mpParentModelPart->Name() << std::endl;

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
    if (IsSubModelPart())
        KRATOS_THROW_ERROR(std::logic_error, "Calling the method of the sub model part ", Name())
        //KRATOS_ERROR << "Calling the CloneTimeStep method of the sub model part " << Name()
        //	<< " please call the one of the parent modelpart : " << mpParentModelPart->Name() << std::endl;

        IndexType new_index = CloneSolutionStep();
    mpProcessInfo->SetAsTimeStepInfo();

    return new_index;
}


ModelPart::IndexType ModelPart::CreateTimeStep(double NewTime)
{
    if (IsSubModelPart())
        KRATOS_THROW_ERROR(std::logic_error, "Calling the method of the sub model part ", Name())
        //KRATOS_ERROR << "Calling the CreateTimeStep method of the sub model part " << Name()
        //	<< " please call the one of the parent modelpart : " << mpParentModelPart->Name() << std::endl;

        IndexType new_index = CreateSolutionStep();
    mpProcessInfo->SetAsTimeStepInfo(NewTime);

    return new_index;
}

ModelPart::IndexType ModelPart::CloneTimeStep(double NewTime)
{
    if (IsSubModelPart())
        KRATOS_THROW_ERROR(std::logic_error, "Calling the CloneSolutionStep method of the sub model part ", Name())
        //	KRATOS_ERROR << "Calling the CloneTimeStep method of the sub model part " << Name()
        //	<< " please call the one of the parent modelpart : " << mpParentModelPart->Name() << std::endl;

        IndexType new_index = CloneSolutionStep();
    mpProcessInfo->SetAsTimeStepInfo(NewTime);

    return new_index;
}

void ModelPart::OverwriteSolutionStepData(IndexType SourceSolutionStepIndex, IndexType DestinationSourceSolutionStepIndex)
{
    if (IsSubModelPart())
        KRATOS_THROW_ERROR(std::logic_error, "Calling the method of the sub model part ", Name())
        //KRATOS_ERROR << "Calling the OverwriteSolutionStepData method of the sub model part " << Name()
        //	<< " please call the one of the parent modelpart : " << mpParentModelPart->Name() << std::endl;

        for (NodeIterator node_iterator = NodesBegin(); node_iterator != NodesEnd(); node_iterator++)
            node_iterator->OverwriteSolutionStepData(SourceSolutionStepIndex, DestinationSourceSolutionStepIndex);

}

void ModelPart::ReduceTimeStep(ModelPart& rModelPart, double NewTime)
{
    KRATOS_TRY

    //ATTENTION: this function does not touch the coordinates of the nodes.
    //It just resets the database values to the values at the beginning of the time step

    if (IsSubModelPart())
        KRATOS_THROW_ERROR(std::logic_error, "Calling the method of the sub model part ", Name())
        //	KRATOS_ERROR << "Calling the OverwriteSolutionStepData method of the sub model part " << Name()
        //				<< " please call the one of the parent modelpart : " << mpParentModelPart->Name() << std::endl;

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

            current_part = current_part->GetParentModelPart();
        }
    }

    KRATOS_CATCH("");
}



/** Inserts a node in the mesh with ThisIndex.
*/
ModelPart::NodeType::Pointer ModelPart::CreateNewNode(int Id, double x, double y, double z, VariablesList* pNewVariablesList, ModelPart::IndexType ThisIndex)
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
        double distance = sqrt( pow( existing_node_it->X() - x,2) + pow(existing_node_it->Y() - y,2) + pow(existing_node_it->Z() - z,2) );

        if(distance > std::numeric_limits<double>::epsilon()*1000)
            KRATOS_ERROR << "trying to create a node with Id " << Id << " however a node with the same Id already exists in the root model part. Existing node coordinates are " << existing_node_it->Coordinates() << " coordinates of the nodes we are attempting to create are :" << x << " " << y << " " << z;

        //if the node we attempt to create is in the same position as the one that is already there, we return the old one
        return *(existing_node_it.base());
    }

    //create a new node
    NodeType::Pointer p_new_node = Kratos::make_shared< NodeType >( Id, x, y, z );

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
        double distance = sqrt( pow( existing_node_it->X() - x,2) + pow(existing_node_it->Y() - y,2) + pow(existing_node_it->Z() - z,2) );

        if(distance > std::numeric_limits<double>::epsilon()*1000)
            KRATOS_ERROR << "trying to create a node with Id " << Id << " however a node with the same Id already exists in the root model part. Existing node coordinates are " << existing_node_it->Coordinates() << " coordinates of the nodes we are attempting to create are :" << x << " " << y << " " << z;

        //if the node we attempt to create is in the same position as the one that is already there, we return the old one
        return *(existing_node_it.base());
    }

    //create a new node
    NodeType::Pointer p_new_node = Kratos::make_shared< NodeType >( Id, x, y, z, mpVariablesList, pThisData, mBufferSize);
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

void ModelPart::RemoveNodes(Flags identifier_flag)
{
    // This method is optimized to free the memory
    //loop over all the meshes
    ModelPart::MeshesContainerType& meshes = this->GetMeshes();
    for(ModelPart::MeshesContainerType::iterator i_mesh = meshes.begin() ; i_mesh != meshes.end() ; i_mesh++)
    {
        //count the nodes to be erase
        const unsigned int nnodes = i_mesh->Nodes().size();
        unsigned int erase_count = 0;
        #pragma omp parallel for reduction(+:erase_count)
        for(int i=0; i<static_cast<int>(nnodes); ++i)
        {
            ModelPart::NodesContainerType::iterator i_node = i_mesh->NodesBegin() + i;

            if( i_node->IsNot(identifier_flag) )
                erase_count++;
        }

        ModelPart::NodesContainerType temp_nodes_container;
        temp_nodes_container.reserve(i_mesh->Nodes().size() - erase_count);

        temp_nodes_container.swap(i_mesh->Nodes());

        for(ModelPart::NodesContainerType::iterator i_node = temp_nodes_container.begin() ; i_node != temp_nodes_container.end() ; i_node++)
        {
            if( i_node->IsNot(identifier_flag) )
                (i_mesh->Nodes()).push_back(std::move(*(i_node.base())));
        }
    }

    //now recursively remove the nodes in the submodelparts
    for (SubModelPartIterator i_sub_model_part = SubModelPartsBegin(); i_sub_model_part != SubModelPartsEnd(); i_sub_model_part++)
        i_sub_model_part->RemoveNodes(identifier_flag);
}

void ModelPart::RemoveNodesFromAllLevels(Flags identifier_flag)
{
    ModelPart& root_model_part = GetRootModelPart();
    root_model_part.RemoveNodes(identifier_flag);
}

ModelPart& ModelPart::GetRootModelPart()
{
    if (IsSubModelPart())
        return mpParentModelPart->GetRootModelPart();
    else
        return *this;
}

void ModelPart::SetNodalSolutionStepVariablesList()
{
    if (IsSubModelPart())
        KRATOS_THROW_ERROR(std::logic_error, "Calling the method of the sub model part ", Name())
        //KRATOS_ERROR << "Calling the SetNodalSolutionStepVariablesList method of the sub model part " << Name()
        //			<< " please call the one of the parent modelpart : " << mpParentModelPart->Name() << std::endl;

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
        if( &(*existing_prop_it) != pNewProperties.get() )
        {
            KRATOS_ERROR << "trying to add a property with existing Id within the model part : " << Name() << ", property Id is :" << pNewProperties->Id();
        }
    }
    else
    {
        GetMesh(ThisIndex).AddProperties(pNewProperties);
    }
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
            if(&(*existing_element_it) != (pNewElement.get()))//check if the pointee coincides
                KRATOS_ERROR << "attempting to add pNewElement with Id :" << pNewElement->Id() << ", unfortunately a (different) element with the same Id already exists" << std::endl;
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

            current_part = current_part->GetParentModelPart();
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
    if(existing_element_iterator != GetMesh(ThisIndex).ElementsEnd() )
        KRATOS_ERROR << "trying to construct an element with ID " << Id << " however an element with the same Id already exists";


    //create the new element
    ElementType const& r_clone_element = KratosComponents<ElementType>::Get(ElementName);
    Element::Pointer p_element = r_clone_element.Create(Id, pElementNodes, pProperties);

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

void ModelPart::RemoveElements(Flags identifier_flag)
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

            if( i_elem->IsNot(identifier_flag) )
                erase_count++;
        }

        ModelPart::ElementsContainerType temp_elements_container;
        temp_elements_container.reserve(i_mesh->Elements().size() - erase_count);

        temp_elements_container.swap(i_mesh->Elements());

        for(ModelPart::ElementsContainerType::iterator i_elem = temp_elements_container.begin() ; i_elem != temp_elements_container.end() ; i_elem++)
        {
            if( i_elem->IsNot(identifier_flag) )
                (i_mesh->Elements()).push_back(std::move(*(i_elem.base())));
        }
    }

    //now recursively remove the elements in the submodelparts
    for (SubModelPartIterator i_sub_model_part = SubModelPartsBegin(); i_sub_model_part != SubModelPartsEnd(); i_sub_model_part++)
        i_sub_model_part->RemoveElements(identifier_flag);
}

void ModelPart::RemoveElementsFromAllLevels(Flags identifier_flag)
{
    ModelPart& root_model_part = GetRootModelPart();
    root_model_part.RemoveElements(identifier_flag);
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
            if(&(*existing_constraint_it) != (pNewMasterSlaveConstraint.get()))//check if the pointee coincides
            {
                KRATOS_ERROR << "attempting to add Master-Slave constraint with Id :" << pNewMasterSlaveConstraint->Id() << ", unfortunately a (different) condition with the same Id already exists" << std::endl;
            }
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

            current_part = current_part->GetParentModelPart();
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
    if(existing_constraint_iterator != GetMesh(ThisIndex).MasterSlaveConstraintsEnd() )
        KRATOS_ERROR << "trying to construct an master-slave constraint with ID " << Id << " however a constraint with the same Id already exists";


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

        if(GetMesh(ThisIndex).HasMasterSlaveConstraint(Id))
            KRATOS_ERROR << "trying to construct an master-slave constraint with ID " << Id << " however a constraint with the same Id already exists";


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


ModelPart::MasterSlaveConstraintType::Pointer ModelPart::CreateNewMasterSlaveConstraint(const std::string& ConstraintName,
                                                                                    ModelPart::IndexType Id,
                                                                                    ModelPart::NodeType& rMasterNode,
                                                                                    const ModelPart::VariableComponentType& rMasterVariable,
                                                                                    ModelPart::NodeType& rSlaveNode,
                                                                                    const ModelPart::VariableComponentType& rSlaveVariable,
                                                                                    double Weight,
                                                                                    double Constant,
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

        auto existing_constraint_iterator = GetMesh(ThisIndex).MasterSlaveConstraints().find(Id);
        if(existing_constraint_iterator != GetMesh(ThisIndex).MasterSlaveConstraintsEnd() )
            KRATOS_ERROR << "trying to construct an master-slave constraint with ID " << Id << " however a constraint with the same Id already exists";


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
            if(&(*existing_condition_it) != (pNewCondition.get()))//check if the pointee coincides
                KRATOS_ERROR << "attempting to add pNewCondition with Id :" << pNewCondition->Id() << ", unfortunately a (different) condition with the same Id already exists" << std::endl;
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

            current_part = current_part->GetParentModelPart();
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
    Geometry< Node < 3 > >::PointsArrayType pConditionNodes;

    for (unsigned int i = 0; i < ConditionNodeIds.size(); i++)
    {
        pConditionNodes.push_back(pGetNode(ConditionNodeIds[i]));
    }

    return CreateNewCondition(ConditionName, Id, pConditionNodes, pProperties, ThisIndex);
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
    if(existing_condition_iterator != GetMesh(ThisIndex).ConditionsEnd() )
    {
        KRATOS_ERROR << "trying to construct a condition with ID " << Id << " however a condition with the same Id already exists";

//         if(pConditionNodes.size() != existing_condition_iterator.size())
//             KRATOS_ERROR << "trying to construct a condition with a different number of nodes. Element Id is ";
//
//         //check if the ids of nodes coincides
//         bool nodes_are_matching = true;
//         for(unsigned int i=0; i<pConditionNodes.size(); i++)
//         {
//
//             if(pConditionNodes[i].Id() != existing_condition_iterator->GetGeometry()[i].Id())
//             {
//                 nodes_are_matching = false;
//                 break;
//             }
//         }
//         if(nodes_are_matching == false)
//             KRATOS_THROW_ERROR << "Error when attempting to construct condition with Id" << Id << " a condition with the same Id exists with the following geometry " << existing_condition_iterator->GetGeometry();
//
//         return *(existing_condition_iterator.base());
    }

    //get the element
    ConditionType const& r_clone_condition = KratosComponents<ConditionType>::Get(ConditionName);
    ConditionType::Pointer p_condition = r_clone_condition.Create(Id, pConditionNodes, pProperties);

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

void ModelPart::RemoveConditions(Flags identifier_flag)
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

            if( i_cond->IsNot(identifier_flag) )
                erase_count++;
        }

        ModelPart::ConditionsContainerType temp_conditions_container;
        temp_conditions_container.reserve(i_mesh->Conditions().size() - erase_count);

        temp_conditions_container.swap(i_mesh->Conditions());

        for(ModelPart::ConditionsContainerType::iterator i_cond = temp_conditions_container.begin() ; i_cond != temp_conditions_container.end() ; i_cond++)
        {
            if( i_cond->IsNot(identifier_flag) )
                (i_mesh->Conditions()).push_back(std::move(*(i_cond.base())));
        }
    }

    //now recursively remove the conditions in the submodelparts
    for (SubModelPartIterator i_sub_model_part = SubModelPartsBegin(); i_sub_model_part != SubModelPartsEnd(); i_sub_model_part++)
        i_sub_model_part->RemoveConditions(identifier_flag);
}

void ModelPart::RemoveConditionsFromAllLevels(Flags identifier_flag)
{
    ModelPart& root_model_part = GetRootModelPart();
    root_model_part.RemoveConditions(identifier_flag);
}


ModelPart&  ModelPart::CreateSubModelPart(std::string const& NewSubModelPartName)
{
    if (mSubModelParts.find(NewSubModelPartName) == mSubModelParts.end())
    {
        Kratos::shared_ptr<ModelPart>  p_model_part = Kratos::make_shared<ModelPart>(NewSubModelPartName);
        p_model_part->SetParentModelPart(this);
        delete p_model_part->mpVariablesList;
        p_model_part->mpVariablesList = mpVariablesList;
        p_model_part->mBufferSize = this->mBufferSize;
        p_model_part->mpProcessInfo = this->mpProcessInfo;
        mSubModelParts.insert(p_model_part);
        return *p_model_part;
    }
    else
        KRATOS_THROW_ERROR(std::logic_error, "There is an already existing sub model part with name ", NewSubModelPartName)
        // Here a warning would be enough. To be disscussed. Pooyan.
        //KRATOS_ERROR << "There is an already existing sub model part with name \"" << NewSubModelPartName << "\" in model part: \"" << Name() << "\"" << std::endl;
    }

void ModelPart::AddSubModelPart(Kratos::shared_ptr<ModelPart> pThisSubModelPart)
{
   KRATOS_ERROR << "cannot add a submodelpart, since submodelparts are univocally owned by their father " << std::endl;
}
/** Remove a sub modelpart with given name.
*/
void  ModelPart::RemoveSubModelPart(std::string const& ThisSubModelPartName)
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
void  ModelPart::RemoveSubModelPart(ModelPart& ThisSubModelPart)
{
    std::string name = ThisSubModelPart.Name();
    // finding the sub model part
    SubModelPartIterator i_sub_model_part = mSubModelParts.find(name);

    if (i_sub_model_part == mSubModelParts.end())
        KRATOS_THROW_ERROR(std::logic_error, "The sub modelpart does not exist", "")
        //KRATOS_ERROR << "The sub modelpart  \"" << name << "\" does not exist in the \"" << Name() << "\" model part to be removed" << std::endl;


    mSubModelParts.erase(name);
}

std::vector<std::string> ModelPart::GetSubModelPartNames()
{
    std::vector<std::string> SubModelPartsNames;

    for(SubModelPartIterator i_sub_model_part = mSubModelParts.begin(); i_sub_model_part != mSubModelParts.end(); i_sub_model_part++)
    {
        SubModelPartsNames.push_back(i_sub_model_part->Name());
    }

    return SubModelPartsNames;
}

void ModelPart::SetBufferSize(ModelPart::IndexType NewBufferSize)
{
    if (IsSubModelPart())
        KRATOS_THROW_ERROR(std::logic_error, "Calling the method of the sub model part ", Name())
        //	KRATOS_ERROR << "Calling the SetBufferSize method of the sub model part " << Name()
        //	<< " please call the one of the parent modelpart : " << mpParentModelPart->Name() << std::endl;

    for(SubModelPartIterator i_sub_model_part = mSubModelParts.begin(); i_sub_model_part != mSubModelParts.end(); i_sub_model_part++)
        {
            i_sub_model_part->mBufferSize = NewBufferSize;
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

/// run input validation
int ModelPart::Check(ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY
    int err = 0;
    for (ElementConstantIterator elem_iterator = ElementsBegin(); elem_iterator != ElementsEnd(); elem_iterator++)
        err = elem_iterator->Check(rCurrentProcessInfo);
    for (ConditionConstantIterator condition_iterator = ConditionsBegin(); condition_iterator != ConditionsEnd(); condition_iterator++)
        err = condition_iterator->Check(rCurrentProcessInfo);
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
    if (!IsSubModelPart())
        rOStream  << "    Buffer Size : " << mBufferSize << std::endl;
    rOStream << "    Number of tables : " << NumberOfTables() << std::endl;
    rOStream << "    Number of sub model parts : " << NumberOfSubModelParts() << std::endl;
    if (!IsSubModelPart())
        mpProcessInfo->PrintData(rOStream);
    rOStream << std::endl;
    for (IndexType i = 0; i < mMeshes.size(); i++)
    {
        rOStream << "    Mesh " << i << " :" << std::endl;
        GetMesh(i).PrintData(rOStream, "    ");
    }
    rOStream << std::endl;
    for (SubModelPartConstantIterator i_sub_model_part = SubModelPartsBegin(); i_sub_model_part != SubModelPartsEnd(); i_sub_model_part++)
    {
        i_sub_model_part->PrintInfo(rOStream, "    ");
        rOStream << std::endl;
        i_sub_model_part->PrintData(rOStream, "    ");
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
    if (!IsSubModelPart())
        rOStream << PrefixString << "    Buffer Size : " << mBufferSize << std::endl;
    rOStream << PrefixString << "    Number of tables : " << NumberOfTables() << std::endl;
    rOStream << PrefixString << "    Number of sub model parts : " << NumberOfSubModelParts() << std::endl;
    if (!IsSubModelPart())
        mpProcessInfo->PrintData(rOStream);
    rOStream << std::endl;
    for (IndexType i = 0; i < mMeshes.size(); i++)
    {
        rOStream << PrefixString << "    Mesh " << i << " :" << std::endl;
        GetMesh(i).PrintData(rOStream, PrefixString + "    ");
    }

    for (SubModelPartConstantIterator i_sub_model_part = SubModelPartsBegin(); i_sub_model_part != SubModelPartsEnd(); i_sub_model_part++)
    {
        i_sub_model_part->PrintInfo(rOStream, PrefixString + "    ");
        rOStream << std::endl;
        i_sub_model_part->PrintData(rOStream, PrefixString + "    ");
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
    //const VariablesList* p_list = &mVariablesList;
    // I'm saving it as pointer so the nodes pointers will point to it as stored pointer. Pooyan.
    rSerializer.save("Variables List", mpVariablesList);
    rSerializer.save("Meshes", mMeshes);
    rSerializer.save("SubModelParts", mSubModelParts);
}

void ModelPart::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DataValueContainer);
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Flags );
    rSerializer.load("Name", mName);
    rSerializer.load("Buffer Size", mBufferSize);
    rSerializer.load("ProcessInfo", mpProcessInfo);
    rSerializer.load("Tables", mTables);
    //VariablesList* p_list = &mVariablesList;
    rSerializer.load("Variables List", mpVariablesList);
    rSerializer.load("Meshes", mMeshes);
    rSerializer.load("SubModelParts", mSubModelParts);

    for (SubModelPartIterator i_sub_model_part = SubModelPartsBegin(); i_sub_model_part != SubModelPartsEnd(); i_sub_model_part++)
        i_sub_model_part->SetParentModelPart(this);
}


/// input stream function
//	inline std::istream & operator >>(std::istream& rIStream,
//		ModelPart& rThis)
//	{
//		return rIStream;
//	}

//	/// output stream function
//	inline std::ostream & operator <<(std::ostream& rOStream,
//		const ModelPart& rThis)
//	{
//		rThis.PrintInfo(rOStream);
//		rOStream << std::endl;
//		rThis.PrintData(rOStream);

//		return rOStream;
//	}


}  // namespace Kratos.
