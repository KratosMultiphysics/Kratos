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
		mMeshes.push_back(mesh.Clone());
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
		mName = NewName;
		MeshType mesh;
		mMeshes.push_back(mesh.Clone());
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
		mName = NewName;
		MeshType mesh;
		mMeshes.push_back(mesh.Clone());
		mpCommunicator->SetLocalMesh(pGetMesh());  // assigning the current mesh to the local mesh of communicator for openmp cases
	}

	// Copy constructor.
	ModelPart::ModelPart(ModelPart const& rOther)
		: DataValueContainer(rOther)
		, Flags(rOther)
		, mName(rOther.mName)
		, mBufferSize(rOther.mBufferSize)
		, mpProcessInfo(rOther.mpProcessInfo)
		, mIndices(rOther.mIndices)
		, mMeshes(rOther.mMeshes)
		, mpVariablesList(new VariablesList(*rOther.mpVariablesList))
		, mpCommunicator(rOther.mpCommunicator)
		, mpParentModelPart(rOther.mpParentModelPart)
		, mSubModelParts(rOther.mSubModelParts)
	{
	}

	/// Destructor.
	ModelPart::~ModelPart()
	{
		for (NodeIterator i_node = NodesBegin(); i_node != NodesEnd(); i_node++)
		{
			if (i_node->pGetVariablesList() == mpVariablesList)
				i_node->ClearSolutionStepsData();
		}

		for (SubModelPartIterator i_sub_model_part = SubModelPartsBegin(); i_sub_model_part != SubModelPartsEnd(); i_sub_model_part++)
			delete i_sub_model_part.base()->second;

		if (!IsSubModelPart())
			delete mpVariablesList;
	}


	/// Assignment operator.
	ModelPart & ModelPart::operator=(ModelPart const& rOther)
	{
		mName = rOther.mName;
		mBufferSize = rOther.mBufferSize;
		mpProcessInfo = rOther.mpProcessInfo;
		mIndices = rOther.mIndices;
		mMeshes = rOther.mMeshes; 
		// I should not set the parent for a model part while it breaks the hierarchy. Pooyan.
		//mpParentModelPart = rOther.mpParentModelPart;
		mSubModelParts = rOther.mSubModelParts;

		//KRATOS_THROW_ERROR(std::logic_error, "This method needs updating and is not working. Pooyan", "")

		*mpVariablesList = *rOther.mpVariablesList;

		return *this;
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

		for (NodeIterator node_iterator = NodesBegin(); node_iterator != NodesEnd(); node_iterator++)
			node_iterator->CloneSolutionStepData();

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
			mpParentModelPart->AddNode(pNewNode, ThisIndex);

		GetMesh(ThisIndex).AddNode(pNewNode);
	}

	/** Inserts a node in the mesh with ThisIndex.
	*/
	ModelPart::NodeType::Pointer ModelPart::CreateNewNode(int Id, double x, double y, double z, VariablesList* pNewVariablesList, ModelPart::IndexType ThisIndex)
	{
		if (IsSubModelPart())
		{
			NodeType::Pointer p_new_node = mpParentModelPart->CreateNewNode(Id, x, y, z, pNewVariablesList, ThisIndex);
			GetMesh(ThisIndex).AddNode(p_new_node);

			return p_new_node;
		}

		//create a new node
		NodeType::Pointer p_new_node = NodeType::Pointer(new NodeType(Id, x, y, z));

		// Giving model part's variables list to the node
		p_new_node->SetSolutionStepVariablesList(pNewVariablesList);

		//set buffer size
		p_new_node->SetBufferSize(mBufferSize);

		//add the new node to the list of nodes
		GetMesh(ThisIndex).AddNode(p_new_node);

		return p_new_node;
	}

	ModelPart::NodeType::Pointer ModelPart::CreateNewNode(ModelPart::IndexType Id, double x, double y, double z, ModelPart::IndexType ThisIndex)
	{
		if (IsSubModelPart())
		{
			NodeType::Pointer p_new_node = mpParentModelPart->CreateNewNode(Id, x, y, z, ThisIndex);
			GetMesh(ThisIndex).AddNode(p_new_node);

			return p_new_node;
		}

		//create a new node
		NodeType::Pointer p_new_node = NodeType::Pointer(new NodeType(Id, x, y, z));

		// Giving model part's variables list to the node
		p_new_node->SetSolutionStepVariablesList(mpVariablesList);

		//set buffer size
		p_new_node->SetBufferSize(mBufferSize);

		//add the new node to the list of nodes
		GetMesh(ThisIndex).AddNode(p_new_node);

		return p_new_node;
	}

	ModelPart::NodeType::Pointer ModelPart::CreateNewNode(ModelPart::IndexType Id, double x, double y, double z, double* pThisData, ModelPart::IndexType ThisIndex)
	{
		if (IsSubModelPart())
		{
			NodeType::Pointer p_new_node = mpParentModelPart->CreateNewNode(Id, x, y, z, pThisData, ThisIndex);
			GetMesh(ThisIndex).AddNode(p_new_node);

			return p_new_node;
		}

		//create a new node
		NodeType::Pointer p_new_node = NodeType::Pointer(new NodeType(Id, x, y, z, mpVariablesList, pThisData, mBufferSize));

		//add the new node to the list of nodes
		GetMesh(ThisIndex).AddNode(p_new_node);

		return p_new_node;

	}

	ModelPart::NodeType::Pointer ModelPart::CreateNewNode(ModelPart::IndexType NodeId, ModelPart::NodeType const& rSourceNode, ModelPart::IndexType ThisIndex)
	{
		if (IsSubModelPart())
		{
			NodeType::Pointer p_new_node = mpParentModelPart->CreateNewNode(NodeId, rSourceNode, ThisIndex);
			GetMesh(ThisIndex).AddNode(p_new_node);

			return p_new_node;
		}

		//create a new node
		NodeType::Pointer p_new_node = NodeType::Pointer(new NodeType(NodeId, rSourceNode.X(), rSourceNode.Y(), rSourceNode.Z()));

		// Giving model part's variables list to the node
		p_new_node->SetSolutionStepVariablesList(mpVariablesList);

		//set buffer size
		p_new_node->SetBufferSize(mBufferSize);

		//add the new node to the list of nodes
		GetMesh(ThisIndex).AddNode(p_new_node);

		return p_new_node;

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
			mpParentModelPart->AddProperties(pNewProperties);

		GetMesh(ThisIndex).AddProperties(pNewProperties);
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
			mpParentModelPart->AddElement(pNewElement, ThisIndex);

		GetMesh(ThisIndex).AddElement(pNewElement);
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

		for (unsigned int i = 0; i < ElementNodeIds.size(); i++) {
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
		if (IsSubModelPart())
		{
			ElementType::Pointer p_new_element = mpParentModelPart->CreateNewElement(ElementName, Id, pElementNodes, pProperties, ThisIndex);
			GetMesh(ThisIndex).AddElement(p_new_element);
			return p_new_element;
		}

		//create the new element
		ElementType const& r_clone_element = KratosComponents<ElementType>::Get(ElementName);
		Element::Pointer p_element = r_clone_element.Create(Id, pElementNodes, pProperties);

		//add the new element
		GetMesh(ThisIndex).AddElement(p_element);

		return p_element;
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

	/** Inserts a condition in the mesh with ThisIndex.
	*/
	void ModelPart::AddCondition(ModelPart::ConditionType::Pointer pNewCondition, ModelPart::IndexType ThisIndex)
	{
		if (IsSubModelPart())
			mpParentModelPart->AddCondition(pNewCondition, ThisIndex);

		GetMesh(ThisIndex).AddCondition(pNewCondition);
	}

	/** Inserts a condition in the mesh with ThisIndex.
	*/
	ModelPart::ConditionType::Pointer ModelPart::CreateNewCondition(std::string ConditionName,
		ModelPart::IndexType Id, std::vector<IndexType> ConditionNodeIds, 
		ModelPart::PropertiesType::Pointer pProperties, ModelPart::IndexType ThisIndex)
	{
		Geometry< Node < 3 > >::PointsArrayType pConditionNodes;

		for (unsigned int i = 0; i < ConditionNodeIds.size(); i++) {
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
		if (IsSubModelPart())
		{
			ConditionType::Pointer p_new_condition = mpParentModelPart->CreateNewCondition(ConditionName, Id, pConditionNodes, pProperties, ThisIndex);
			GetMesh(ThisIndex).AddCondition(p_new_condition);
			return p_new_condition;
		}

		//get the element
		ConditionType const& r_clone_condition = KratosComponents<ConditionType>::Get(ConditionName);
		ConditionType::Pointer p_condition = r_clone_condition.Create(Id, pConditionNodes, pProperties);

		//add the new element
		GetMesh(ThisIndex).AddCondition(p_condition);

		return p_condition;
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


	ModelPart&  ModelPart::CreateSubModelPart(std::string const& NewSubModelPartName)
	{
		if (mSubModelParts.find(NewSubModelPartName) == mSubModelParts.end())
		{
			ModelPart* p_model_part = new ModelPart(NewSubModelPartName);
			p_model_part->SetParentModelPart(this);
			delete p_model_part->mpVariablesList;
			p_model_part->mpVariablesList = mpVariablesList;
                        p_model_part->mBufferSize = this->mBufferSize;
			return *(mSubModelParts.insert(p_model_part));
		}
		else
                    KRATOS_THROW_ERROR(std::logic_error, "There is an already existing sub model part with name ", NewSubModelPartName)
			// Here a warning would be enough. To be disscussed. Pooyan.
			//KRATOS_ERROR << "There is an already existing sub model part with name \"" << NewSubModelPartName << "\" in model part: \"" << Name() << "\"" << std::endl;
	}

	void ModelPart::AddSubModelPart(ModelPart& rThisSubModelPart)
	{
		if (mSubModelParts.find(rThisSubModelPart.Name()) != mSubModelParts.end())
			// Here a warning would be enough. To be disscussed. Pooyan.
			KRATOS_ERROR << "There is an already existing sub model part with name \"" << rThisSubModelPart.Name() << "\" in model part: \"" << Name() << "\"" << std::endl;
			
		if (IsSubModelPart())
		{
			mpParentModelPart->AddSubModelPart(rThisSubModelPart);
			return;
		}

		rThisSubModelPart.SetParentModelPart(this);
	}
	/** Remove a sub modelpart with given name.
	*/
	void  ModelPart::RemoveSubModelPart(std::string const& ThisSubModelPartName)
	{
		// finding the sub model part
		SubModelPartIterator i_sub_model_part = mSubModelParts.find(ThisSubModelPartName);

		if (i_sub_model_part == mSubModelParts.end())
			return; // TODO: send a warning here. Pooyan.

		// deallocate the sub model part
		delete i_sub_model_part.base()->second;

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

					// deallocate the sub model part
		delete i_sub_model_part.base()->second;

		mSubModelParts.erase(name);
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

		for (NodeIterator node_iterator = NodesBegin(); node_iterator != NodesEnd(); node_iterator++)
			node_iterator->SetBufferSize(mBufferSize);

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
		return mName + " model part";
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
            rOStream << "    Mesh " << i << " : " << std::endl;
            GetMesh(i).PrintData(rOStream, "    ");
        }

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
            rOStream << PrefixString << "    Mesh " << i << " : " << std::endl;
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
		rSerializer.save("Name", mName);
		rSerializer.save("Buffer Size", mBufferSize);
		rSerializer.save("ProcessInfo", mpProcessInfo);
		//const VariablesList* p_list = &mVariablesList;
		// I'm saving it as pointer so the nodes pointers will point to it as stored pointer. Pooyan.
		rSerializer.save("Variables List", mpVariablesList);
		rSerializer.save("Meshes", mMeshes);
	}

	void ModelPart::load(Serializer& rSerializer)
	{
		KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DataValueContainer);
		rSerializer.load("Name", mName);
		rSerializer.load("Buffer Size", mBufferSize);
		rSerializer.load("ProcessInfo", mpProcessInfo);
		//VariablesList* p_list = &mVariablesList;
		rSerializer.load("Variables List", mpVariablesList);
		rSerializer.load("Meshes", mMeshes);
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


