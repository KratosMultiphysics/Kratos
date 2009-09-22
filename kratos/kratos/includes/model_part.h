/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
 */

//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: nelson $
//   Date:                $Date: 2008-12-09 15:23:36 $
//   Revision:            $Revision: 1.13 $
//
//


#if !defined(KRATOS_MODEL_PART_H_INCLUDED )
#define  KRATOS_MODEL_PART_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <sstream>
#include <cstddef>


// External includes 
/* #include <boost/serialization/base_object.hpp> */
/* #include <boost/serialization/utility.hpp> */
/* #include <boost/serialization/list.hpp> */



// Project includes
#include "includes/define.h"
#include "includes/process_info.h"
#include "containers/data_value_container.h"
//#include "containers/fix_data_value_container.h"
#include "includes/mesh.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/communicator.h"


namespace Kratos
{

    ///@name Kratos Globals
    ///@{

    ///@}
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name  Enum's
    ///@{

    ///@}
    ///@name  Functions
    ///@{

    ///@}
    ///@name Kratos Classes
    ///@{

    /// ModelPart class.

    /** Detail class definition.
     */
    class ModelPart : public DataValueContainer
    {

        struct GlobalIndex
        {
        };
    public:
        ///@name  Enum's
        ///@{

        enum OwnershipType
        {
            Kratos_All,
            Kratos_Local,
            Kratos_Ghost,
            Kratos_Ownership_Size
        };

        ///@}
        ///@name Type Definitions
        ///@{

        /// Pointer definition of ModelPart
        KRATOS_CLASS_POINTER_DEFINITION(ModelPart);

        typedef unsigned int IndexType;

        typedef unsigned int SizeType;

        typedef Node < 3 > NodeType;
        typedef Properties PropertiesType;
        typedef Element ElementType;
        typedef Condition ConditionType;

        typedef Mesh<NodeType, PropertiesType, ElementType, ConditionType> MeshType;

        typedef PointerVector<MeshType> MeshesContainerType;

        /// Nodes container. Which is a vector set of nodes with their Id's as key.
        typedef MeshType::NodesContainerType NodesContainerType;

        /** Iterator over the nodes. This iterator is an indirect
            iterator over Node::Pointer which turn back a reference to
            node by * operator and not a pointer for more convenient
            usage. */
        typedef MeshType::NodeIterator NodeIterator;

        /** Const iterator over the nodes. This iterator is an indirect
            iterator over Node::Pointer which turn back a reference to
            node by * operator and not a pointer for more convenient
            usage. */
        typedef MeshType::NodeConstantIterator NodeConstantIterator;

        /** Iterator over the properties. This iterator is an indirect
            iterator over Properties::Pointer which turn back a reference to
            properties by * operator and not a pointer for more convenient
            usage. */

        /// Properties container. Which is a vector set of Properties with their Id's as key.
        typedef MeshType::PropertiesContainerType PropertiesContainerType;

        /** Iterator over the Properties. This iterator is an indirect
            iterator over Node::Pointer which turn back a reference to
            node by * operator and not a pointer for more convenient
            usage. */
        typedef MeshType::PropertiesIterator PropertiesIterator;

        /** Const iterator over the Properties. This iterator is an indirect
            iterator over Properties::Pointer which turn back a reference to
            Properties by * operator and not a pointer for more convenient
            usage. */
        typedef MeshType::PropertiesConstantIterator PropertiesConstantIterator;

        /** Iterator over the properties. This iterator is an indirect
            iterator over Properties::Pointer which turn back a reference to
            properties by * operator and not a pointer for more convenient
            usage. */

        /// Element container. A vector set of Elements with their Id's as key.
        typedef MeshType::ElementsContainerType ElementsContainerType;

        /** Iterator over the Elements. This iterator is an indirect
            iterator over Elements::Pointer which turn back a reference to
            Element by * operator and not a pointer for more convenient
            usage. */
        typedef MeshType::ElementIterator ElementIterator;

        /** Const iterator over the Elements. This iterator is an indirect
            iterator over Elements::Pointer which turn back a reference to
            Element by * operator and not a pointer for more convenient
            usage. */
        typedef MeshType::ElementConstantIterator ElementConstantIterator;

        /// Condintions container. A vector set of Conditions with their Id's as key.
        typedef MeshType::ConditionsContainerType ConditionsContainerType;

        /** Iterator over the Conditions. This iterator is an indirect
           iterator over Conditions::Pointer which turn back a reference to
           Condition by * operator and not a pointer for more convenient
           usage. */
        typedef MeshType::ConditionIterator ConditionIterator;

        /** Const iterator over the Conditions. This iterator is an indirect
            iterator over Conditions::Pointer which turn back a reference to
            Condition by * operator and not a pointer for more convenient
            usage. */
        typedef MeshType::ConditionConstantIterator ConditionConstantIterator;

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.

        ModelPart()
        : mBufferSize(1)
        , mCurrentIndex(0)
        , mProcessInfo()
        , mIndices(1, 0)
        , mpCommunicator(new Communicator)
        {
            mName = "Default";
            MeshType mesh;
            for (IndexType i = 0; i < mBufferSize; i++)
                mMeshes.push_back(mesh.Clone());
        }

        ModelPart(std::string const& NewName)
        : mBufferSize(1)
        , mCurrentIndex(0)
        , mProcessInfo()
        , mIndices(1, 0)
        , mpCommunicator(new Communicator)
        {
            mName = NewName;
            MeshType mesh;
            for (IndexType i = 0; i < mBufferSize; i++)
                mMeshes.push_back(mesh.Clone());
        }

        ModelPart(std::string const& NewName, IndexType NewBufferSize)
        : mBufferSize(NewBufferSize)
        , mCurrentIndex(0)
        , mProcessInfo()
        , mIndices(NewBufferSize, 0)
        , mpCommunicator(new Communicator)
        {
            mName = NewName;
            MeshType mesh;
            for (IndexType i = 0; i < mBufferSize; i++)
                mMeshes.push_back(mesh.Clone());
        }

        /// Copy constructor.

        ModelPart(ModelPart const& rOther)
        : mName(rOther.mName)
        , mBufferSize(rOther.mBufferSize)
        , mCurrentIndex(rOther.mCurrentIndex)
        , mProcessInfo(rOther.mProcessInfo)
        , mIndices(rOther.mIndices)
        , mMeshes(rOther.mMeshes)
        , mVariablesList(rOther.mVariablesList)
        , mpCommunicator(rOther.mpCommunicator)
        {
        }


        /// Destructor.

        virtual ~ModelPart()
        {
            for (NodeIterator i_node = NodesBegin(); i_node != NodesEnd(); i_node++)
            {
                if (i_node->pGetVariablesList() == &mVariablesList)
                    i_node->ClearSolutionStepsData();
            }
            // 	  for(IndexType i = 0 ; i < mBufferSize ; ++i)
            // 	    RemoveSolutionStepData(mIndices[i], mMeshes[i]);
        }


        ///@}
        ///@name Operators
        ///@{

        /// Assignment operator.

        ModelPart & operator=(ModelPart const& rOther)
        {
            mName = rOther.mName;
            mBufferSize = rOther.mBufferSize;
            mCurrentIndex = rOther.mCurrentIndex;
            mProcessInfo = rOther.mProcessInfo;
            mIndices = rOther.mIndices;
            mMeshes = rOther.mMeshes;

            return *this;
        }

        ///@}
        ///@name Solution Steps
        ///@{

        //       IndexType CreateSolutionStep()
        // 	{
        // 	  IndexType new_index = Counter<GlobalIndex>::Increment();
        // 	  mCurrentIndex++;

        // 	  if(mCurrentIndex >= mBufferSize)
        // 	    mCurrentIndex = 0;

        // 	  RemoveSolutionStepData(mIndices[mCurrentIndex], mMeshes[mCurrentIndex]);

        // 	  mMeshes(mCurrentIndex) = MeshType::Pointer(new MeshType);
        // 	  mIndices[mCurrentIndex] = new_index;

        // 	  mProcessInfo.CreateSolutionStepInfo(new_index);
        // 	  mProcessInfo.ClearHistory(mBufferSize);

        // 	  return new_index;
        // 	}

        IndexType CreateSolutionStep()
        {
            KRATOS_ERROR(std::logic_error, "This method needs updating and is not working. Pooyan", "")
            IndexType new_index = Counter<GlobalIndex>::Increment();
            mCurrentIndex++;

            if (mCurrentIndex >= mBufferSize)
                mCurrentIndex = 0;

            mMeshes(mCurrentIndex) = MeshType::Pointer(new MeshType);
            mIndices[mCurrentIndex] = new_index;

            mProcessInfo.CreateSolutionStepInfo(new_index);
            mProcessInfo.ReIndexBuffer(mBufferSize);
            mProcessInfo.ClearHistory(mBufferSize);

            return new_index;
        }

        IndexType CloneSolutionStep()
        {
            // 	  IndexType new_index = Counter<GlobalIndex>::Increment();
            //new_index = 0; //CHAPUZA!!!! TODO!! controlar!!

            // 	  IndexType old_index = mIndices[mCurrentIndex];
            // 	  IndexType current_index = mCurrentIndex + 1;
            // 	  if(current_index >= mBufferSize)
            // 	    current_index = 0;

            for (NodeIterator node_iterator = NodesBegin(); node_iterator != NodesEnd(); node_iterator++)
                node_iterator->CloneSolutionStepData();
            /* 	  CopySolutionStepData(new_index, old_index); */

            // 	  mMeshes(current_index) = mMeshes(mCurrentIndex);
            // 	  mIndices[current_index] = new_index;
            mCurrentIndex++;
            mProcessInfo.CloneSolutionStepInfo();
            // //	  mProcessInfo.ReIndexBuffer(mBufferSize);

            mProcessInfo.ClearHistory(mBufferSize);


            // 	  mCurrentIndex = current_index;

            return 0;
        }

        // This method works better than clone but just for unchanged mesh.
        //       IndexType OverwriteSolutionStep()
        // 	{
        // 	  IndexType new_index = Counter<GlobalIndex>::Increment();
        // 	  IndexType old_index = mIndices[mCurrentIndex];
        // 	  IndexType current_index = mCurrentIndex + 1;

        // 	  if(current_index >= mBufferSize)
        // 	    current_index = 0;

        // 	  if(!mMeshes[current_index].Nodes().empty())
        // 	    OverwriteSolutionStepData(new_index, old_index, mIndices[current_index]);
        // 	  else
        // 	    CopySolutionStepData(new_index, old_index);

        // 	  mMeshes[current_index] = mMeshes[mCurrentIndex];
        // 	  mIndices[current_index] = new_index;

        // 	  mProcessInfo.CloneSolutionStepInfo(new_index);
        // 	  mProcessInfo.ClearHistory(mBufferSize);

        // 	  mCurrentIndex = current_index;

        // 	  return new_index;
        // 	}

        // commented due to a bug, Pooyan.
        //       IndexType CreateTimeStep()
        // 	{
        // 	  IndexType new_index = CreateSolutionStep();
        // 	  mProcessInfo.SetAsTimeStepInfo();

        // 	  return new_index;
        // 	}

        IndexType CloneTimeStep()
        {
            IndexType new_index = CloneSolutionStep();
            mProcessInfo.SetAsTimeStepInfo();

            return new_index;
        }

        //       IndexType OverwriteTimeStep()
        // 	{
        // 	  IndexType new_index = OverwriteSolutionStep();
        // 	  mProcessInfo.SetAsTimeStepInfo();

        // 	  return new_index;
        // 	}

        IndexType CreateTimeStep(double NewTime)
        {
            IndexType new_index = CreateSolutionStep();
            mProcessInfo.SetAsTimeStepInfo(NewTime);

            return new_index;
        }

        IndexType CloneTimeStep(double NewTime)
        {
            IndexType new_index = CloneSolutionStep();
            mProcessInfo.SetAsTimeStepInfo(NewTime);


            return new_index;
        }

        //       IndexType OverwriteTimeStep(double NewTime)
        // 	{
        // 	  IndexType new_index = OverwriteSolutionStep();
        // 	  mProcessInfo.SetAsTimeStepInfo(NewTime);

        // 	  return new_index;
        // 	}


        //       void CopySolutionStepData(IndexType SolutionStepIndex, IndexType SourceSolutionStepIndex)
        //       {
        // 	  for(NodeIterator node_iterator = NodesBegin() ; node_iterator != NodesEnd() ; node_iterator++)
        // 	    node_iterator->CloneSolutionStepData(SourceSolutionStepIndex);
        //       }

        void OverwriteSolutionStepData(IndexType SourceSolutionStepIndex, IndexType DestinationSourceSolutionStepIndex)
        {
            for (NodeIterator node_iterator = NodesBegin(); node_iterator != NodesEnd(); node_iterator++)
                node_iterator->OverwriteSolutionStepData(SourceSolutionStepIndex, DestinationSourceSolutionStepIndex);
        }

        void ReduceTimeStep(ModelPart& rModelPart, double NewTime)
        {
            KRATOS_TRY

            //ATTENTION: this function does not touch the coordinates of the nodes.
            //It just resets the database values to the values at the beginning of the time step
            rModelPart.OverwriteSolutionStepData(1, 0);
            rModelPart.GetProcessInfo().SetCurrentTime(NewTime);

            KRATOS_CATCH("error in reducing the time step")

        }

        ///@}
        ///@name Nodes
        ///@{

        SizeType NumberOfNodes(IndexType ThisIndex = 0) const
        {
            return GetMesh(ThisIndex).NumberOfNodes();
        }

        /** Inserts a node in the current mesh.
         */
        void AddNode(NodeType::Pointer pNewNode, IndexType ThisIndex = 0)
        {
            GetMesh(ThisIndex).AddNode(pNewNode);
        }

        /** Inserts a node in the current mesh.
         */
        NodeType::Pointer CreateNewNode(int Id, double x, double y, double z, VariablesList* pNewVariablesList, IndexType ThisIndex = 0)
        {
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

        NodeType::Pointer CreateNewNode(IndexType Id, double x, double y, double z, IndexType ThisIndex = 0)
        {
            //create a new node
            NodeType::Pointer p_new_node = NodeType::Pointer(new NodeType(Id, x, y, z));

            // Giving model part's variables list to the node
            p_new_node->SetSolutionStepVariablesList(&mVariablesList);

            //set buffer size
            p_new_node->SetBufferSize(mBufferSize);

            //add the new node to the list of nodes
            GetMesh(ThisIndex).AddNode(p_new_node);

            return p_new_node;
        }

        NodeType::Pointer CreateNewNode(IndexType Id, double x, double y, double z, double* pThisData, IndexType ThisIndex = 0)
        {
            //create a new node
            NodeType::Pointer p_new_node = NodeType::Pointer(new NodeType(Id, x, y, z, &mVariablesList, pThisData, mBufferSize));

            //add the new node to the list of nodes
            GetMesh(ThisIndex).AddNode(p_new_node);

            return p_new_node;

        }

        NodeType::Pointer CreateNewNode(IndexType NodeId, NodeType const& rSourceNode, IndexType ThisIndex = 0)
        {
            //create a new node
            NodeType::Pointer p_new_node = NodeType::Pointer(new NodeType(NodeId, rSourceNode.X(), rSourceNode.Y(), rSourceNode.Z()));

            // Giving model part's variables list to the node
            p_new_node->SetSolutionStepVariablesList(&mVariablesList);

            //set buffer size
            p_new_node->SetBufferSize(mBufferSize);

            //add the new node to the list of nodes
            GetMesh(ThisIndex).AddNode(p_new_node);

            return p_new_node;

        }

        void AssignNode(NodeType::Pointer pThisNode, IndexType ThisIndex = 0)
        {
            // Giving model part's variables list to the node
            pThisNode->SetSolutionStepVariablesList(&mVariablesList);

            //set buffer size
            pThisNode->SetBufferSize(mBufferSize);

            //add the new node to the list of nodes
            GetMesh(ThisIndex).AddNode(pThisNode);

        }

        /** Returns the Node::Pointer  corresponding to it's identifier */
        NodeType::Pointer pGetNode(IndexType NodeId, IndexType ThisIndex = 0)
        {
            return GetMesh(ThisIndex).pGetNode(NodeId);
        }

        /** Returns a reference node corresponding to it's identifier */
        NodeType& GetNode(IndexType NodeId, IndexType ThisIndex = 0)
        {
            return GetMesh(ThisIndex).GetNode(NodeId);
        }

        /** Remove the node with given Id from current mesh.
         */
        void RemoveNode(IndexType NodeId, IndexType ThisIndex = 0)
        {
            GetMesh(ThisIndex).RemoveNode(NodeId);
        }

        /** Remove given node from current mesh.
         */
        void RemoveNode(NodeType& ThisNode, IndexType ThisIndex = 0)
        {
            GetMesh(ThisIndex).RemoveNode(ThisNode);
        }

        /** Remove given node from current mesh.
         */
        void RemoveNode(NodeType::Pointer pThisNode, IndexType ThisIndex = 0)
        {
            GetMesh(ThisIndex).RemoveNode(pThisNode);
        }

        NodeIterator NodesBegin(IndexType ThisIndex = 0)
        {
            return GetMesh(ThisIndex).NodesBegin();
        }

        NodeConstantIterator NodesBegin(IndexType ThisIndex = 0) const
        {
            return GetMesh(ThisIndex).NodesBegin();
        }

        NodeIterator NodesEnd(IndexType ThisIndex = 0)
        {
            return GetMesh(ThisIndex).NodesEnd();
        }

        NodeConstantIterator NodesEnd(IndexType ThisIndex = 0) const
        {
            return GetMesh(ThisIndex).NodesEnd();
        }

        NodesContainerType& Nodes(IndexType ThisIndex = 0)
        {
            return GetMesh(ThisIndex).Nodes();
        }

        NodesContainerType::Pointer pNodes(IndexType ThisIndex = 0)
        {
            return GetMesh(ThisIndex).pNodes();
        }

        void SetNodes(NodesContainerType::Pointer pOtherNodes, IndexType ThisIndex = 0)
        {
            GetMesh(ThisIndex).SetNodes(pOtherNodes);
        }

        NodesContainerType::ContainerType& NodesArray(IndexType ThisIndex = 0)
        {
            return GetMesh(ThisIndex).NodesArray();
        }

        template<class TDataType>
        void AddNodalSolutionStepVariable(Variable<TDataType> const& ThisVariable)
        {
            mVariablesList.Add(ThisVariable);
        }

        VariablesList& GetNodalSolutionStepVariablesList()
        {
            return mVariablesList;
        }

        void SetNodalSolutionStepVariablesList()
        {
            for (NodeIterator i_node = NodesBegin(); i_node != NodesEnd(); ++i_node)
                i_node->SetSolutionStepVariablesList(&mVariablesList);
        }

        SizeType GetNodalSolutionStepDataSize()
        {
            return mVariablesList.DataSize();
        }

        SizeType GetNodalSolutionStepTotalDataSize()
        {
            return mVariablesList.DataSize() * mBufferSize;
        }


        ///@}
        ///@name Properties
        ///@{

        SizeType NumberOfProperties(IndexType ThisIndex = 0) const
        {
            return GetMesh(ThisIndex).NumberOfProperties();
        }

        /** Inserts a properties in the current mesh.
         */
        void AddProperties(PropertiesType::Pointer pNewProperties, IndexType ThisIndex = 0)
        {
            GetMesh(ThisIndex).AddProperties(pNewProperties);
        }

        /** Returns the Properties::Pointer  corresponding to it's identifier */
        PropertiesType::Pointer pGetProperties(IndexType PropertiesId, IndexType ThisIndex = 0)
        {
            return GetMesh(ThisIndex).pGetProperties(PropertiesId);
        }

        /** Returns a reference Properties corresponding to it's identifier */
        PropertiesType& GetProperties(IndexType PropertiesId, IndexType ThisIndex = 0)
        {
            return GetMesh(ThisIndex).GetProperties(PropertiesId);
        }

        /** Remove the Properties with given Id from current mesh.
         */
        void RemoveProperties(IndexType PropertiesId, IndexType ThisIndex = 0)
        {
            GetMesh(ThisIndex).RemoveProperties(PropertiesId);
        }

        /** Remove given Properties from current mesh.
         */
        void RemoveProperties(PropertiesType& ThisProperties, IndexType ThisIndex = 0)
        {
            GetMesh(ThisIndex).RemoveProperties(ThisProperties);
        }

        /** Remove given Properties from current mesh.
         */
        void RemoveProperties(PropertiesType::Pointer pThisProperties, IndexType ThisIndex = 0)
        {
            GetMesh(ThisIndex).RemoveProperties(pThisProperties);
        }

        PropertiesIterator PropertiesBegin(IndexType ThisIndex = 0)
        {
            return GetMesh(ThisIndex).PropertiesBegin();
        }

        PropertiesConstantIterator PropertiesBegin(IndexType ThisIndex = 0) const
        {
            return GetMesh(ThisIndex).PropertiesBegin();
        }

        PropertiesIterator PropertiesEnd(IndexType ThisIndex = 0)
        {
            return GetMesh(ThisIndex).PropertiesEnd();
        }

        PropertiesConstantIterator PropertiesEnd(IndexType ThisIndex = 0) const
        {
            return GetMesh(ThisIndex).PropertiesEnd();
        }

        /**
         * temporarily renamed Properties() function because the declaration of Properties()
         * here violates the declaration of Properties() in properties.h
         * (janosch, in agreement with pooyan)
         */
        PropertiesContainerType& rProperties(IndexType ThisIndex = 0)
        {
            return GetMesh(ThisIndex).Properties();
        }

        PropertiesContainerType::Pointer pProperties(IndexType ThisIndex = 0)
        {
            return GetMesh(ThisIndex).pProperties();
        }

        void SetProperties(PropertiesContainerType::Pointer pOtherProperties, IndexType ThisIndex = 0)
        {
            GetMesh(ThisIndex).SetProperties(pOtherProperties);
        }

        PropertiesContainerType::ContainerType& PropertiesArray(IndexType ThisIndex = 0)
        {
            return GetMesh(ThisIndex).PropertiesArray();
        }

        ///@}
        ///@name Elements
        ///@{

        SizeType NumberOfElements(IndexType ThisIndex = 0) const
        {
            return GetMesh(ThisIndex).NumberOfElements();
        }

        /** Inserts a element in the current mesh.
         */
        void AddElement(ElementType::Pointer pNewElement, IndexType ThisIndex = 0)
        {
            GetMesh(ThisIndex).AddElement(pNewElement);
        }

        /** Returns the Element::Pointer  corresponding to it's identifier */
        ElementType::Pointer pGetElement(IndexType ElementId, IndexType ThisIndex = 0)
        {
            return GetMesh(ThisIndex).pGetElement(ElementId);
        }

        /** Returns a reference element corresponding to it's identifier */
        ElementType& GetElement(IndexType ElementId, IndexType ThisIndex = 0)
        {
            return GetMesh(ThisIndex).GetElement(ElementId);
        }

        /** Remove the element with given Id from current mesh.
         */
        void RemoveElement(IndexType ElementId, IndexType ThisIndex = 0)
        {
            GetMesh(ThisIndex).RemoveElement(ElementId);
        }

        /** Remove given element from current mesh.
         */
        void RemoveElement(ElementType& ThisElement, IndexType ThisIndex = 0)
        {
            GetMesh(ThisIndex).RemoveElement(ThisElement);
        }

        /** Remove given element from current mesh.
         */
        void RemoveElement(ElementType::Pointer pThisElement, IndexType ThisIndex = 0)
        {
            GetMesh(ThisIndex).RemoveElement(pThisElement);
        }

        ElementIterator ElementsBegin(IndexType ThisIndex = 0)
        {
            return GetMesh(ThisIndex).ElementsBegin();
        }

        ElementConstantIterator ElementsBegin(IndexType ThisIndex = 0) const
        {
            return GetMesh(ThisIndex).ElementsBegin();
        }

        ElementIterator ElementsEnd(IndexType ThisIndex = 0)
        {
            return GetMesh(ThisIndex).ElementsEnd();
        }

        ElementConstantIterator ElementsEnd(IndexType ThisIndex = 0) const
        {
            return GetMesh(ThisIndex).ElementsEnd();
        }

        ElementsContainerType& Elements(IndexType ThisIndex = 0)
        {
            return GetMesh(ThisIndex).Elements();
        }

        ElementsContainerType::Pointer pElements(IndexType ThisIndex = 0)
        {
            return GetMesh(ThisIndex).pElements();
        }

        void SetElements(ElementsContainerType::Pointer pOtherElements, IndexType ThisIndex = 0)
        {
            GetMesh(ThisIndex).SetElements(pOtherElements);
        }

        ElementsContainerType::ContainerType& ElementsArray(IndexType ThisIndex = 0)
        {
            return GetMesh(ThisIndex).ElementsArray();
        }


        ///@}
        ///@name Conditions
        ///@{

        SizeType NumberOfConditions(IndexType ThisIndex = 0) const
        {
            return GetMesh(ThisIndex).NumberOfConditions();
        }

        /** Inserts a condition in the current mesh.
         */
        void AddCondition(ConditionType::Pointer pNewCondition, IndexType ThisIndex = 0)
        {
            GetMesh(ThisIndex).AddCondition(pNewCondition);
        }

        /** Returns the Condition::Pointer  corresponding to it's identifier */
        ConditionType::Pointer pGetCondition(IndexType ConditionId, IndexType ThisIndex = 0)
        {
            return GetMesh(ThisIndex).pGetCondition(ConditionId);
        }

        /** Returns a reference condition corresponding to it's identifier */
        ConditionType& GetCondition(IndexType ConditionId, IndexType ThisIndex = 0)
        {
            return GetMesh(ThisIndex).GetCondition(ConditionId);
        }

        /** Remove the condition with given Id from current mesh.
         */
        void RemoveCondition(IndexType ConditionId, IndexType ThisIndex = 0)
        {
            GetMesh(ThisIndex).RemoveCondition(ConditionId);
        }

        /** Remove given condition from current mesh.
         */
        void RemoveCondition(ConditionType& ThisCondition, IndexType ThisIndex = 0)
        {
            GetMesh(ThisIndex).RemoveCondition(ThisCondition);
        }

        /** Remove given condition from current mesh.
         */
        void RemoveCondition(ConditionType::Pointer pThisCondition, IndexType ThisIndex = 0)
        {
            GetMesh(ThisIndex).RemoveCondition(pThisCondition);
        }

        ConditionIterator ConditionsBegin(IndexType ThisIndex = 0)
        {
            return GetMesh(ThisIndex).ConditionsBegin();
        }

        ConditionConstantIterator ConditionsBegin(IndexType ThisIndex = 0) const
        {
            return GetMesh(ThisIndex).ConditionsBegin();
        }

        ConditionIterator ConditionsEnd(IndexType ThisIndex = 0)
        {
            return GetMesh(ThisIndex).ConditionsEnd();
        }

        ConditionConstantIterator ConditionsEnd(IndexType ThisIndex = 0) const
        {
            return GetMesh(ThisIndex).ConditionsEnd();
        }

        ConditionsContainerType& Conditions(IndexType ThisIndex = 0)
        {
            return GetMesh(ThisIndex).Conditions();
        }

        ConditionsContainerType::Pointer pConditions(IndexType ThisIndex = 0)
        {
            return GetMesh(ThisIndex).pConditions();
        }

        void SetConditions(ConditionsContainerType::Pointer pOtherConditions, IndexType ThisIndex = 0)
        {
            GetMesh(ThisIndex).SetConditions(pOtherConditions);
        }

        ConditionsContainerType::ContainerType& ConditionsArray(IndexType ThisIndex = 0)
        {
            return GetMesh(ThisIndex).ConditionsArray();
        }

        ///@}
        ///@name Access
        ///@{

        ProcessInfo& GetProcessInfo()
        {
            return mProcessInfo;
        }

        ProcessInfo const& GetProcessInfo() const
        {
            return mProcessInfo;
        }

        void SetProcessInfo(ProcessInfo& NewProcessInfo)
        {
            mProcessInfo = NewProcessInfo;
        }

        MeshType::Pointer pGetMesh(IndexType ThisIndex = 0)
        {
            return mMeshes(ThisIndex);
        }

        const MeshType::Pointer pGetMesh(IndexType ThisIndex = 0) const
        {
            return mMeshes(ThisIndex);
        }

        MeshType& GetMesh(IndexType ThisIndex = 0)
        {
            return mMeshes[ThisIndex];
        }

        MeshType const& GetMesh(IndexType ThisIndex = 0) const
        {
            return mMeshes[ThisIndex];
        }

        MeshesContainerType& GetMeshes()
        {
            return mMeshes;
        }

        MeshesContainerType const& GetMeshes() const
        {
            return mMeshes;
        }

        std::string& Name()
        {
            return mName;
        }

        std::string const& Name() const
        {
            return mName;
        }

        Communicator& GetCommunicator()
        {
            return *mpCommunicator;
        }

        Communicator const& GetCommunicator() const
        {
            return *mpCommunicator;
        }

        void SetCommunicator(Communicator::Pointer pNewCommunicator)
        {
            mpCommunicator = pNewCommunicator;
        }

        ///@}
        ///@name Operations
        ///@{

        void SetBufferSize(IndexType NewBufferSize)
        {
            mBufferSize = NewBufferSize;

            for (NodeIterator node_iterator = NodesBegin(); node_iterator != NodesEnd(); node_iterator++)
                node_iterator->SetBufferSize(mBufferSize);

        }

        IndexType GetBufferSize()
        {
            return mBufferSize;
        }


        ///@}
        ///@name Access
        ///@{


        ///@}
        ///@name Inquiry
        ///@{


        ///@}
        ///@name Input and output
        ///@{

        /// Turn back information as a string.

        virtual std::string Info() const
        {
            return mName + " model part";
        }

        /// Print information about this object.

        virtual void PrintInfo(std::ostream& rOStream) const
        {
            rOStream << Info();
        }

        /// Print object's data.

        virtual void PrintData(std::ostream& rOStream) const
        {
            rOStream << "    Buffer Size : " << mBufferSize << std::endl;
            mProcessInfo.PrintData(rOStream);
            rOStream << std::endl;
            for (IndexType i = 0; i < mMeshes.size(); i++)
            {
                rOStream << "    Mesh " << i << " : " << std::endl;
                GetMesh(i).PrintData(rOStream);
            }
        }


        ///@}
        ///@name Friends
        ///@{


        ///@}

    protected:
        ///@name Protected static Member Variables
        ///@{


        ///@}
        ///@name Protected member Variables
        ///@{


        ///@}
        ///@name Protected Operators
        ///@{


        ///@}
        ///@name Protected Operations
        ///@{


        ///@}
        ///@name Protected  Access
        ///@{


        ///@}
        ///@name Protected Inquiry
        ///@{


        ///@}
        ///@name Protected LifeCycle
        ///@{


        ///@}

    private:
        ///@name Static Member Variables
        ///@{


        ///@}
        ///@name Member Variables
        ///@{

        std::string mName;

        IndexType mBufferSize;

        IndexType mCurrentIndex;

        ProcessInfo mProcessInfo;

        std::vector<IndexType> mIndices;

        MeshesContainerType mMeshes;

        VariablesList mVariablesList;

        Communicator::Pointer mpCommunicator;

        ///@}
        ///@name Private Operators
        ///@{


        ///@}
        ///@name Private Operations
        ///@{

        //       friend class boost::serialization::access;

        //       template<class TArchive>
        // 	  void serialize(TArchive & ThisArchive, const unsigned int ThisVersion)
        // 	  {
        // /* 	      ThisArchive & mName & mBufferSize & mCurrentIndex; */
        // 	  }

        //       void RemoveSolutionStepData(IndexType SolutionStepIndex, MeshType& ThisMesh)
        // 	{
        // 	  for(NodeIterator i_node = ThisMesh.NodesBegin() ; i_node != ThisMesh.NodesEnd() ; ++i_node)
        // 	    i_node->RemoveSolutionStepNodalData(SolutionStepIndex);
        // 	}

        ///@}
        ///@name Private  Access
        ///@{


        ///@}
        ///@name Private Inquiry
        ///@{


        ///@}
        ///@name Un accessible methods
        ///@{


        ///@}

    }; // Class ModelPart 

    ///@}

    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    /// input stream function
    inline std::istream & operator >>(std::istream& rIStream,
            ModelPart& rThis);

    /// output stream function

    inline std::ostream & operator <<(std::ostream& rOStream,
            const ModelPart& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }
    ///@}


} // namespace Kratos.

#endif // KRATOS_MODEL_PART_H_INCLUDED  defined 


