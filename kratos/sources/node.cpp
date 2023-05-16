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
//                   Riccardo Rossi
//

// System includes

// External includes

// Project includes
#include "includes/node.h"

namespace Kratos
{

Node::Node()
    : BaseType()
    , Flags()
    , mNodalData(0)
    , mDofs()
    , mData()
    , mInitialPosition()
    , mNodeLock()
{
    CreateSolutionStepData();
}

/***********************************************************************************/
/***********************************************************************************/

Node::Node(const IndexType NewId )
    : BaseType()
    , Flags()
    , mNodalData(NewId)
    , mDofs()
    , mData()
    , mInitialPosition()
    , mNodeLock()
{
    KRATOS_ERROR <<  "Calling the default constructor for the node ... illegal operation!!" << std::endl;
    CreateSolutionStepData();
}

/***********************************************************************************/
/***********************************************************************************/

Node::Node(const IndexType NewId, const double NewX)
    : BaseType(NewX)
    , Flags()
    , mNodalData(NewId)
    , mDofs()
    , mData()
    , mInitialPosition(NewX)
    , mNodeLock()
{
    CreateSolutionStepData();
}

/***********************************************************************************/
/***********************************************************************************/

Node::Node(const IndexType NewId, const double NewX, const double NewY)
    : BaseType(NewX, NewY)
    , Flags()
    , mNodalData(NewId)
    , mDofs()
    , mData()
    , mInitialPosition(NewX, NewY)
    , mNodeLock()
{
    CreateSolutionStepData();
}

/***********************************************************************************/
/***********************************************************************************/

Node::Node(const IndexType NewId, const double NewX, const double NewY, const double NewZ)
    : BaseType(NewX, NewY, NewZ)
    , Flags()
    , mNodalData(NewId)
    , mDofs()
    , mData()
    , mInitialPosition(NewX, NewY, NewZ)
    , mNodeLock()
{
    CreateSolutionStepData();
}

/***********************************************************************************/
/***********************************************************************************/

Node::Node(IndexType NewId, PointType const& rThisPoint)
    : BaseType(rThisPoint)
    , Flags()
    , mNodalData(NewId)
    , mDofs()
    , mData()
    , mInitialPosition(rThisPoint)
    , mNodeLock()
{
    CreateSolutionStepData();
}

/***********************************************************************************/
/***********************************************************************************/

Node::Node(const IndexType NewId, std::vector<double> const&  rOtherCoordinates)
    : BaseType(rOtherCoordinates)
    , Flags()
    , mNodalData(NewId)
    , mDofs()
    , mData()
    , mInitialPosition()
    , mNodeLock()
{
    CreateSolutionStepData();
}

/***********************************************************************************/
/***********************************************************************************/

Node::Node(const IndexType NewId, const double NewX, const double NewY, const double NewZ, VariablesList::Pointer  pVariablesList, BlockType const * ThisData, SizeType NewQueueSize)
    : BaseType(NewX, NewY, NewZ)
    , Flags()
    , mNodalData(NewId, pVariablesList,ThisData,NewQueueSize)
    , mDofs()
    , mData()
    , mInitialPosition(NewX, NewY, NewZ)
    , mNodeLock()
{
}

/***********************************************************************************/
/***********************************************************************************/

typename Node::Pointer Node::Clone()
{
    Node::Pointer p_new_node = Kratos::make_intrusive<Node>( this->Id(), (*this)[0], (*this)[1], (*this)[2]);
    p_new_node->mNodalData = this->mNodalData;

    Node::DofsContainerType& my_dofs = (this)->GetDofs();
    for (auto it_dof = my_dofs.begin(); it_dof != my_dofs.end(); it_dof++) {
        p_new_node->pAddDof(**it_dof);
    }

    p_new_node->mData = this->mData;
    p_new_node->mInitialPosition = this->mInitialPosition;

    p_new_node->Set(Flags(*this));

    return p_new_node;
}

/***********************************************************************************/
/***********************************************************************************/

Node::~Node()
{
    ClearSolutionStepsData();
}

/***********************************************************************************/
/***********************************************************************************/

Node& Node::operator=(const Node& rOther)
{
    BaseType::operator=(rOther);
    Flags::operator =(rOther);

    mNodalData = rOther.mNodalData;

    // Deep copying the dofs
    for(auto it_dof = rOther.mDofs.begin() ; it_dof != rOther.mDofs.end() ; it_dof++) {
        pAddDof(**it_dof);
    }

    mData = rOther.mData;
    mInitialPosition = rOther.mInitialPosition;

    return *this;
}

/***********************************************************************************/
/***********************************************************************************/

bool Node::operator==(const Node& rOther)
{
    return PointType::operator ==(rOther);
}

/***********************************************************************************/
/***********************************************************************************/

double& Node::operator[](const IndexType ThisIndex)
{
    return BaseType::operator[](ThisIndex);
}

/***********************************************************************************/
/***********************************************************************************/

double Node::operator[](const IndexType ThisIndex) const
{
    return BaseType::operator[](ThisIndex);
}

/***********************************************************************************/
/***********************************************************************************/

void Node::SetInitialPosition(const PointType& NewInitialPosition)
{
    mInitialPosition.X() = NewInitialPosition.X();
    mInitialPosition.Y() = NewInitialPosition.Y();
    mInitialPosition.Z() = NewInitialPosition.Z();
}

/***********************************************************************************/
/***********************************************************************************/

void Node::SetInitialPosition(
    const double X,
    const double Y,
    const double Z
    )
{
    mInitialPosition.X() = X;
    mInitialPosition.Y() = Y;
    mInitialPosition.Z() = Z;
}

/***********************************************************************************/
/***********************************************************************************/

VariablesList::Pointer Node::pGetVariablesList()
{
    return SolutionStepData().pGetVariablesList();
}

/***********************************************************************************/
/***********************************************************************************/

const VariablesList::Pointer Node::pGetVariablesList() const
{
    return SolutionStepData().pGetVariablesList();
}

/***********************************************************************************/
/***********************************************************************************/

void Node::CreateSolutionStepData()
{
    SolutionStepData().PushFront();
}

/***********************************************************************************/
/***********************************************************************************/

void Node::CloneSolutionStepData()
{
    SolutionStepData().CloneFront();
}

/***********************************************************************************/
/***********************************************************************************/

Node::IndexType Node::GetBufferSize() const
{
    return SolutionStepData().QueueSize();
}

/***********************************************************************************/
/***********************************************************************************/

void Node::SetBufferSize(IndexType NewBufferSize)
{
    SolutionStepData().Resize(NewBufferSize);
}

/***********************************************************************************/
/***********************************************************************************/

std::string Node::Info() const
{
    std::stringstream buffer;
    buffer << "Node #" << Id();
    return buffer.str();
}

/***********************************************************************************/
/***********************************************************************************/

void Node::PrintInfo(std::ostream& rOStream) const
{
    rOStream << Info();
}

/***********************************************************************************/
/***********************************************************************************/

void Node::PrintData(std::ostream& rOStream) const
{
    BaseType::PrintData(rOStream);
    if(!mDofs.empty())
        rOStream << std::endl << "    Dofs :" << std::endl;
    for(typename DofsContainerType::const_iterator i = mDofs.begin() ; i != mDofs.end() ; i++)
        rOStream << "        " << (*i)->Info() << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void Node::SortDofs(){
    std::sort(mDofs.begin(), mDofs.end(), [](Kratos::unique_ptr<DofType> const& First, Kratos::unique_ptr<DofType> const& Second)->bool{
        return First->GetVariable().Key() < Second->GetVariable().Key();
    });
}

/***********************************************************************************/
/***********************************************************************************/

void Node::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Point );
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Flags );
    rSerializer.save("NodalData", &mNodalData); // Storing it as pointer to be shared by Dof pointer
    rSerializer.save("Data", mData);
    rSerializer.save("Initial Position", mInitialPosition);
    rSerializer.save("Data", mDofs);

}

/***********************************************************************************/
/***********************************************************************************/

void Node::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Point );
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Flags );
    NodalData* p_nodal_data = &mNodalData;
    rSerializer.load("NodalData", p_nodal_data);
    rSerializer.load("Data", mData);
    rSerializer.load("Initial Position", mInitialPosition);
    rSerializer.load("Data", mDofs);
}

}  // namespace Kratos.
