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

Node::Node(IndexType NewId, const double NewX)
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

Node::Node(IndexType NewId, const double NewX, const double NewY)
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

Node::Node(IndexType NewId, const double NewX, const double NewY, const double NewZ)
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

Node::Node(IndexType NewId, Point const& rThisPoint)
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

Node::Node(IndexType NewId, std::vector<double> const&  rOtherCoordinates)
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

Node::Node(IndexType NewId, const double NewX, const double NewY, const double NewZ, VariablesList::Pointer  pVariablesList, BlockType const * ThisData, SizeType NewQueueSize)
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
    Node::Pointer p_new_node = Kratos::make_intrusive<Node >( this->Id(), (*this)[0], (*this)[1], (*this)[2]);
    p_new_node->mNodalData = this->mNodalData;

    Node::DofsContainerType& my_dofs = (this)->GetDofs();
    for (typename DofsContainerType::const_iterator it_dof = my_dofs.begin(); it_dof != my_dofs.end(); it_dof++) {
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

Node::IndexType Node::Id() const
{
    return mNodalData.Id();
}

/***********************************************************************************/
/***********************************************************************************/

Node::IndexType Node::GetId() const
{
    return mNodalData.Id();
}

/***********************************************************************************/
/***********************************************************************************/

void Node::SetId(IndexType NewId)
{
    mNodalData.SetId(NewId);
}

/***********************************************************************************/
/***********************************************************************************/

LockObject& Node::GetLock()
{
    return mNodeLock;
}

/***********************************************************************************/
/***********************************************************************************/

Node& Node::operator=(const Node& rOther)
{
    BaseType::operator=(rOther);
    Flags::operator =(rOther);

    mNodalData = rOther.mNodalData;

    // Deep copying the dofs
    for(typename DofsContainerType::const_iterator it_dof = rOther.mDofs.begin() ; it_dof != rOther.mDofs.end() ; it_dof++) {
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
    return Point::operator ==(rOther);
}

/***********************************************************************************/
/***********************************************************************************/

double& Node::operator[](IndexType ThisIndex)
{
    return BaseType::operator[](ThisIndex);
}

/***********************************************************************************/
/***********************************************************************************/

double Node::operator[](IndexType ThisIndex) const
{
    return BaseType::operator[](ThisIndex);
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

void Node::OverwriteSolutionStepData(IndexType SourceSolutionStepIndex, IndexType DestinationSourceSolutionStepIndex)
{
    SolutionStepData().AssignData(SolutionStepData().Data(SourceSolutionStepIndex), DestinationSourceSolutionStepIndex);
}

/***********************************************************************************/
/***********************************************************************************/

void Node::ClearSolutionStepsData()
{
    SolutionStepData().Clear();
}

/***********************************************************************************/
/***********************************************************************************/

void Node::SetSolutionStepVariablesList(VariablesList::Pointer pVariablesList)
{
    SolutionStepData().SetVariablesList(pVariablesList);
}

/***********************************************************************************/
/***********************************************************************************/

VariablesListDataValueContainer& Node::SolutionStepData()
{
    return mNodalData.GetSolutionStepData();
}

/***********************************************************************************/
/***********************************************************************************/

const VariablesListDataValueContainer& Node::SolutionStepData() const
{
    return mNodalData.GetSolutionStepData();
}

/***********************************************************************************/
/***********************************************************************************/

DataValueContainer& Node::Data()
{
    return mData;
}

/***********************************************************************************/
/***********************************************************************************/

DataValueContainer& Node::GetData()
{
    return mData;
}

/***********************************************************************************/
/***********************************************************************************/

const DataValueContainer& Node::GetData() const
{
    return mData;
}

/***********************************************************************************/
/***********************************************************************************/

bool Node::SolutionStepsDataHas(const VariableData& rThisVariable) const
{
    return SolutionStepData().Has(rThisVariable);
}

/***********************************************************************************/
/***********************************************************************************/

IndexType Node::GetBufferSize() const
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

const Point& Node::GetInitialPosition() const
{
    return mInitialPosition;
}

/***********************************************************************************/
/***********************************************************************************/

Point& Node::GetInitialPosition()
{
    return mInitialPosition;
}

/***********************************************************************************/
/***********************************************************************************/

double& Node::X0()
{
    return mInitialPosition.X();
}

/***********************************************************************************/
/***********************************************************************************/

double& Node::Y0()
{
    return mInitialPosition.Y();
}

/***********************************************************************************/
/***********************************************************************************/

double& Node::Z0()
{
    return mInitialPosition.Z();
}

/***********************************************************************************/
/***********************************************************************************/

double Node::X0() const
{
    return mInitialPosition.X();
}

/***********************************************************************************/
/***********************************************************************************/

double Node::Y0() const
{
    return mInitialPosition.Y();
}

/***********************************************************************************/
/***********************************************************************************/

double Node::Z0() const
{
    return mInitialPosition.Z();
}

/***********************************************************************************/
/***********************************************************************************/

void Node::SetInitialPosition(const Point& NewInitialPosition)
{
    mInitialPosition.X() = NewInitialPosition.X();
    mInitialPosition.Y() = NewInitialPosition.Y();
    mInitialPosition.Z() = NewInitialPosition.Z();
}

/***********************************************************************************/
/***********************************************************************************/

void Node::SetInitialPosition(double X,double Y, double Z)
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

Node::DofsContainerType& Node::GetDofs()
{
    return mDofs;
}

/***********************************************************************************/
/***********************************************************************************/

const Node::DofsContainerType& Node::GetDofs() const
{
    return mDofs;
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

void Node::SortDofs()
{
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
