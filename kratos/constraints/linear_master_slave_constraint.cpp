//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aditya Ghantasala
//

// System includes

// External includes

// Project includes
#include "linear_master_slave_constraint.h"
#include "utilities/atomic_utilities.h"
#include "includes/kratos_flags.h"

namespace Kratos
{

LinearMasterSlaveConstraint::LinearMasterSlaveConstraint(
    IndexType Id,
    NodeType& rMasterNode,
    const VariableType& rMasterVariable,
    NodeType& rSlaveNode,
    const VariableType& rSlaveVariable,
    const double Weight,
    const double Constant
    ) : MasterSlaveConstraint(Id)
{
    // Resizing the memeber variables
    mRelationMatrix.resize(1,1,false);
    mConstantVector.resize(1,false);

    // Obtaining the dofs from the variables
    mSlaveDofsVector.push_back(rSlaveNode.pGetDof(rSlaveVariable));
    mMasterDofsVector.push_back(rMasterNode.pGetDof(rMasterVariable));

    mRelationMatrix(0,0) = Weight;
    mConstantVector(0) = Constant;

    // Setting the slave flag on the node
    rSlaveNode.Set(SLAVE);
}

void LinearMasterSlaveConstraint::EquationIdVector(
    EquationIdVectorType& rSlaveEquationIds,
    EquationIdVectorType& rMasterEquationIds,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    if (rSlaveEquationIds.size() != mSlaveDofsVector.size())
        rSlaveEquationIds.resize(mSlaveDofsVector.size());

    if (rMasterEquationIds.size() != mMasterDofsVector.size())
        rMasterEquationIds.resize(mMasterDofsVector.size());

    for(IndexType i=0; i<rSlaveEquationIds.size(); ++i)
        rSlaveEquationIds[i] = mSlaveDofsVector[i]->EquationId();

    for(IndexType i=0; i<rMasterEquationIds.size(); ++i)
        rMasterEquationIds[i] = mMasterDofsVector[i]->EquationId();
}

void LinearMasterSlaveConstraint::ResetSlaveDofs(const ProcessInfo& rCurrentProcessInfo)
{
    for (IndexType i = 0; i < mSlaveDofsVector.size(); ++i) {
        AtomicMult(mSlaveDofsVector[i]->GetSolutionStepValue(), 0.0);
    }
}

void LinearMasterSlaveConstraint::Apply(const ProcessInfo& rCurrentProcessInfo)
{
    // Saving the master dofs values
    Vector master_dofs_values(mMasterDofsVector.size());

    for (IndexType i = 0; i < mMasterDofsVector.size(); ++i) {
        master_dofs_values[i] = mMasterDofsVector[i]->GetSolutionStepValue();
    }

    // Apply the constraint to the slave dofs
    for (IndexType i = 0; i < mRelationMatrix.size1(); ++i) {
        double aux = mConstantVector[i];
        for(IndexType j = 0; j < mRelationMatrix.size2(); ++j) {
            aux += mRelationMatrix(i,j) * master_dofs_values[j];
        }

        AtomicAdd(mSlaveDofsVector[i]->GetSolutionStepValue(), aux);
    }
}

void LinearMasterSlaveConstraint::SetLocalSystem(
    const MatrixType& rRelationMatrix,
    const VectorType& rConstantVector,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (mRelationMatrix.size1() != rRelationMatrix.size1() || mRelationMatrix.size2() != rRelationMatrix.size2())
        mRelationMatrix.resize(rRelationMatrix.size1(), rRelationMatrix.size2(), false);

    noalias(mRelationMatrix) = rRelationMatrix;

    if (mConstantVector.size() != rConstantVector.size())
        mConstantVector.resize(rConstantVector.size(), false);
    noalias(mConstantVector) = rConstantVector;
}

void LinearMasterSlaveConstraint::CalculateLocalSystem(
    MatrixType& rRelationMatrix,
    VectorType& rConstantVector,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    if (rRelationMatrix.size1() != mRelationMatrix.size1() || rRelationMatrix.size2() != mRelationMatrix.size2())
        rRelationMatrix.resize(mRelationMatrix.size1(), mRelationMatrix.size2(), false);
    noalias(rRelationMatrix) = mRelationMatrix;

    if (rConstantVector.size() != mConstantVector.size())
        rConstantVector.resize(mConstantVector.size(), false);
    noalias(rConstantVector) = mConstantVector;
}

void LinearMasterSlaveConstraint::save(Serializer &rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, MasterSlaveConstraint);
    rSerializer.save("SlaveDofVec", mSlaveDofsVector);
    rSerializer.save("MasterDofVec", mMasterDofsVector);
    rSerializer.save("RelationMat", mRelationMatrix);
    rSerializer.save("ConstantVec", mConstantVector);
}

void LinearMasterSlaveConstraint::load(Serializer &rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, MasterSlaveConstraint);
    rSerializer.load("SlaveDofVec", mSlaveDofsVector);
    rSerializer.load("MasterDofVec", mMasterDofsVector);
    rSerializer.load("RelationMat", mRelationMatrix);
    rSerializer.load("ConstantVec", mConstantVector);
}

///@}

} // namespace Kratos
