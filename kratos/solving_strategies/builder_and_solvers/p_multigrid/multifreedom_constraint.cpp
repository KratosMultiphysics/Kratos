//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//

// Project includes
#include "solving_strategies/builder_and_solvers/p_multigrid/multifreedom_constraint.hpp" // MultifreedomConstraint
#include "includes/variables.h" // CONSTRAINT_LABELS, CONSTRAINT_DOFS

// System includes
#include <algorithm> // std::copy, std::transform


namespace Kratos {


MultifreedomConstraint::MultifreedomConstraint(const IndexType Id,
                                               DofArray&& rDofs,
                                               const std::vector<std::size_t>& rConstraintLabels)
    : MasterSlaveConstraint(Id),
      mDofs(std::move(rDofs))
{
    Vector constraint_labels(rConstraintLabels.size());
    std::copy(rConstraintLabels.begin(), rConstraintLabels.end(), constraint_labels.begin());
    this->SetValue(CONSTRAINT_LABELS, constraint_labels);
}


void MultifreedomConstraint::InitializeSolutionStep(const ProcessInfo& rProcessInfo)
{
    this->SetConstraintDofs();
}


void MultifreedomConstraint::InitializeNonLinearIteration(const ProcessInfo& rProcessInfo)
{
    this->SetConstraintDofs();
}


void MultifreedomConstraint::GetDofList(DofPointerVectorType& rSlaveDofs,
                                        DofPointerVectorType& rMasterDofs,
                                        const ProcessInfo&) const
{
    rSlaveDofs.clear();
    rMasterDofs.resize(mDofs.size());
    std::copy(mDofs.begin(),
              mDofs.end(),
              rMasterDofs.begin());
}


void MultifreedomConstraint::SetDofList(const DofPointerVectorType& rSlaveDofs,
                                        const DofPointerVectorType& rMasterDofs,
                                        const ProcessInfo&)
{
    mDofs.resize(rSlaveDofs.size() + rMasterDofs.size());
    std::copy(rSlaveDofs.begin(),
              rSlaveDofs.end(),
              mDofs.begin());
    std::copy(rMasterDofs.begin(),
              rMasterDofs.end(),
              mDofs.begin() + rSlaveDofs.size());
}


void MultifreedomConstraint::EquationIdVector(EquationIdVectorType& rSlaveDofs,
                                              EquationIdVectorType& rMasterDofs,
                                              const ProcessInfo& rProcessInfo) const
{
    rSlaveDofs.clear();
    rMasterDofs.resize(mDofs.size());
    std::transform(mDofs.begin(),
                   mDofs.end(),
                   rMasterDofs.begin(),
                   [](const Dof<double>* p_dof) {return p_dof->EquationId();});
}


void MultifreedomConstraint::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, MasterSlaveConstraint);
    rSerializer.save("dofs", mDofs);
}


void MultifreedomConstraint::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, MasterSlaveConstraint);
    rSerializer.load("dofs", mDofs);
}


void MultifreedomConstraint::SetConstraintDofs()
{
    Vector constraint_dofs(mDofs.size());
    std::transform(mDofs.begin(),
                   mDofs.end(),
                   constraint_dofs.begin(),
                   [](Dof<double>* p_dof){return p_dof->EquationId();});
    this->SetValue(CONSTRAINT_DOFS, constraint_dofs);
}


} // namespace Kratos
