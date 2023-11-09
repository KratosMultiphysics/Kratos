// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//
#include "spy_condition.h"

namespace Kratos::Testing {

void SpyCondition::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    mSolutionStepFinalized = true;
}

void SpyCondition::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    mSolutionStepInitialized = true;
}

void SpyCondition::InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    mNonLinIterationInitialized = true;
}

void SpyCondition::FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    mNonLinIterationFinalized = true;
}

void SpyCondition::EquationIdVector(Condition::EquationIdVectorType& rResult,
                                    const ProcessInfo& rCurrentProcessInfo) const
{
    mIsEquationIdSet = true;
}

void SpyCondition::GetDofList(Condition::DofsVectorType& rElementalDofList,
                              const ProcessInfo& rCurrentProcessInfo) const
{
    mIsGetDofListCalled = true;
}

bool SpyCondition::IsSolutionStepFinalized() const
{
    return mSolutionStepFinalized;
}

bool SpyCondition::IsSolutionStepInitialized() const
{
    return mSolutionStepInitialized;
}

bool SpyCondition::IsNonLinIterationInitialized() const
{
    return mNonLinIterationInitialized;
}

bool SpyCondition::IsNonLinIterationFinalized() const
{
    return mNonLinIterationFinalized;
}

bool SpyCondition::IsEquationIdSet() const
{
    return mIsEquationIdSet;
}

bool SpyCondition::IsGetDofListCalled() const
{
    return mIsGetDofListCalled;
}

} // namespace Kratos::Testing
