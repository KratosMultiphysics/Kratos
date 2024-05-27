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

void SpyCondition::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    mSolutionStepInitialized = true;
}

void SpyCondition::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    mSolutionStepFinalized = true;
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
    mIsEquationIdRetrieved = true;
}

void SpyCondition::GetDofList(Condition::DofsVectorType& rElementalDofList,
                              const ProcessInfo& rCurrentProcessInfo) const
{
    mIsGetDofListCalled = true;
}

bool SpyCondition::IsSolutionStepInitialized() const
{
    return mSolutionStepInitialized;
}

bool SpyCondition::IsSolutionStepFinalized() const
{
    return mSolutionStepFinalized;
}

bool SpyCondition::IsNonLinIterationInitialized() const
{
    return mNonLinIterationInitialized;
}

bool SpyCondition::IsNonLinIterationFinalized() const
{
    return mNonLinIterationFinalized;
}

bool SpyCondition::IsEquationIdRetrieved() const
{
    return mIsEquationIdRetrieved;
}

bool SpyCondition::IsGetDofListCalled() const
{
    return mIsGetDofListCalled;
}

} // namespace Kratos::Testing
