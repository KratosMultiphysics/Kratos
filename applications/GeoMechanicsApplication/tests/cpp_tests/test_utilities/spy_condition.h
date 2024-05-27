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

#pragma once

#include "includes/condition.h"

namespace Kratos::Testing {

class SpyCondition : public Condition {
public:
    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;
    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;
    void InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) override;
    void FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(EquationIdVectorType& rResult,
                          const ProcessInfo& rCurrentProcessInfo) const override;

    void GetDofList(DofsVectorType& rElementalDofList,
                    const ProcessInfo& rCurrentProcessInfo) const override;

    bool IsSolutionStepInitialized() const;
    bool IsSolutionStepFinalized() const;
    bool IsNonLinIterationInitialized() const;
    bool IsNonLinIterationFinalized() const;
    bool IsEquationIdRetrieved() const;
    bool IsGetDofListCalled() const;

private:
    bool mSolutionStepInitialized = false;
    bool mSolutionStepFinalized = false;
    bool mNonLinIterationInitialized = false;
    bool mNonLinIterationFinalized = false;
    mutable bool mIsEquationIdRetrieved = false;
    mutable bool mIsGetDofListCalled = false;
};

} // namespace Kratos
