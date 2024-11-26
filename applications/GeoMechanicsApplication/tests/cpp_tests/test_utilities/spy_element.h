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

#include "includes/element.h"
#include <gmock/gmock.h>

namespace Kratos::Testing
{

class SpyElement : public Element
{
public:
    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;
    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;
    void InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) override;
    void FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) override;

    MOCK_METHOD(void, EquationIdVector, (EquationIdVectorType&, const ProcessInfo&), (const, override));
    MOCK_METHOD(void, GetDofList, (DofsVectorType&, const ProcessInfo&), (const, override));

    bool IsSolutionStepInitialized() const;
    bool IsSolutionStepFinalized() const;
    bool IsNonLinIterationInitialized() const;
    bool IsNonLinIterationFinalized() const;

private:
    bool mSolutionStepInitialized    = false;
    bool mSolutionStepFinalized      = false;
    bool mNonLinIterationInitialized = false;
    bool mNonLinIterationFinalized   = false;
};

} // namespace Kratos::Testing
