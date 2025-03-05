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
#include "spy_element.h"

namespace Kratos::Testing
{

void SpyElement::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    mSolutionStepInitialized = true;
}

void SpyElement::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    mSolutionStepFinalized = true;
}

void SpyElement::InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    mNonLinIterationInitialized = true;
}

void SpyElement::FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    mNonLinIterationFinalized = true;
}

bool SpyElement::IsSolutionStepInitialized() const { return mSolutionStepInitialized; }

bool SpyElement::IsSolutionStepFinalized() const { return mSolutionStepFinalized; }

bool SpyElement::IsNonLinIterationInitialized() const { return mNonLinIterationInitialized; }

bool SpyElement::IsNonLinIterationFinalized() const { return mNonLinIterationFinalized; }

} // namespace Kratos::Testing
