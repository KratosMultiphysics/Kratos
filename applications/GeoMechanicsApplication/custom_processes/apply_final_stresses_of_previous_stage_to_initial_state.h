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

#include "includes/constitutive_law.h"
#include "processes/process.h"

namespace Kratos
{

class ModelPart;
class Parameters;
class Element;

class KRATOS_API(GEO_MECHANICS_APPLICATION) ApplyFinalStressesOfPreviousStageToInitialState : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ApplyFinalStressesOfPreviousStageToInitialState);
    ApplyFinalStressesOfPreviousStageToInitialState(ModelPart& rModelPart, const Parameters&);
    ~ApplyFinalStressesOfPreviousStageToInitialState() override = default;

    ApplyFinalStressesOfPreviousStageToInitialState(const ApplyFinalStressesOfPreviousStageToInitialState&) = delete;
    ApplyFinalStressesOfPreviousStageToInitialState& operator=(const ApplyFinalStressesOfPreviousStageToInitialState&) = delete;

    void ExecuteInitialize() override;
    void ExecuteBeforeSolutionLoop() override;

private:
    ModelPart& mrModelPart;

    static void CheckRetrievedElementData(const std::vector<ConstitutiveLaw::Pointer>& rConstitutiveLaws,
                                          const std::vector<Vector>& rStressesOnIntegrationPoints,
                                          IndexType                  ElementId);

    std::map<IndexType, std::vector<Vector>> mStressesByElementId;
};

} // namespace Kratos