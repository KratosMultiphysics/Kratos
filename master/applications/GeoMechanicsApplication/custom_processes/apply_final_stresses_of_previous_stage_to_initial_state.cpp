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

#include "apply_final_stresses_of_previous_stage_to_initial_state.h"
#include "includes/initial_state.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

#include <vector>

namespace Kratos
{

ApplyFinalStressesOfPreviousStageToInitialState::ApplyFinalStressesOfPreviousStageToInitialState(ModelPart& rModelPart,
                                                                                                 const Parameters&)
    : mrModelPart(rModelPart)
{
}

void ApplyFinalStressesOfPreviousStageToInitialState::ExecuteInitialize()
{
    block_for_each(mrModelPart.Elements(), [this](Element& rElement) {
        std::vector<Vector> stresses_on_integration_points;
        rElement.CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, stresses_on_integration_points,
                                              mrModelPart.GetProcessInfo());
        if (stresses_on_integration_points.empty()) {
            rElement.CalculateOnIntegrationPoints(
                CAUCHY_STRESS_VECTOR, stresses_on_integration_points, mrModelPart.GetProcessInfo());
        }
        std::vector<ConstitutiveLaw::Pointer> constitutive_laws;
        rElement.CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, constitutive_laws, mrModelPart.GetProcessInfo());

        CheckRetrievedElementData(constitutive_laws, stresses_on_integration_points, rElement.GetId());
        mStressesByElementId[rElement.GetId()] = stresses_on_integration_points;
    });
}

void ApplyFinalStressesOfPreviousStageToInitialState::ExecuteBeforeSolutionLoop()
{
    block_for_each(mrModelPart.Elements(), [this](Element& rElement) {
        std::vector<ConstitutiveLaw::Pointer> constitutive_laws;
        rElement.CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, constitutive_laws, mrModelPart.GetProcessInfo());
        const auto stresses_on_integration_points = mStressesByElementId.at(rElement.GetId());
        for (auto i = std::size_t{0}; i < constitutive_laws.size(); ++i) {
            auto p_initial_state = make_intrusive<InitialState>();
            p_initial_state->SetInitialStressVector(stresses_on_integration_points[i]);
            p_initial_state->SetInitialStrainVector(ZeroVector{constitutive_laws[i]->GetStrainSize()});
            constitutive_laws[i]->SetInitialState(p_initial_state);
            constitutive_laws[i]->InitializeMaterial(rElement.GetProperties(), rElement.GetGeometry(), {});
        }
    });
    mStressesByElementId.clear();
}

void ApplyFinalStressesOfPreviousStageToInitialState::CheckRetrievedElementData(
    const std::vector<ConstitutiveLaw::Pointer>& rConstitutiveLaws,
    const std::vector<Vector>&                   rStressesOnIntegrationPoints,
    IndexType                                    ElementId)
{
    KRATOS_ERROR_IF(rConstitutiveLaws.empty())
        << "The constitutive laws on the integration points could not be retrieved for element "
        << ElementId << std::endl;
    KRATOS_ERROR_IF(rStressesOnIntegrationPoints.empty())
        << "The stress vectors on the integration points could not be retrieved for element "
        << ElementId << std::endl;
    KRATOS_ERROR_IF(rStressesOnIntegrationPoints.size() != rConstitutiveLaws.size())
        << "Number of retrieved stress vectors (" << rStressesOnIntegrationPoints.size()
        << ") does not match the number of constitutive laws (" << rConstitutiveLaws.size()
        << ") for element " << ElementId << std::endl;
}

} // namespace Kratos