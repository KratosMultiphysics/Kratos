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
#include "containers/model.h"
#include "custom_utilities/process_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "includes/initial_state.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

#include <string>
#include <vector>

namespace Kratos
{
using namespace std::string_literals;

ApplyFinalStressesOfPreviousStageToInitialState::ApplyFinalStressesOfPreviousStageToInitialState(
    Model& rModel, const Parameters& rProcessSettings)
{
    mrModelParts = ProcessUtilities::GetModelPartsFromSettings(
        rModel, rProcessSettings, ApplyFinalStressesOfPreviousStageToInitialState::Info());
}

void ApplyFinalStressesOfPreviousStageToInitialState::ExecuteInitialize()
{
    for (const auto& r_model_part : mrModelParts) {
        block_for_each(r_model_part.get().Elements(), [&r_model_part, this](Element& rElement) {
            std::vector<Vector> stresses_on_integration_points;
            rElement.CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, stresses_on_integration_points,
                                                  r_model_part.get().GetProcessInfo());
            if (stresses_on_integration_points.empty()) {
                rElement.CalculateOnIntegrationPoints(GEO_EFFECTIVE_TRACTION_VECTOR, stresses_on_integration_points,
                                                      r_model_part.get().GetProcessInfo());
            }
            std::vector<ConstitutiveLaw::Pointer> constitutive_laws;
            rElement.CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, constitutive_laws,
                                                  r_model_part.get().GetProcessInfo());

            CheckRetrievedElementData(constitutive_laws, stresses_on_integration_points, rElement.GetId());
            mStressesByElementId[rElement.GetId()] = stresses_on_integration_points;
        });
    }
}

void ApplyFinalStressesOfPreviousStageToInitialState::ExecuteBeforeSolutionLoop()
{
    for (const auto& r_model_part : mrModelParts) {
        block_for_each(r_model_part.get().Elements(), [&r_model_part, this](Element& rElement) {
            std::vector<ConstitutiveLaw::Pointer> constitutive_laws;
            rElement.CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, constitutive_laws,
                                                  r_model_part.get().GetProcessInfo());
            const auto stresses_on_integration_points = mStressesByElementId.at(rElement.GetId());
            for (auto i = std::size_t{0}; i < constitutive_laws.size(); ++i) {
                auto p_initial_state = make_intrusive<InitialState>();
                p_initial_state->SetInitialStressVector(stresses_on_integration_points[i]);
                p_initial_state->SetInitialStrainVector(ZeroVector{constitutive_laws[i]->GetStrainSize()});
                constitutive_laws[i]->SetInitialState(p_initial_state);
                constitutive_laws[i]->InitializeMaterial(rElement.GetProperties(), rElement.GetGeometry(), {});
            }
        });
    }
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

std::string ApplyFinalStressesOfPreviousStageToInitialState::Info() const
{
    return "ApplyFinalStressesOfPreviousStageToInitialState"s;
}

} // namespace Kratos