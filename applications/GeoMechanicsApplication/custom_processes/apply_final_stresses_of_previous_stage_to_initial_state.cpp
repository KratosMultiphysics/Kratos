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
#include "custom_constitutive/reference_3d_mohr_coulomb_plane_strain_law.h"
#include "custom_utilities/process_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "includes/initial_state.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

#include <algorithm>
#include <vector>

namespace Kratos
{
using namespace std::string_literals;

namespace
{
bool ContainsOnlyNonEmptyVectors(const std::vector<Vector>& rVectors)
{
    return !rVectors.empty() &&
           std::ranges::all_of(rVectors, [](const Vector& rVector) { return !rVector.empty(); });
}
} // namespace

ApplyFinalStressesOfPreviousStageToInitialState::ApplyFinalStressesOfPreviousStageToInitialState(
    Model& rModel, const Parameters& rProcessSettings)
{
    mrModelParts = ProcessUtilities::GetModelPartsFromSettings(
        rModel, rProcessSettings, ApplyFinalStressesOfPreviousStageToInitialState::Info());
}

void ApplyFinalStressesOfPreviousStageToInitialState::ExecuteInitialize()
{
    for (const auto& r_model_part : mrModelParts) {
        for (auto& r_element : r_model_part.get().Elements()) {
            std::vector<Vector> stresses_on_integration_points;
            auto                retrieved_stress_variable = "PK2_STRESS_VECTOR"s;
            r_element.CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, stresses_on_integration_points,
                                                   r_model_part.get().GetProcessInfo());
            if (!ContainsOnlyNonEmptyVectors(stresses_on_integration_points)) {
                retrieved_stress_variable = "CAUCHY_STRESS_VECTOR"s;
                r_element.CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, stresses_on_integration_points,
                                                       r_model_part.get().GetProcessInfo());
            }
            if (!ContainsOnlyNonEmptyVectors(stresses_on_integration_points)) {
                retrieved_stress_variable = "GEO_EFFECTIVE_TRACTION_VECTOR"s;
                r_element.CalculateOnIntegrationPoints(
                    GEO_EFFECTIVE_TRACTION_VECTOR, stresses_on_integration_points,
                    r_model_part.get().GetProcessInfo());
            }
            std::vector<ConstitutiveLaw::Pointer> constitutive_laws;
            r_element.CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, constitutive_laws,
                                                   r_model_part.get().GetProcessInfo());

            CheckRetrievedElementData(constitutive_laws, stresses_on_integration_points, r_element.GetId());
            if (r_element.GetId() == 4) {
                KRATOS_WATCH("Element 4: stresses captured at the start of stage 2")
                KRATOS_WATCH(retrieved_stress_variable)
                for (IndexType integration_point_index = 0;
                     integration_point_index < stresses_on_integration_points.size();
                     ++integration_point_index) {
                    KRATOS_WATCH(integration_point_index)
                    KRATOS_WATCH(stresses_on_integration_points[integration_point_index])
                }
            }
            mStressesByElementId[r_element.GetId()] = stresses_on_integration_points;
        }
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
                if (rElement.GetId() == 4) {
                    if (auto p_reference_law =
                            dynamic_cast<Reference3DMohrCoulombPlaneStrainLaw*>(
                                constitutive_laws[i].get())) {
                        p_reference_law->SetInitialStateTracing(true);
                    }
                    KRATOS_WATCH("Element 4: InitialState assigned to the stage 2 wrapper")
                    KRATOS_WATCH(i)
                    KRATOS_WATCH(p_initial_state->GetInitialStressVector())
                    KRATOS_WATCH(p_initial_state->GetInitialStrainVector())
                }
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
    KRATOS_ERROR_IF(std::ranges::any_of(rStressesOnIntegrationPoints,
                                       [](const Vector& rVector) { return rVector.empty(); }))
        << "One or more stress vectors retrieved on the integration points are empty for element "
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
