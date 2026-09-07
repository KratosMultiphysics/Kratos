// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//

#include "custom_processes/apply_incremental_c_phi_reduction_process.h"

#include "containers/model.h"
#include "custom_constitutive/mohr_coulomb_law.h"
#include "custom_constitutive/mohr_coulomb_with_tension_cutoff_elastoplastic_tangent_matrix.h"
#include "custom_utilities/check_utilities.hpp"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/process_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "includes/model_part.h"
#include "utilities/math_utils.h"

#include <cmath>

namespace Kratos
{
namespace
{

bool UsesUmatParametersForIncrementalReduction(const Properties& rProperties)
{
    KRATOS_ERROR_IF_NOT(rProperties.Has(CONSTITUTIVE_LAW))
        << "Properties do not have CONSTITUTIVE_LAW" << std::endl;

    const auto law_name = rProperties[CONSTITUTIVE_LAW]->Info();
    return law_name.find("UMAT") != std::string::npos || law_name.find("UDSM") != std::string::npos;
}

bool UsesInternalMohrCoulombModelForIncrementalReduction(const Element& rElement)
{
    KRATOS_ERROR_IF_NOT(rElement.GetProperties().Has(CONSTITUTIVE_LAW))
        << "Properties do not have CONSTITUTIVE_LAW" << std::endl;

    const auto* p_constitutive_law = rElement.GetProperties()[CONSTITUTIVE_LAW].get();
    return dynamic_cast<const MohrCoulombLaw*>(p_constitutive_law) != nullptr ||
           dynamic_cast<const MohrCoulombWithTensionCutOffElastoPlasticTangentMatrix*>(
               p_constitutive_law) != nullptr;
}

} // namespace

ApplyIncrementalCPhiReductionProcess::ApplyIncrementalCPhiReductionProcess(
    Model& rModel, const Parameters& rProcessSettings)
    : mrModelParts(ProcessUtilities::GetModelPartsFromSettings(
          rModel, rProcessSettings, ApplyIncrementalCPhiReductionProcess::Info())),
      mIncrementStrategy(rProcessSettings["increment_strategy"].GetString()),
      mTrialFactor(1.0),
      mLastConvergedFactor(std::nullopt),
      mFactorIncrement(rProcessSettings["factor_increment"].GetDouble()),
      mTrialNumber(0),
      mMaxTrials(0)
{
    const auto max_trials = rProcessSettings["max_trials"].GetInt();
    KRATOS_ERROR_IF(max_trials <= 0) << "max_trials must be greater than zero." << std::endl;
    mMaxTrials = static_cast<std::size_t>(max_trials);
}

void ApplyIncrementalCPhiReductionProcess::ExecuteInitializeSolutionStep()
{
    KRATOS_ERROR_IF(mTrialNumber >= mMaxTrials)
        << "The maximum number of c-phi reduction trials (" << mMaxTrials << ") has been reached."
        << std::endl;

    KRATOS_INFO("ApplyIncrementalCPhiReductionProcess")
        << "Trial " << mTrialNumber + 1 << " with safety factor F = " << mTrialFactor << std::endl;
    ApplyReduction();
}

void ApplyIncrementalCPhiReductionProcess::ExecuteFinalizeSolutionStep()
{
    mLastConvergedFactor = mTrialFactor;
    ++mTrialNumber;

    KRATOS_INFO("ApplyIncrementalCPhiReductionProcess")
        << "Converged with safety factor F = " << *mLastConvergedFactor << std::endl;

    if (mTrialNumber < mMaxTrials) {
        AdvanceFactor();
    }
}

void ApplyIncrementalCPhiReductionProcess::ExecuteFinalize()
{
    if (mLastConvergedFactor.has_value()) {
        KRATOS_INFO("ApplyIncrementalCPhiReductionProcess")
            << "Last converged safety factor F = " << *mLastConvergedFactor << std::endl;
    }
}

int ApplyIncrementalCPhiReductionProcess::Check()
{
    KRATOS_ERROR_IF(mIncrementStrategy != "fixed")
        << "Unsupported increment strategy: " << mIncrementStrategy
        << ". The only available strategy is 'fixed'." << std::endl;
    KRATOS_ERROR_IF_NOT(std::isfinite(mFactorIncrement) && mFactorIncrement > 0.0)
        << "factor_increment must be finite and greater than zero." << std::endl;
    KRATOS_ERROR_IF(mMaxTrials == 0) << "max_trials must be greater than zero." << std::endl;
    KRATOS_ERROR_IF(std::ranges::all_of(mrModelParts, [](const auto& r_model_part) {
        return r_model_part.get().Elements().empty();
    })) << "None of the provided model parts contains an element.\n";

    for (const auto& r_model_part : mrModelParts) {
        for (const auto& r_element : r_model_part.get().Elements()) {
            const auto&           r_properties = r_element.GetProperties();
            const CheckProperties check_properties(r_properties, "model part property",
                                                   CheckProperties::Bounds::AllInclusive);
            if (UsesUmatParametersForIncrementalReduction(r_properties)) {
                check_properties.CheckAvailability(UMAT_PARAMETERS);
                check_properties.Check(INDEX_OF_UMAT_PHI_PARAMETER, 1,
                                       static_cast<int>(r_properties[UMAT_PARAMETERS].size()));
                check_properties.Check(INDEX_OF_UMAT_C_PARAMETER, 1,
                                       static_cast<int>(r_properties[UMAT_PARAMETERS].size()));
            } else {
                check_properties.Check(GEO_COHESION);
                check_properties.Check(GEO_FRICTION_ANGLE);
            }
        }
    }
    return 0;
}

std::string ApplyIncrementalCPhiReductionProcess::Info() const
{
    return "ApplyIncrementalCPhiReductionProcess";
}

double ApplyIncrementalCPhiReductionProcess::GetTrialFactor() const { return mTrialFactor; }

bool ApplyIncrementalCPhiReductionProcess::HasLastConvergedFactor() const
{
    return mLastConvergedFactor.has_value();
}

double ApplyIncrementalCPhiReductionProcess::GetLastConvergedFactor() const
{
    KRATOS_ERROR_IF_NOT(mLastConvergedFactor.has_value())
        << "No c-phi reduction trial has converged yet." << std::endl;
    return *mLastConvergedFactor;
}

std::size_t ApplyIncrementalCPhiReductionProcess::GetTrialNumber() const
{
    return mTrialNumber;
}

bool ApplyIncrementalCPhiReductionProcess::IsFinished() const
{
    return mTrialNumber >= mMaxTrials;
}

double ApplyIncrementalCPhiReductionProcess::CalculateNewFactor() const
{
    KRATOS_ERROR_IF(mIncrementStrategy != "fixed")
        << "Increment strategy '" << mIncrementStrategy << "' is not implemented." << std::endl;
    return mTrialFactor + mFactorIncrement;
}

void ApplyIncrementalCPhiReductionProcess::AdvanceFactor()
{
    mTrialFactor = CalculateNewFactor();
}

void ApplyIncrementalCPhiReductionProcess::ApplyReduction()
{
    KRATOS_TRY
    for (auto& r_model_part : mrModelParts) {
        block_for_each(r_model_part.get().Elements(), [this, &r_model_part](Element& rElement) {
            const auto  property_id           = rElement.GetProperties().Id();
            const auto& r_original_properties = r_model_part.get().GetProperties(property_id);
            const auto  friction_angle = GetAndCheckFrictionAngle(r_original_properties, property_id);
            const auto  cohesion       = GetAndCheckCohesion(r_original_properties);

            SetCPhiAtElement(rElement,
                             CalculateReducedFrictionAngle(friction_angle, mTrialFactor),
                             CalculateReducedCohesion(cohesion, mTrialFactor));
        });
    }
    KRATOS_CATCH("")
}

double ApplyIncrementalCPhiReductionProcess::GetAndCheckCohesion(const Properties& rProperties)
{
    const auto cohesion = ConstitutiveLawUtilities::GetCohesion(rProperties);
    KRATOS_ERROR_IF(cohesion < 0.0) << "Cohesion must not be negative: " << cohesion << std::endl;
    return cohesion;
}

double ApplyIncrementalCPhiReductionProcess::GetAndCheckFrictionAngle(
    const Properties& rProperties, IndexType PropertyId)
{
    const auto friction_angle = ConstitutiveLawUtilities::GetFrictionAngleInDegrees(rProperties);
    KRATOS_ERROR_IF(friction_angle < 0.0 || friction_angle > 90.0)
        << "Friction angle in properties " << PropertyId << " is outside [0, 90] degrees: "
        << friction_angle << std::endl;
    return friction_angle;
}

double ApplyIncrementalCPhiReductionProcess::CalculateReducedCohesion(double Cohesion, double Factor)
{
    KRATOS_ERROR_IF_NOT(std::isfinite(Factor) && Factor >= 1.0)
        << "Safety factor F must be finite and greater than or equal to 1." << std::endl;
    return Cohesion / Factor;
}

double ApplyIncrementalCPhiReductionProcess::CalculateReducedFrictionAngle(double FrictionAngle,
                                                                           double Factor)
{
    KRATOS_ERROR_IF_NOT(std::isfinite(Factor) && Factor >= 1.0)
        << "Safety factor F must be finite and greater than or equal to 1." << std::endl;
    const auto tangent = std::tan(MathUtils<>::DegreesToRadians(FrictionAngle)) / Factor;
    return std::atan(tangent) * 180.0 / Globals::Pi;
}

void ApplyIncrementalCPhiReductionProcess::SetCPhiAtElement(Element& rElement,
                                                             double ReducedFrictionAngle,
                                                             double ReducedCohesion)
{
    const auto& r_properties         = rElement.GetProperties();
    auto        p_reduced_properties = Kratos::make_shared<Properties>(r_properties);

    if (UsesUmatParametersForIncrementalReduction(r_properties)) {
        auto& r_umat_parameters = p_reduced_properties->GetValue(UMAT_PARAMETERS);
        r_umat_parameters[r_properties[INDEX_OF_UMAT_PHI_PARAMETER] - 1] = ReducedFrictionAngle;
        r_umat_parameters[r_properties[INDEX_OF_UMAT_C_PARAMETER] - 1]   = ReducedCohesion;
    } else {
        p_reduced_properties->SetValue(GEO_FRICTION_ANGLE, ReducedFrictionAngle);
        p_reduced_properties->SetValue(GEO_COHESION, ReducedCohesion);
    }
    rElement.SetProperties(p_reduced_properties);

    if (UsesInternalMohrCoulombModelForIncrementalReduction(rElement)) {
        InitializeParametersForInternalMohrCoulombModel(rElement);
    }
}

void ApplyIncrementalCPhiReductionProcess::InitializeParametersForInternalMohrCoulombModel(
    Element& rElement)
{
    std::vector<ConstitutiveLaw::Pointer> constitutive_laws;
    const ProcessInfo                     dummy_process_info;
    rElement.CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, constitutive_laws, dummy_process_info);
    const auto& r_properties   = rElement.GetProperties();
    const auto  dummy_geometry = Geometry<Node>{};
    const auto  dummy_vector   = Vector();

    for (const auto& p_law : constitutive_laws) {
        if (const auto p_mohr_coulomb = dynamic_cast<MohrCoulombLaw*>(p_law.get())) {
            p_mohr_coulomb->InitializeMaterial(r_properties, dummy_geometry, dummy_vector);
        } else if (const auto p_mohr_coulomb =
                       dynamic_cast<MohrCoulombWithTensionCutOffElastoPlasticTangentMatrix*>(
                           p_law.get())) {
            p_mohr_coulomb->InitializeMaterial(r_properties, dummy_geometry, dummy_vector);
        }
    }
}

} // namespace Kratos
