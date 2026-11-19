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
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

#include <algorithm>
#include <cmath>
#include <functional>

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

Parameters GetAdaptiveIncrementSettings(const Parameters& rProcessSettings)
{
    auto settings = rProcessSettings.Has("adaptive_increment")
                        ? rProcessSettings["adaptive_increment"].Clone()
                        : Parameters("{}");
    settings.ValidateAndAssignDefaults(Parameters(R"(
    {
        "displacement_measure": "max_norm",
        "slope_averaging_steps": 2,
        "min_increment": 0.001,
        "max_increment": 0.2,
        "max_growth_factor": 1.5,
        "target_factor": 0.0,
        "target_approach_fraction": 0.25
    })"));
    return settings;
}

Parameters GetRetrySettings(const Parameters& rProcessSettings)
{
    auto settings = rProcessSettings.Has("retry") ? rProcessSettings["retry"].Clone() : Parameters("{}");
    settings.ValidateAndAssignDefaults(Parameters(R"(
    {
        "max_cuts": 0,
        "cut_factor": 0.5,
        "bracket_tolerance": 0.0,
        "bracket_max_increment": 0.05,
        "bracket_max_stiffness_ratio": 0.3
    })"));
    return settings;
}

Parameters GetFailureCriterionSettings(const Parameters& rProcessSettings)
{
    auto settings = rProcessSettings.Has("failure_criterion")
                        ? rProcessSettings["failure_criterion"].Clone()
                        : Parameters("{}");
    settings.ValidateAndAssignDefaults(Parameters(R"(
    {
        "type": "non_convergence",
        "stop_stiffness_ratio": 0.02,
        "extrapolation_points": 3
    })"));
    return settings;
}

std::vector<double> GetFactorList(const Parameters& rProcessSettings)
{
    std::vector<double> result;
    if (rProcessSettings.Has("factor_list")) {
        const auto factors = rProcessSettings["factor_list"];
        for (std::size_t i = 0; i < factors.size(); ++i) {
            result.push_back(factors[i].GetDouble());
        }
    }
    return result;
}

std::vector<std::size_t> GetMonitoredNodeIds(const Parameters& rProcessSettings)
{
    std::vector<std::size_t> result;
    if (rProcessSettings.Has("monitored_node_ids")) {
        const auto node_ids = rProcessSettings["monitored_node_ids"];
        for (std::size_t i = 0; i < node_ids.size(); ++i) {
            result.push_back(static_cast<std::size_t>(node_ids[i].GetInt()));
        }
    }
    return result;
}

} // namespace

ApplyIncrementalCPhiReductionProcess::ApplyIncrementalCPhiReductionProcess(
    Model& rModel, const Parameters& rProcessSettings)
    : mrModelParts(ProcessUtilities::GetModelPartsFromSettings(
          rModel, rProcessSettings, ApplyIncrementalCPhiReductionProcess::Info())),
      mIncrementStrategy(rProcessSettings["increment_strategy"].GetString()),
      mTrialFactor(1.0),
      mLastConvergedFactor(std::nullopt),
      mInitialFactorIncrement(rProcessSettings["factor_increment"].GetDouble()),
      mFactorIncrement(mInitialFactorIncrement),
      mFactorList(GetFactorList(rProcessSettings)),
      mTrialNumber(0),
      mMaxTrials(0),
      mMonitoredNodeIds(GetMonitoredNodeIds(rProcessSettings))
{
    if (mIncrementStrategy == "list" && !mFactorList.empty()) mTrialFactor = mFactorList.front();

    const auto max_trials = rProcessSettings["max_trials"].GetInt();
    KRATOS_ERROR_IF(max_trials <= 0) << "max_trials must be greater than zero." << std::endl;
    mMaxTrials = static_cast<std::size_t>(max_trials);

    const auto adaptive_settings = GetAdaptiveIncrementSettings(rProcessSettings);
    mDisplacementMeasure         = adaptive_settings["displacement_measure"].GetString();
    const auto averaging_steps   = adaptive_settings["slope_averaging_steps"].GetInt();
    KRATOS_ERROR_IF(averaging_steps <= 0)
        << "adaptive_increment.slope_averaging_steps must be greater than zero." << std::endl;
    mSlopeAveragingSteps = static_cast<std::size_t>(averaging_steps);
    mMinFactorIncrement  = adaptive_settings["min_increment"].GetDouble();
    mMaxFactorIncrement  = adaptive_settings["max_increment"].GetDouble();
    mMaxGrowthFactor     = adaptive_settings["max_growth_factor"].GetDouble();
    mTargetFactor            = adaptive_settings["target_factor"].GetDouble();
    mTargetApproachFraction  = adaptive_settings["target_approach_fraction"].GetDouble();
    KRATOS_ERROR_IF(mTargetFactor > 0.0 && (mTargetApproachFraction <= 0.0 || mTargetApproachFraction > 1.0))
        << "adaptive_increment.target_approach_fraction must be in (0, 1]." << std::endl;

    const auto retry_settings = GetRetrySettings(rProcessSettings);
    const auto max_cuts       = retry_settings["max_cuts"].GetInt();
    KRATOS_ERROR_IF(max_cuts < 0) << "retry.max_cuts must not be negative." << std::endl;
    mMaxCuts   = static_cast<std::size_t>(max_cuts);
    mCutFactor = retry_settings["cut_factor"].GetDouble();
    mBracketTolerance    = retry_settings["bracket_tolerance"].GetDouble();
    mBracketMaxIncrement = retry_settings["bracket_max_increment"].GetDouble();
    mBracketMaxStiffnessRatio = retry_settings["bracket_max_stiffness_ratio"].GetDouble();

    const auto failure_settings = GetFailureCriterionSettings(rProcessSettings);
    mFailureCriterion           = failure_settings["type"].GetString();
    mStopStiffnessRatio         = failure_settings["stop_stiffness_ratio"].GetDouble();
    const auto points           = failure_settings["extrapolation_points"].GetInt();
    KRATOS_ERROR_IF(points < 2) << "failure_criterion.extrapolation_points must be at least 2." << std::endl;
    mExtrapolationPoints = static_cast<std::size_t>(points);
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
    mCutsInCurrentStep   = 0;
    ++mTrialNumber;

    KRATOS_INFO("ApplyIncrementalCPhiReductionProcess")
        << "Converged with safety factor F = " << *mLastConvergedFactor << std::endl;

    UpdateSlope();

    if (IsBracketClosed()) {
        mStoppedByBracket = true;
        KRATOS_INFO("ApplyIncrementalCPhiReductionProcess")
            << "Factor of safety bracketed in [" << *mLastConvergedFactor << ", " << *mUpperBound
            << "] (tolerance " << mBracketTolerance << "); stopping the c-phi reduction" << std::endl;
        return;
    }

    if (mTrialNumber < mMaxTrials && !mStoppedByStiffness) {
        AdvanceFactor();
    }
}

void ApplyIncrementalCPhiReductionProcess::ExecuteFinalize()
{
    if (mLastConvergedFactor.has_value()) {
        KRATOS_INFO("ApplyIncrementalCPhiReductionProcess")
            << "Last converged safety factor F = " << *mLastConvergedFactor << std::endl;
        KRATOS_INFO("ApplyIncrementalCPhiReductionProcess")
            << "Stop reason: " << GetStopReason() << ", stiffness ratio k/k0 = " << GetStiffnessRatio()
            << ", factor of safety F = " << GetFactorOfSafety() << std::endl;
    }
}

int ApplyIncrementalCPhiReductionProcess::Check()
{
    KRATOS_ERROR_IF(mIncrementStrategy != "fixed" && mIncrementStrategy != "adaptive" &&
                    mIncrementStrategy != "list")
        << "Unsupported increment strategy: " << mIncrementStrategy
        << ". The available strategies are 'fixed', 'adaptive' and 'list'." << std::endl;
    if (mIncrementStrategy == "list") {
        KRATOS_ERROR_IF(mFactorList.empty()) << "The 'list' strategy requires factor_list." << std::endl;
        KRATOS_ERROR_IF(mFactorList.front() < 1.0) << "factor_list values must be at least 1." << std::endl;
        KRATOS_ERROR_IF(std::ranges::adjacent_find(mFactorList, std::greater_equal<>{}) != mFactorList.end())
            << "factor_list must be strictly increasing." << std::endl;
        KRATOS_ERROR_IF(mMaxCuts > 0) << "The 'list' strategy does not support retry.max_cuts." << std::endl;
    }
    KRATOS_ERROR_IF_NOT(std::isfinite(mInitialFactorIncrement) && mInitialFactorIncrement > 0.0)
        << "factor_increment must be finite and greater than zero." << std::endl;
    KRATOS_ERROR_IF(mMaxTrials == 0) << "max_trials must be greater than zero." << std::endl;
    KRATOS_ERROR_IF(mFailureCriterion != "non_convergence" && mFailureCriterion != "stiffness_ratio")
        << "Unsupported failure_criterion.type: " << mFailureCriterion
        << ". The available types are 'non_convergence' and 'stiffness_ratio'." << std::endl;
    KRATOS_ERROR_IF_NOT(mStopStiffnessRatio > 0.0 && mStopStiffnessRatio < 1.0)
        << "failure_criterion.stop_stiffness_ratio must be in the open interval (0, 1)." << std::endl;
    KRATOS_ERROR_IF_NOT(mCutFactor > 0.0 && mCutFactor < 1.0)
        << "retry.cut_factor must be in the open interval (0, 1)." << std::endl;
    if (mIncrementStrategy == "adaptive") {
        KRATOS_ERROR_IF(mDisplacementMeasure != "max_norm" && mDisplacementMeasure != "monitored_nodes")
            << "Unsupported adaptive_increment.displacement_measure: " << mDisplacementMeasure
            << ". The available measures are 'max_norm' and 'monitored_nodes'." << std::endl;
        KRATOS_ERROR_IF_NOT(mMinFactorIncrement > 0.0 && mMinFactorIncrement <= mMaxFactorIncrement)
            << "adaptive_increment requires 0 < min_increment <= max_increment." << std::endl;
        KRATOS_ERROR_IF(mMaxGrowthFactor < 1.0)
            << "adaptive_increment.max_growth_factor must be at least 1." << std::endl;
        if (mDisplacementMeasure == "monitored_nodes") {
            KRATOS_ERROR_IF(mMonitoredNodeIds.empty())
                << "The 'monitored_nodes' displacement measure requires monitored_node_ids." << std::endl;
            const auto& r_root_model_part = mrModelParts.front().get().GetRootModelPart();
            for (const auto node_id : mMonitoredNodeIds) {
                KRATOS_ERROR_IF_NOT(r_root_model_part.HasNode(node_id))
                    << "Monitored node " << node_id << " does not exist." << std::endl;
            }
        }
    }
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

double ApplyIncrementalCPhiReductionProcess::GetFactorIncrement() const { return mFactorIncrement; }

std::size_t ApplyIncrementalCPhiReductionProcess::GetTotalNumberOfCuts() const
{
    return mTotalNumberOfCuts;
}

std::string ApplyIncrementalCPhiReductionProcess::GetStopReason() const
{
    if (mStoppedByStiffness) return "stiffness_ratio";
    if (mStoppedByBracket) return "bracket";
    if (mStoppedByNonConvergence) return "non_convergence";
    if (IsListCompleted()) return "list_completed";
    if (mTrialNumber >= mMaxTrials) return "max_trials";
    return "running";
}

double ApplyIncrementalCPhiReductionProcess::GetStiffnessRatio() const
{
    return mFactorsAndStiffnessRatios.empty() ? 1.0 : mFactorsAndStiffnessRatios.back().second;
}

double ApplyIncrementalCPhiReductionProcess::GetFactorOfSafety() const
{
    const auto last_converged_factor = GetLastConvergedFactor();
    if (mFailureCriterion == "non_convergence") return last_converged_factor;

    // With the stiffness criterion the factor of safety never comes from non-convergence; when the
    // reduction stops before reaching the criterion, the extrapolation of k -> 0 is still the best estimate
    const auto extrapolated_factor =
        mExtrapolatedFactor.has_value() ? mExtrapolatedFactor : ExtrapolateFactorAtZeroStiffness();
    return extrapolated_factor.value_or(last_converged_factor);
}

bool ApplyIncrementalCPhiReductionProcess::CutIncrementAfterNonConvergence()
{
    const auto is_small_increment = mLastConvergedFactor.has_value() &&
                                    mTrialFactor - *mLastConvergedFactor <= mBracketMaxIncrement + 1e-12;
    const auto is_softening = GetStiffnessRatio() < mBracketMaxStiffnessRatio;
    if (mBracketTolerance > 0.0 && is_small_increment && !is_softening) {
        KRATOS_INFO("ApplyIncrementalCPhiReductionProcess")
            << "Failed trial at F = " << mTrialFactor << " not taken as upper bound: k/k0 = "
            << GetStiffnessRatio() << " >= " << mBracketMaxStiffnessRatio << " (slope still stiff)" << std::endl;
    }
    if (mBracketTolerance > 0.0 && is_small_increment && is_softening) {
        mUpperBound = mUpperBound.has_value() ? std::min(*mUpperBound, mTrialFactor) : mTrialFactor;
        KRATOS_INFO("ApplyIncrementalCPhiReductionProcess")
            << "Upper bound of the factor of safety: F < " << *mUpperBound << std::endl;
        if (IsBracketClosed()) {
            mStoppedByBracket = true;
            KRATOS_INFO("ApplyIncrementalCPhiReductionProcess")
                << "Factor of safety bracketed in [" << *mLastConvergedFactor << ", " << *mUpperBound
                << "] (tolerance " << mBracketTolerance << "); stopping the c-phi reduction" << std::endl;
            return false;
        }
    }

    if (!mLastConvergedFactor.has_value() || mCutsInCurrentStep >= mMaxCuts) {
        mStoppedByNonConvergence = true;
        KRATOS_INFO("ApplyIncrementalCPhiReductionProcess")
            << "No convergence at F = " << mTrialFactor << " after " << mCutsInCurrentStep
            << " cuts of the factor increment; stopping the c-phi reduction" << std::endl;
        return false;
    }

    ++mCutsInCurrentStep;
    ++mTotalNumberOfCuts;
    // cut the increment that was actually tried (it may have been limited by the target or the bracket)
    mFactorIncrement = (mTrialFactor - *mLastConvergedFactor) * mCutFactor;
    mTrialFactor = *mLastConvergedFactor + LimitIncrement(*mLastConvergedFactor, mFactorIncrement);
    KRATOS_INFO("ApplyIncrementalCPhiReductionProcess")
        << "No convergence; cut " << mCutsInCurrentStep << " of " << mMaxCuts
        << ": retrying from F = " << *mLastConvergedFactor << " with factor increment "
        << mFactorIncrement << std::endl;
    return true;
}

std::size_t ApplyIncrementalCPhiReductionProcess::GetTrialNumber() const
{
    return mTrialNumber;
}

bool ApplyIncrementalCPhiReductionProcess::IsFinished() const
{
    return mTrialNumber >= mMaxTrials || mStoppedByNonConvergence || mStoppedByStiffness ||
           mStoppedByBracket || IsListCompleted();
}

bool ApplyIncrementalCPhiReductionProcess::IsBracketClosed() const
{
    return mBracketTolerance > 0.0 && mUpperBound.has_value() && mLastConvergedFactor.has_value() &&
           *mUpperBound - *mLastConvergedFactor <= mBracketTolerance + 1e-12;
}

double ApplyIncrementalCPhiReductionProcess::LimitIncrement(double StartFactor, double Increment) const
{
    const auto factor = StartFactor;
    auto       result = Increment;
    // approach the target factor with steps proportional to the remaining distance
    if (mTargetFactor > 0.0 && factor < mTargetFactor) {
        result = std::min(result, std::max(mTargetApproachFraction * (mTargetFactor - factor), mMinFactorIncrement));
    }
    // never reach a known upper bound: bisect the bracket instead
    if (mUpperBound.has_value()) {
        result = std::min(result, 0.5 * (*mUpperBound - factor));
    }
    return result;
}

bool ApplyIncrementalCPhiReductionProcess::IsListCompleted() const
{
    return mIncrementStrategy == "list" && mTrialNumber >= mFactorList.size();
}

double ApplyIncrementalCPhiReductionProcess::CalculateNewFactor() const
{
    KRATOS_ERROR_IF(mIncrementStrategy != "fixed" && mIncrementStrategy != "adaptive")
        << "Increment strategy '" << mIncrementStrategy << "' is not implemented." << std::endl;
    return mTrialFactor + LimitIncrement(mTrialFactor, mFactorIncrement);
}

double ApplyIncrementalCPhiReductionProcess::CalculateDisplacementMeasure() const
{
    const auto displacement_norm = [](const Node& rNode) {
        return norm_2(rNode.FastGetSolutionStepValue(DISPLACEMENT));
    };

    if (mDisplacementMeasure == "monitored_nodes") {
        const auto& r_root_model_part = mrModelParts.front().get().GetRootModelPart();
        auto        result            = 0.0;
        for (const auto node_id : mMonitoredNodeIds) {
            result = std::max(result, displacement_norm(r_root_model_part.GetNode(node_id)));
        }
        return result;
    }

    auto result = 0.0;
    for (const auto& r_model_part : mrModelParts) {
        const auto local_max = block_for_each<MaxReduction<double>>(
            r_model_part.get().GetCommunicator().LocalMesh().Nodes(), displacement_norm);
        result = std::max(
            result, r_model_part.get().GetCommunicator().GetDataCommunicator().MaxAll(local_max));
    }
    return result;
}

void ApplyIncrementalCPhiReductionProcess::UpdateSlope()
{
    mConvergedFactorsAndDisplacements.emplace_back(mTrialFactor, CalculateDisplacementMeasure());
    if (mConvergedFactorsAndDisplacements.size() <= mSlopeAveragingSteps) return;

    // Secant slope du/dF over the last mSlopeAveragingSteps converged steps
    const auto& [factor, displacement] = mConvergedFactorsAndDisplacements.back();
    const auto& [previous_factor, previous_displacement] =
        mConvergedFactorsAndDisplacements[mConvergedFactorsAndDisplacements.size() - 1 - mSlopeAveragingSteps];
    if (factor <= previous_factor) return;
    const auto slope = (displacement - previous_displacement) / (factor - previous_factor);

    if (!mReferenceSlope.has_value()) {
        if (slope > 0.0) {
            mReferenceSlope = slope;
            KRATOS_INFO("ApplyIncrementalCPhiReductionProcess")
                << "Reference slope du/dF = " << slope << std::endl;
        }
        return;
    }

    // k/k0 = slope_ref / slope, assigned to the middle of the secant interval
    if (slope > 0.0) {
        const auto stiffness_ratio = *mReferenceSlope / slope;
        mFactorsAndStiffnessRatios.emplace_back(0.5 * (factor + previous_factor), stiffness_ratio);
        KRATOS_INFO("ApplyIncrementalCPhiReductionProcess")
            << "Stiffness ratio k/k0 = " << stiffness_ratio << std::endl;

        if (mFailureCriterion == "stiffness_ratio" && stiffness_ratio < mStopStiffnessRatio) {
            mStoppedByStiffness = true;
            mExtrapolatedFactor = ExtrapolateFactorAtZeroStiffness();
            KRATOS_INFO("ApplyIncrementalCPhiReductionProcess")
                << "k/k0 = " << stiffness_ratio << " < " << mStopStiffnessRatio
                << ": failure reached, factor of safety (k -> 0) F = " << GetFactorOfSafety() << std::endl;
            return;
        }
    }

    if (mIncrementStrategy == "adaptive") {
        UpdateAdaptiveFactorIncrement(slope);
    }
}

std::optional<double> ApplyIncrementalCPhiReductionProcess::ExtrapolateFactorAtZeroStiffness() const
{
    // Least-squares line k/k0 = a + b F through the last points; F at k/k0 = 0 if k decreases
    if (mFactorsAndStiffnessRatios.size() < 2) return std::nullopt;
    const auto n     = std::min(mExtrapolationPoints, mFactorsAndStiffnessRatios.size());
    const auto first = mFactorsAndStiffnessRatios.end() - static_cast<std::ptrdiff_t>(n);
    auto mean_factor = 0.0, mean_ratio = 0.0;
    for (auto it = first; it != mFactorsAndStiffnessRatios.end(); ++it) {
        mean_factor += it->first / static_cast<double>(n);
        mean_ratio += it->second / static_cast<double>(n);
    }
    auto covariance = 0.0, variance = 0.0;
    for (auto it = first; it != mFactorsAndStiffnessRatios.end(); ++it) {
        covariance += (it->first - mean_factor) * (it->second - mean_ratio);
        variance += (it->first - mean_factor) * (it->first - mean_factor);
    }
    if (variance <= 0.0 || covariance >= 0.0) return std::nullopt;
    const auto extrapolated_factor = mean_factor - mean_ratio * variance / covariance;
    return std::max(extrapolated_factor, GetLastConvergedFactor());
}

void ApplyIncrementalCPhiReductionProcess::UpdateAdaptiveFactorIncrement(double Slope)
{
    // Same displacement increment per step as at the start: dF = dF_0 * slope_ref / slope
    const auto target_increment = Slope > 0.0
                                      ? mInitialFactorIncrement * (*mReferenceSlope / Slope)
                                      : mFactorIncrement * mMaxGrowthFactor;
    mFactorIncrement = std::clamp(std::min(target_increment, mFactorIncrement * mMaxGrowthFactor),
                                  mMinFactorIncrement, mMaxFactorIncrement);

    KRATOS_INFO("ApplyIncrementalCPhiReductionProcess")
        << "Slope du/dF = " << Slope << " (reference " << *mReferenceSlope
        << "), next factor increment = " << mFactorIncrement << std::endl;
}

void ApplyIncrementalCPhiReductionProcess::AdvanceFactor()
{
    if (mIncrementStrategy == "list") {
        if (!IsListCompleted()) mTrialFactor = mFactorList[mTrialNumber];
        return;
    }
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
