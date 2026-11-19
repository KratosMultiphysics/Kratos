// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//

#pragma once

#include "includes/element.h"
#include "includes/kratos_export_api.h"
#include "processes/process.h"

#include <cstddef>
#include <functional>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace Kratos
{


class Model;
class Parameters;

class KRATOS_API(GEO_MECHANICS_APPLICATION) ApplyIncrementalCPhiReductionProcess : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ApplyIncrementalCPhiReductionProcess);

    ApplyIncrementalCPhiReductionProcess(Model& rModel, const Parameters& rProcessSettings);
    void                      ExecuteInitializeSolutionStep() override;
    void                      ExecuteFinalizeSolutionStep() override;
    void                      ExecuteFinalize() override;
    int                       Check() override;
    [[nodiscard]] std::string Info() const override;
    [[nodiscard]] double      GetTrialFactor() const;
    [[nodiscard]] bool        HasLastConvergedFactor() const;
    [[nodiscard]] double      GetLastConvergedFactor() const;
    [[nodiscard]] double      GetFactorIncrement() const;
    [[nodiscard]] std::size_t GetTotalNumberOfCuts() const;
    [[nodiscard]] std::string GetStopReason() const;
    [[nodiscard]] double      GetStiffnessRatio() const;
    [[nodiscard]] double      GetFactorOfSafety() const;
    bool                      CutIncrementAfterNonConvergence();
    [[nodiscard]] std::size_t GetTrialNumber() const;
    [[nodiscard]] bool        IsFinished() const;

private:
    std::vector<std::reference_wrapper<ModelPart>> mrModelParts;

    std::string mIncrementStrategy;
    double mTrialFactor;
    std::optional<double> mLastConvergedFactor;
    double mInitialFactorIncrement;
    double mFactorIncrement;

    // "list" strategy: the factors to apply, in order
    std::vector<double> mFactorList;

    std::size_t mTrialNumber;
    std::size_t mMaxTrials;

    // On non-convergence the step is retried from the last converged factor with a smaller increment
    std::size_t mMaxCuts;
    double      mCutFactor;
    std::size_t mCutsInCurrentStep = 0;
    std::size_t mTotalNumberOfCuts = 0;
    bool        mStoppedByNonConvergence = false;

    // Bracketing: a failed trial whose increment was at most mBracketMaxIncrement is taken as an upper
    // bound of the factor of safety (larger failed increments may fail for numerical reasons only).
    // Later trials stay below that bound (bisection) and the reduction stops once the bound is within
    // mBracketTolerance of the last converged factor. Disabled when mBracketTolerance <= 0.
    double                mBracketTolerance;
    double                mBracketMaxIncrement;
    // A failed trial is only an upper bound when the slope is already softening (last k/k0 below this
    // value): failures while the slope is still stiff are failures of the iteration, not of the slope.
    double                mBracketMaxStiffnessRatio;
    std::optional<double> mUpperBound;
    bool                  mStoppedByBracket = false;

    // Target factor (e.g. a limit equilibrium estimate such as Bishop's): below it, the increment is
    // limited to mTargetApproachFraction times the remaining distance, so the steps shrink when approaching
    // the expected failure. Disabled when mTargetFactor <= 0.
    double mTargetFactor;
    double mTargetApproachFraction;

    // "stiffness_ratio" failure criterion: k = dF/du (inverse of the secant slope of the displacement
    // measure), stop when k/k0 < mStopStiffnessRatio and take as factor of safety the factor where the
    // linear fit of k/k0 versus F over the last mExtrapolationPoints points reaches zero
    std::string mFailureCriterion;
    double      mStopStiffnessRatio;
    std::size_t mExtrapolationPoints;
    std::vector<std::pair<double, double>> mFactorsAndStiffnessRatios;
    bool                  mStoppedByStiffness = false;
    std::optional<double> mExtrapolatedFactor;

    // "adaptive" strategy: the increment follows the slope of a displacement measure versus the
    // factor, so that every step produces roughly the same displacement increment as the first ones
    std::string              mDisplacementMeasure;
    std::vector<std::size_t> mMonitoredNodeIds;
    std::size_t              mSlopeAveragingSteps;
    double                   mMinFactorIncrement;
    double                   mMaxFactorIncrement;
    double                   mMaxGrowthFactor;
    std::vector<std::pair<double, double>> mConvergedFactorsAndDisplacements;
    std::optional<double>    mReferenceSlope;

    [[nodiscard]] bool IsListCompleted() const;

    [[nodiscard]] double CalculateNewFactor() const;

    [[nodiscard]] double LimitIncrement(double StartFactor, double Increment) const;

    [[nodiscard]] bool IsBracketClosed() const;

    [[nodiscard]] double CalculateDisplacementMeasure() const;

    void UpdateSlope();

    void UpdateAdaptiveFactorIncrement(double Slope);

    [[nodiscard]] std::optional<double> ExtrapolateFactorAtZeroStiffness() const;

    void ApplyReduction();

    void AdvanceFactor();

    [[nodiscard]] static double GetAndCheckCohesion(const Properties& rProperties);

    [[nodiscard]] static double GetAndCheckFrictionAngle(const Properties& rProperties,
                                                         IndexType         PropertyId);

    [[nodiscard]] static double CalculateReducedCohesion(
        double Cohesion,
        double Factor);

    [[nodiscard]] static double CalculateReducedFrictionAngle(
        double FrictionAngle,
        double Factor);

    static void SetCPhiAtElement(Element& rElement, double ReducedFrictionAngle, double ReducedCohesion);

    static void InitializeParametersForInternalMohrCoulombModel(Element& rElement);
};

} // namespace Kratos
