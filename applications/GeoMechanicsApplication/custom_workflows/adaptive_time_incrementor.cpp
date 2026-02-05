// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Wijtze Pieter Kikstra
//                   Anne van de Graaf
//

#include "adaptive_time_incrementor.h"
#include "time_step_end_state.hpp"

#include "includes/exception.h"

namespace Kratos
{

AdaptiveTimeIncrementor::AdaptiveTimeIncrementor(double                StartTime,
                                                 double                EndTime,
                                                 double                StartIncrement,
                                                 std::size_t           MaxNumOfCycles,
                                                 double                ReductionFactor,
                                                 double                IncreaseFactor,
                                                 std::optional<double> UserMinDeltaTime,
                                                 double                MaxTimeStepFactor,
                                                 std::size_t           MinNumOfIterations,
                                                 std::size_t           MaxNumOfIterations)
    : TimeIncrementor(),
      mEndTime(EndTime),
      mTimeSpan(EndTime - StartTime),
      mDeltaTime(std::min(StartIncrement, mTimeSpan)), // avoid exceeding the end time
      mInitialDeltaTime(mDeltaTime),
      mMaxNumOfCycles(MaxNumOfCycles),
      mReductionFactor(ReductionFactor),
      mIncreaseFactor(IncreaseFactor),
      mUserMinDeltaTime(std::move(UserMinDeltaTime)),
      mMaxDeltaTime(MaxTimeStepFactor * mDeltaTime),
      mMinNumOfIterations(MinNumOfIterations),
      mMaxNumOfIterations(MaxNumOfIterations)
{
    KRATOS_ERROR_IF(StartTime >= mEndTime)
        << "Start time (" << StartTime << ") must be smaller than end time (" << mEndTime << ")";
    KRATOS_ERROR_IF(mDeltaTime <= 0.0) << "Start increment must be positive, but got " << mDeltaTime;
    KRATOS_ERROR_IF(mMaxNumOfCycles == std::size_t(0))
        << "Maximum number of cycles must be positive";
    KRATOS_ERROR_IF(mMinNumOfIterations >= mMaxNumOfIterations)
        << "Minimum number of iterations (" << mMinNumOfIterations
        << ") is not less than maximum number of iterations (" << mMaxNumOfIterations << ")";
    KRATOS_ERROR_IF(mReductionFactor > 1.0)
        << "Reduction factor must not be greater than 1, but got " << mReductionFactor;
    KRATOS_ERROR_IF(mReductionFactor <= 0.0) << "Reduction factor must be positive, but got " << mReductionFactor;
    KRATOS_ERROR_IF(mIncreaseFactor < 1.0)
        << "Increase factor must be greater than or equal to 1, but got " << mIncreaseFactor;
    KRATOS_ERROR_IF(MaxTimeStepFactor < 1.0)
        << "Max_delta_time_factor must be greater than or equal to 1, but got " << MaxTimeStepFactor;
}

bool AdaptiveTimeIncrementor::WantNextStep(const TimeStepEndState& rPreviousState) const
{
    return rPreviousState.time < mEndTime;
}

bool AdaptiveTimeIncrementor::WantRetryStep(std::size_t CycleNumber, const TimeStepEndState& rPreviousState) const
{
    if (CycleNumber == 0) return true; // always carry out a first attempt

    if (rPreviousState.Converged()) return false; // the time step is done

    return CycleNumber < mMaxNumOfCycles; // stopping criterion
}

double AdaptiveTimeIncrementor::GetIncrement(double PreviousTime) const
{
    constexpr auto delta_time_as_fraction_of_time_span = 0.0001;
    const auto     default_min_delta_time =
        std::min(mInitialDeltaTime, delta_time_as_fraction_of_time_span * mTimeSpan);
    if (mEndTime - (PreviousTime + mDeltaTime) < mUserMinDeltaTime.value_or(default_min_delta_time)) {
        KRATOS_ERROR_IF(mEndTime - PreviousTime < mUserMinDeltaTime.value_or(default_min_delta_time))
            << "Delta time (" << mEndTime - PreviousTime << ") is smaller than "
            << (mUserMinDeltaTime ? "given" : "default") << " minimum allowable value "
            << mUserMinDeltaTime.value_or(default_min_delta_time) << std::endl;

        // Up-scaling to reach end_time without small increment
        return mEndTime - PreviousTime;
    }

    KRATOS_ERROR_IF(mDeltaTime < mUserMinDeltaTime.value_or(default_min_delta_time))
        << "Delta time (" << mDeltaTime << ") is smaller than "
        << (mUserMinDeltaTime ? "given" : "default") << " minimum allowable value "
        << mUserMinDeltaTime.value_or(default_min_delta_time) << std::endl;

    return mDeltaTime;
}

void AdaptiveTimeIncrementor::PostTimeStepExecution(const TimeStepEndState& rResultantState)
{
    constexpr auto delta_time_as_fraction_of_time_span = 0.0001;
    const auto     default_min_delta_time =
        std::min(mInitialDeltaTime, delta_time_as_fraction_of_time_span * mTimeSpan);

    if (rResultantState.Converged()) {
        if (rResultantState.num_of_iterations < mMinNumOfIterations) {
            mDeltaTime = std::min(mDeltaTime * mIncreaseFactor, mMaxDeltaTime);
        } else if (rResultantState.num_of_iterations == mMaxNumOfIterations) {
            mDeltaTime *= mReductionFactor;
        }
    }

    else {
        // non converged, scale down step and restart
        mDeltaTime *= mReductionFactor;
        KRATOS_ERROR_IF(mDeltaTime < mUserMinDeltaTime.value_or(default_min_delta_time))
            << "Delta time (" << mDeltaTime << ") is smaller than "
            << (mUserMinDeltaTime ? "given" : "default") << " minimum allowable value "
            << mUserMinDeltaTime.value_or(default_min_delta_time) << std::endl;
    }
}

} // namespace Kratos
