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

AdaptiveTimeIncrementor::AdaptiveTimeIncrementor(double      StartTime,
                                                 double      EndTime,
                                                 double      StartIncrement,
                                                 std::size_t MaxNumOfCycles,
                                                 double      ReductionFactor,
                                                 double      IncreaseFactor,
                                                 double      MaxTimeStepFactor,
                                                 std::size_t MinNumOfIterations,
                                                 std::size_t MaxNumOfIterations)
    : TimeIncrementor(),
      mEndTime(EndTime),
      mDeltaTime(std::min(StartIncrement, EndTime - StartTime)), // avoid exceeding the end time
      mMaxNumOfCycles(MaxNumOfCycles),
      mReductionFactor(ReductionFactor),
      mIncreaseFactor(IncreaseFactor),
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

double AdaptiveTimeIncrementor::GetIncrement() const { return mDeltaTime; }

void AdaptiveTimeIncrementor::PostTimeStepExecution(const TimeStepEndState& rResultantState)
{
    if (rResultantState.NonConverged() ||
        (rResultantState.Converged() && (rResultantState.num_of_iterations == mMaxNumOfIterations))) {
        mDeltaTime *= mReductionFactor;
    } else if (rResultantState.Converged() && (rResultantState.num_of_iterations < mMinNumOfIterations)) {
        mDeltaTime = std::min(mDeltaTime * mIncreaseFactor, mMaxDeltaTime);
    }

    // Avoid incrementing the time beyond the end time
    mDeltaTime = std::min(mDeltaTime, mEndTime - rResultantState.time);

    // Avoid very small remaining time steps
    const auto small_time_increment = 1.E-3 * mDeltaTime;
    if ((mEndTime - (rResultantState.time + mDeltaTime)) < small_time_increment) {
        mDeltaTime = mEndTime - rResultantState.time;
    }
}

} // namespace Kratos
